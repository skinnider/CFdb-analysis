#' Score pairwise interactions based on cross-validation. 
predict_interactions = function(features, meta, splits,
                                classifier = c("NB", "RF", "SVM", "LR"),
                                multi_map = TRUE,
                                n_threads = 1) {
  library(naivebayes) 
  library(LiblineaR)
  library(ranger)
  library(randomForest)
  library(speedglm)
  
  # set up CV scores holder
  n_folds = length(splits)
  n_interactions = nrow(features)
  interaction_names = with(features, paste0(protein_A, '|', protein_B))
  clf_scores = matrix(NA, ncol = n_folds, nrow = n_interactions,
                      dimnames = list(interaction_names))
  
  # iterate through splits
  for (split_idx in seq_along(splits)) {
    message("working on split ", split_idx, " of ", length(splits), " ...")
    split = splits[[split_idx]]
    train_labels = split$train
    test_labels = split$test
    
    # make labels, accounting for multi-mapping
    if (multi_map) {
      labels = get_labels(features, train_labels, test_labels, meta)
      
      # catch datasets without positive examples
      if (n_distinct(na.omit(labels$train_label)) < 2) {
        out = data.frame()
        attr(out, 'error') = 'error'
        message("no labelled examples in split ", split_idx, "; exiting")
        return(out)
      }
    } else {
      train_labels = features %>%
        mutate_at(vars(protein_A, protein_B), ~ gsub(":.*$", "", .)) %>%
        left_join(train_labels, by = c('protein_A', 'protein_B')) %>% 
        pull(label)
      test_labels = features %>%
        mutate_at(vars(protein_A, protein_B), ~ gsub(":.*$", "", .)) %>%
        left_join(test_labels, by = c('protein_A', 'protein_B')) %>% 
        pull(label)
      labels = list(train_label = train_labels, test_label = test_labels)
    }
    
    # train classifier on labelled training data
    clf_idxs = which(!is.na(labels$train_label))
    clf_labels = labels$train_label[clf_idxs] %>%
      as.factor()
    clf_data = features[clf_idxs, ] %>%
      dplyr::select(-protein_A, -protein_B)
    clf_data_labeled = clf_data %>%
      mutate(label = clf_labels)
    clf = switch(classifier,
                 NB = naive_bayes(clf_data, clf_labels),
                 SVM = LiblineaR(clf_data, clf_labels, type = 2),
                 RF = ranger(data = clf_data_labeled, 
                             dependent.variable.name = "label",
                             probability = TRUE,
                             num.trees = 100,
                             num.threads = n_threads),
                 RF2 = randomForest(x = clf_data,
                                    y = clf_labels,
                                    ntree = 100),
                 LR = speedglm(label ~ ., clf_data_labeled, 
                               family = stats::binomial()))
    
    # predict interactions on the withheld set
    test_data = features[-clf_idxs, ] %>%
      dplyr::select(-protein_A, -protein_B)
    predictions = switch(
      classifier, 
      NB = predict(clf, test_data, type = 'prob', threshold = 5e-324),
      SVM = predict(clf, test_data, decisionValues = TRUE),
      RF = predict(clf, test_data, num.threads = n_threads),
      RF2 = predict(clf, test_data, type = 'prob'),
      LR = predict(clf, test_data, type = 'response'))
    
    ## extract predictions as numeric vector
    predictions = switch(
      classifier,
      NB = predictions[, "1"],
      SVM = -1.0 * predictions$decisionValues[, "0"],
      RF = predictions[[1]][, "1"],
      RF2 = predictions[, "1"],
      LR = predictions)
    # assign scores
    clf_scores[-clf_idxs, split_idx] = predictions
  }
  
  # average over folds 
  mean_scores = setNames(rowMeans(clf_scores, na.rm = TRUE),
                         rownames(clf_scores))
  
  # collate the labels from the last train/test split
  merged_labels = pmin(labels$train_label, labels$test_label, na.rm = TRUE)
  
  # create ranked data frame
  interactions = features %>%
    dplyr::select(protein_A, protein_B) %>%
    mutate(score = mean_scores, label = merged_labels) %>%
    arrange(desc(score)) %>%
    # calculate precision
    mutate(precision = calculate_precision(label))
  return(interactions)
}

get_labels = function(features, train_labels, test_labels, meta) {
  # first, map protein groups to gene IDs, accounting for multi-mapping
  gene_map = meta %>% dplyr::select(`Gene names`, `Majority protein IDs`)
  ## we can save some time by prefiltering to labelled genes
  labelled_genes = c(
    with(train_labels, c(protein_A, protein_B)),
    with(test_labels, c(protein_A, protein_B))
  ) %>% unique()
  gene_map %<>%
    mutate(`Gene names` = strsplit(`Gene names`, ';')) %>%
    unnest(`Gene names`) %>%
    filter(`Gene names` %in% labelled_genes) %>%
    group_by(`Majority protein IDs`) %>%
    summarise(`Gene names` = paste0(`Gene names`, collapse = ';')) %>%
    ungroup()
  gene_features = features %>%
    dplyr::select(protein_A, protein_B) %>%
    dplyr::rename(`Majority protein IDs` = protein_A) %>%
    left_join(gene_map, by = 'Majority protein IDs') %>%
    dplyr::rename(protein_A = `Majority protein IDs`,
                  gene_A = `Gene names`,
                  `Majority protein IDs` = protein_B) %>%
    left_join(gene_map, by = 'Majority protein IDs') %>%
    dplyr::rename(protein_B = `Majority protein IDs`,
                  gene_B = `Gene names`) %>%
    # flag original row and filter to pairs in the reference
    mutate(idx = row_number()) %>%
    filter(!is.na(gene_A), !is.na(gene_B)) %>%
    # split gene names
    mutate(gene_A = strsplit(gene_A, ';'),
           gene_B = strsplit(gene_B, ';')) %>%
    unnest(gene_A) %>% unnest(gene_B)
  
  # now, label gene pairs
  train_df = train_labels %>%
    dplyr::rename(gene_A = protein_A, gene_B = protein_B, train_label = label)
  test_df = test_labels %>%
    dplyr::rename(gene_A = protein_A, gene_B = protein_B, test_label = label)
  gene_labels = gene_features %>%
    left_join(train_df, by = c('gene_A', 'gene_B')) %>%
    left_join(test_df, by = c('gene_A', 'gene_B'))
  
  # consolidate to protein IDs
  labels = features %>%
    mutate(idx = row_number()) %>%
    dplyr::select(protein_A, protein_B, idx) %>%
    left_join(gene_labels, by = c('protein_A', 'protein_B', 'idx')) %>%
    group_by(idx) %>%
    # ignore pairs with >1 potential label
    summarise(train_label = ifelse(n_distinct(na.omit(train_label)) > 1, 
                                   NA, dplyr::first(train_label)),
              test_label = ifelse(n_distinct(na.omit(test_label)) > 1, 
                                  NA, dplyr::first(test_label))) %>%
    ungroup()
  return(labels)
}
