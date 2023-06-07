split_by_proteins = function(reference, percent_in_train = 0.7) {
  # convert to adjacency matrix
  adj = reference %>%
    distinct(protein_A, protein_B) %>%
    PrInCE::adjacency_matrix_from_data_frame()
  # shuffle the order of the rows and columns
  shuff = sample(rownames(adj))
  adj %<>% extract(shuff, shuff)
  # now, split at 70/30 border
  border = floor(nrow(adj) * percent_in_train)
  train_proteins = rownames(adj)[seq_len(border)]
  test_proteins = setdiff(rownames(adj), train_proteins)
  # filter the original reference
  train = reference %>%
    filter(protein_A %in% train_proteins, protein_B %in% train_proteins)
  test = reference %>%
    filter(protein_A %in% test_proteins, protein_B %in% test_proteins)
  return(list('train' = train, 'test' = test))
}

cv_by_proteins = function(reference, n_folds = 5) {
  # convert to adjacency matrix
  adj = reference %>%
    distinct(protein_A, protein_B) %>%
    PrInCE::adjacency_matrix_from_data_frame()
  
  # split the matrix into folds
  borders = seq(1, nrow(adj), nrow(adj) / n_folds) %>%
    floor()
  fold_idxs = seq_len(nrow(adj)) %>%
    sample() %>%
    split(borders)
  folds = map(fold_idxs, ~ extract(adj, ., .)) %>%
    setNames(seq_along(.))
  
  # convert back to long
  folds %<>%
    map(~ reshape2::melt(., varnames = c('protein_A', 'protein_B'), 
                         value.name = 'label', as.is = TRUE) %>% 
          filter(protein_A < protein_B))
  
  return(folds)
}

cv_by_pairs = function(reference, n_folds = 5) {
  # convert to adjacency matrix to label true-negatives
  adj = reference %>%
    distinct(protein_A, protein_B) %>%
    adjacency_matrix_from_data_frame()
  
  # convert back to long
  pairs = adj %>%
    reshape2::melt(varnames = c('protein_A', 'protein_B'), 
                   value.name = 'label', as.is = TRUE) %>% 
    filter(protein_A < protein_B)
  
  # now, split the pairs into folds
  borders = seq(1, nrow(pairs), nrow(pairs) / n_folds) %>%
    floor()
  fold_idxs = seq_len(nrow(pairs)) %>%
    sample() %>%
    split(borders)
  folds = map(fold_idxs, ~ extract(pairs, ., )) %>%
    setNames(seq_along(.))
  return(folds)
}

get_cv_labels = function(folds, mode = 'pairs',
                         remove_overlapping_labels = TRUE) {
  splits = list()
  if (mode == "pairs") {
    for (test_idx in seq_along(folds)) {
      test = folds[[test_idx]]
      train = folds[-test_idx] %>% bind_rows()
      # add to split
      split = list(train = train, test = test)
      splits[[test_idx]] = split
    }
  }
  splits %<>% setNames(seq_along(.))
  return(splits)
}
