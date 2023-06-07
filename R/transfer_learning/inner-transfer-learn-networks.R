setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-transfer-learn-networks.R')
grid = read.delim("sh/grids/transfer-learn-networks.txt")
for (param_name in colnames(grid))
  parser$add_argument(paste0('--', param_name),
                      type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(PrInCE)
library(flavin)
library(AUC)
library(ontologyIndex)
library(dismay)

source("R/functions.R")
source("R/functions/split_by_complex.R")
source("R/functions/get_features.R")
source("R/functions/predict_interactions.R")
source("R/functions/detect_system.R")

# read chromatograms
chrom = readRDS("data/bee/chromatograms.rds")

# read CORUM (mapped to bee with InParanoid)
complexes = read.csv("data/bee/complexes_bee.csv") %>% 
  as_annotation_list('gene_name', 'complex')
complex_proteins = unique(unlist(complexes))
n_subunits = lengths(complexes)

# create pairwise interactions 
reference = complexes %>% 
  map_dfr(~ tidyr::crossing(protein_A = ., protein_B = .), .id = 'complex') %>%
  filter(protein_A < protein_B) %>%
  distinct() %>%
  # for interactions found in more than one complex, pick the smaller one
  mutate(n_subunits = n_subunits[complex] %>% unname()) %>%
  group_by(protein_A, protein_B) %>%
  arrange(n_subunits) %>%
  mutate(keep = row_number() == 1) %>%
  ungroup() %>%
  filter(keep)

# set up cross-validation folds in the reference set
percent_in_train = (args$n_folds - 1) / args$n_folds
if (args$split_by == 'proteins') {
  folds = cv_by_proteins(reference, n_folds = args$n_folds)
  splits = get_cv_labels(folds, mode = 'pairs')
} else if (args$split_by == 'pairs') {
  folds = cv_by_pairs(reference, n_folds = args$n_folds)
  splits = get_cv_labels(folds, mode = 'pairs')
} else {
  stop("not sure how to split by: ", args$split_by)
}

# filter chromatogram matrices
mats = map(chrom, ~ {
  mat = .x
  # filter by n_fractions
  keep = rowSums(!is.na(mat) & is.finite(mat) & mat != 0) >= args$min_fractions
  mat %<>% extract(keep, )
  return(mat)
})

# calculate features on each matrix
feats = list()
for (mat_idx in seq_along(mats)) {
  message("calculating features for matrix ", mat_idx, " of ", length(mats), 
          " ...")
  mat = mats[[mat_idx]]
  
  # read pre-calculated feature matrices, if we can
  matrix_dir = file.path(base_dir, "matrices", names(chrom)[mat_idx])
  quant_mode = 'matrix'
  
  # calculate features
  if (args$feature_select == "best_first") {
    best_features = readRDS("data/analysis_grid/best_features.rds") %>%
      extract2('combined') %>%
      # keep only one metric
      group_by(metric) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      arrange(desc(median)) %>%
      # spread 'transformation' back out
      mutate(transform = ifelse(transformation == 'log-transform',
                                'log', 'none'),
             normalize = ifelse(transformation == 'quantile\nnormalization',
                                'quantile', 'none'))
    grid = head(best_features, args$n_features)
    features = get_features(mat,
                            feature_grid = grid,
                            n_features = args$n_features, 
                            matrix_dir = matrix_dir, 
                            quant_mode = quant_mode)
  } else if (args$feature_select == "random") {
    features = get_features(mat,
                            feature_grid = NULL,
                            n_features = args$n_features, 
                            matrix_dir = matrix_dir,
                            quant_mode = quant_mode)
  } else if (args$feature_select == "PrInCE") {
    features = PrInCE::calculate_features(mat, gaussians = NULL, 
                                          co_apex = FALSE)
  }
  
  # filter to proteins that passed the filter
  features %<>% filter(protein_A %in% rownames(mat),
                       protein_B %in% rownames(mat))
  
  # label features
  adj = adjacency_matrix_from_list(complexes)
  labels = make_labels(adj, features)
  features %<>% mutate(label = labels)
  
  # add to list
  feats[[mat_idx]] = features
}

# now, add a fixed proportion of labelled features to each replicate
if (args$prop_augmentation > 0) {
  feature_table_file = file.path(base_dir, 'transfer_learning', 
                                 'features-table.rds')
  feature_table = readRDS(feature_table_file)
  set.seed(0)
  for (mat_idx in seq_along(mats)) {
    features = feats[[mat_idx]]
    
    # count number of labelled pairs
    n_labels = sum(!is.na(features$label))
    # how many do we want to add?
    n_to_add = round(n_labels * args$prop_augmentation)
    n_pos = round(n_to_add / 2)
    n_neg = round(n_to_add / 2)
    
    # add labelled pairs to the dataset
    feature_table_pos = feature_table %>% 
      filter(label == 1) %>% 
      sample_n(n_pos)
    feature_table_neg = feature_table %>% 
      filter(label == 0) %>% 
      sample_n(n_neg)
    feature_table0 = bind_rows(feature_table_pos, feature_table_neg)
    
    # keep only same features as the bee matrix
    feature_table0 %<>% extract(, colnames(features))
    
    # assign names to transfer features based on labelled data in this replicate
    # (this way, we merge transferred proteins across replicates the same way as
    # real proteins)
    labelled_features = filter(features, !is.na(label)) %>% 
      # only keep as many as there are rows in feature_table0
      sample_n(nrow(feature_table0)) %>% 
      # tag these as 'transferred' proteins
      mutate(protein_A = paste0('transfer_', protein_A),
             protein_B = paste0('transfer_', protein_B))
    # now, replace proteins in feature_table0 with these
    feature_table0 %<>%
      dplyr::select(-protein_A, -protein_B) %>% 
      mutate(protein_A = labelled_features$protein_A,
             protein_B = labelled_features$protein_B) %>% 
      dplyr::select(protein_A, protein_B, everything())
    
    # add to existing features
    features %<>% bind_rows(feature_table0)
    feats[[mat_idx]] = features
  }
  ## clear table from memory
  rm(feature_table)
}

# combine features from each matrix
message("concatenating ", length(feats), " features ...")
features = feats %>% 
  # remove labels
  map(~ dplyr::select(.x, -label)) %>% 
  concatenate_features()
message("  done concatenating features")
## clear matrices and feature matrices from memory
rm(mats)
rm(feats)

# replace missing data
message("replacing missing data ...")
features %<>% replace_missing_data()

# predict interactions
message("predicting interactions ...")
network = predict_interactions(features, meta, splits,
                               classifier = args$classifier,
                               multi_map = FALSE)

# save the network itself
output_dir = dirname(args$output_file)
if (!dir.exists(output_dir)) 
  dir.create(output_dir, recursive = TRUE)
saveRDS(network, args$output_file)
