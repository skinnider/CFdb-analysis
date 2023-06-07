# Predict gene-level consensus CF-MS interactomes for each species in the
# resource.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-species-interactomes.R')
grid = read.delim("sh/grids/species-interactomes.txt")
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

# un-escape species from argparse
args$species %<>% chartr('_', ' ', .)

# read experiments from this species
expts = read.csv("~/git/CFdb-searches/data/experiments.csv") %>%
  filter(Species == args$species) %>% 
  mutate(Replicate = fct_recode(Replicate,
                                'BN3to12_RAW' = 'BN3to12',
                                'BN4to16_RAW' = 'BN4to16',
                                'Ghosts_DDM_BiobasicSEC-2' = 'Ghosts_6061_DDM_BiobasicSEC-2',
                                'Hemolysate_old_9330_SEC' = 'Hemolysate_9330_old_SEC',
                                'whole_rbc_IP_lysis_6840' = 'whole_rbc_IP_lysis',
                                'Hemolysate_8994_IEX-3' = 'Hemolysate_8994_IEX'))
# filter by version
if (args$version == 'V1')
  expts %<>% filter(Version == 'V1')

# read CORUM
if (args$species == 'Homo sapiens') {
  complexes = read.delim("data/complex/CORUM/complexes_human.txt") %>%
    as_annotation_list('gene_name', 'complex')
} else {
  species_map = read.csv("data/GO/species-map.csv")
  proteome = filter(species_map, species == args$species) %>% pull(proteome)
  complexes = paste0("data/euk_interactomes/CORUM/complexes_", 
                     proteome, ".txt") %>% 
    read.csv() %>%
    as_annotation_list('target', 'complex')
}
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

# now, read the gene-level input datasets
input_dirs = file.path("~/git/CFdb-searches/data/chromatograms_gene", 
                       expts$Accession, 
                       expts$Replicate, 
                       'default')
input_files = file.path(input_dirs, 'iBAQ.rds')
input_files[!file.exists(input_files)] %<>% gsub("iBAQ", 'Reporter_intensity', .)

# read chromatogram matrices, already mapped to genes
mats = map(input_files, ~ {
  mat = readRDS(.)
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
  matrix_dir = file.path(base_dir, "matrices")
  chrom_dir = input_dirs[mat_idx] %>%
    gsub("^.*chromatograms_gene\\/", "", .)
  matrix_dir %<>% file.path(chrom_dir)
  
  # quant mode: iBAQ or reporter intensity
  matrix_files = list.files(matrix_dir, pattern = 'rds')
  quant_mode = ifelse(!any(grepl('^iBAQ', matrix_files)), 
                      'Reporter_intensity', 'iBAQ')
  
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
  feats[[mat_idx]] = features
}

# combine features from each matrix
message("concatenating ", length(feats), " features ...")
features = concatenate_features(feats)
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
