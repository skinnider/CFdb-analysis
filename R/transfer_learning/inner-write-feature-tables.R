setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-write-feature-tables.R')
grid = read.delim("sh/grids/write-feature-tables.txt")
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
library(nlme)
## other deps: preprocessCore; dismay

# read list of metrics from functions.R
source("R/functions.R")
source("R/functions/detect_system.R")
source("R/functions/get_features.R")

# un-escape species from argparse
args$Species %<>% chartr('_', ' ', .)

# define features to write for this dataset
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
                            'quantile', 'none')) %>% 
  dplyr::select(-transformation, -median)
top_n = head(best_features, 4)

# read chromatogram matrix
mat = readRDS(args$input_file)

# pre-filter any rows with <4 proteins
empty = rowSums(!is.na(mat) & is.finite(mat) & mat != 0) < 4
mat %<>% extract(!empty, )
# catch a single infinite value in one chromatogram
mat[is.infinite(mat)] = NA
# catch 'NaN' values in iBAQ/ratio
mat[is.nan(mat)] = NA

# filter to CORUM proteins
if (args$Species == 'Homo sapiens') {
  complexes = read.delim("data/complex/CORUM/complexes_human.txt") %>%
    as_annotation_list('gene_name', 'complex') %>% 
    extract(lengths(.) > 1)
} else {
  species_map = read.csv("data/GO/species-map.csv")
  proteome = filter(species_map, species == args$Species) %>% pull(proteome)
  complexes = paste0("data/euk_interactomes/CORUM/complexes_", 
                     proteome, ".txt") %>% 
    read.csv() %>%
    as_annotation_list('target', 'complex') %>% 
    extract(lengths(.) > 1)
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

# get features
matrix_dir = file.path(base_dir, "matrices", args$Accession, args$Replicate,
                       'default')
quant_mode = ifelse(args$Quantitation == 'TMT', 
                    'Reporter_intensity', 'iBAQ')
features = get_features(mat,
                        feature_grid = top_n,
                        n_features = 4,
                        matrix_dir = matrix_dir, 
                        quant_mode = quant_mode)
## don't impute NAs

# get labels
adj = adjacency_matrix_from_list(complexes)
labels = make_labels(adj, features)

# drop unlabelled proteins
features %<>%
  mutate(label = labels) %>% 
  drop_na(label)

# save feature table
output_dir = dirname(args$output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
saveRDS(features, args$output_file)
