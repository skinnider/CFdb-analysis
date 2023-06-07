# Write a more complete E. coli PR curve, out to ~13,000 interactions,
# to match past HTS networks.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(PrInCE)
library(flavin)
args = list(); source("R/functions/detect_system.R")

# read the network
ppi = file.path(base_dir, "prok_interactomes",
                "network-species=Escherichia_coli-n_replicates=5-version=V2-n_features=4-classifier=RF-feature_select=best_first-split_by=proteins-n_folds=10-min_fractions=4.rds") %>% 
  readRDS()

# match # of PPIs to HTS
## Babu2017 = 12,801
pr = head(ppi, 15e3)

# calculate precision
complexes = read.csv("data/prok_interactomes/EcoCyc/complexes_UP000000625.csv") %>% 
  as_annotation_list('gene_name', 'complex')
labels = complexes %>%
  adjacency_matrix_from_list() %>%
  make_labels(pr)
precision = calculate_precision(labels)
pr$label = labels
pr$precision = precision

# save the network
output_file = "data/prok_interactomes/e-coli-PR.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir)) 
  dir.create(output_dir, recursive = TRUE)
saveRDS(pr, output_file)
