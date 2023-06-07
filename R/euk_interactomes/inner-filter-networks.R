setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-filter-networks.R')
grid = read.delim("sh/grids/filter-networks.txt")
for (param_name in colnames(grid))
  parser$add_argument(paste0('--', param_name),
                      type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(PrInCE)
library(flavin)

# un-escape species from argparse
args$species %<>% chartr('_', ' ', .)

# read the network
ppi = readRDS(args$input_file)

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

# re-label and recalculate precision
labels = complexes %>%
  adjacency_matrix_from_list() %>%
  make_labels(ppi)
precision = calculate_precision(labels)
ppi$label = labels
ppi$precision = precision

# cut at nominal 50% precision
network = threshold_precision(ppi, threshold = 0.5)
if (is.null(network)) 
  network = head(ppi, 0)

# save the network
output_dir = dirname(args$output_file)
if (!dir.exists(output_dir)) 
  dir.create(output_dir, recursive = TRUE)
saveRDS(network, args$output_file)
