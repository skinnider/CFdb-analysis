setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-filter-transfer-networks.R')
grid = read.delim("sh/grids/filter-transfer-networks.txt")
for (param_name in colnames(grid))
  parser$add_argument(paste0('--', param_name),
                      type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(PrInCE)
library(flavin)

# read the network
ppi = readRDS(args$input_file) %>% 
  # remove transferred features
  filter(!grepl('transfer_', protein_A),
         !grepl('transfer_', protein_B))

# read CORUM (mapped to bee with InParanoid)
complexes = read.csv("data/bee/complexes_bee.csv") %>% 
  as_annotation_list('gene_name', 'complex')

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
