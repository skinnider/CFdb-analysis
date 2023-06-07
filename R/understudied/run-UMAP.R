setwd("~/git/CFdb-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(Matrix)
library(uwot)
args = list(); source("R/functions/detect_system.R")

# load input matrix
input_file = file.path(base_dir, 'UMAP/features-human.rds')
mat = readRDS(input_file)
# impute missing values with their median
mat %<>% apply(1, function(x) replace(x, is.na(x), median(x, na.rm = TRUE)))

# run UMAP on the integrated data
emb = uwot::umap(X = mat, n_neighbors = 100, spread = 2, min_dist = 0.05)

# save results
output_file = file.path(base_dir, 'UMAP/coords-human.rds')
saveRDS(emb, output_file)
