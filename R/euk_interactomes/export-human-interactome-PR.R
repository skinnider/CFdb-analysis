setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
args = list(); source("R/functions/detect_system.R")

# read all ~47m scored pairs
input_file = file.path(base_dir,
                       'euk_interactomes',
                       "network-species=Homo_sapiens-n_replicates=166-version=V2-n_features=1-classifier=RF-feature_select=best_first-split_by=proteins-n_folds=10-min_fractions=4.rds")
pairs = readRDS(input_file)

# extract ~same number of interactions as BioPlex 3
net = head(pairs, 125e3)

# save into git
saveRDS(net, "data/euk_interactomes/human-interactome-v2-PR.rds")