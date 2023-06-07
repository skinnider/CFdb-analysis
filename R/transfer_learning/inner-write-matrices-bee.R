setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-write-matrices-bee.R')
grid = read.delim("sh/grids/write_matrices_bee.txt")
for (param_name in colnames(grid))
  parser$add_argument(paste0('--', param_name),
                      type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

# read list of metrics from functions.R
source("R/functions.R")

# set up output filepath
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = T)

library(tidyverse)
library(magrittr)
library(PrInCE)
library(flavin)
library(AUC)
library(ontologyIndex)
library(dismay)
library(nlme)
## other deps: preprocessCore; dismay

# fix boolean
args$scale %<>% as.logical()

# read chromatogram matrix
mat = readRDS(args$input_file)[[args$replicate]]
# pre-filter any empty rows
empty = rowSums(!is.na(mat) & is.finite(mat) & mat != 0) == 0
mat %<>% extract(!empty, )

# catch a single infinite value in one chromatogram
mat[is.infinite(mat)] = NA
# catch 'NaN' values in iBAQ/ratio
mat[is.nan(mat)] = NA

# step 1: impute missing values (optional)
if (args$missing == 'zero') {
  mat[is.na(mat)] = 0
} else if (args$missing == 'noise') {
  mat[mat == 0] = NA
  mat %<>% clean_profiles(impute_NA = TRUE, smooth = FALSE)
} else if (args$missing == 'NA') {
  mat[mat == 0] = NA
}

# step 2: log-transform (optional)
if (args$transform == "log") {
  mat %<>% log()
  if (args$missing == 'zero') {
    mat[mat == -Inf] = 0
  }
}

# step 3: scale to range [0, 1] (optional)
if (args$scale) {
  scale01 = function(vec) (vec - min(vec, na.rm = T)) / 
    (max(vec, na.rm = T) - min(vec, na.rm = T))
  mat %<>% apply(1, scale01) %>% t()
  ## NaNs are produced from vectors that have only a single observation:
  ## convert these to 1s
  mat %<>% replace(is.nan(.), 1)
}

# step 4: quantile normalize (optional)
if (args$normalize == 'quantile') {
  dims = dimnames(mat)
  mat %<>% preprocessCore::normalize.quantiles()
  dimnames(mat) = dims
}

# step 5: calculate pairwise interactions
cor = get_coexpr(t(mat), args$metric)

# save matrix
saveRDS(cor, args$output_file)
