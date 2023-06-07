# Filter networks to 50% nominal precision.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-filter-transfer-networks.R')
parser$add_argument('--allocation', type = 'character', 
                    default = "rrg-ljfoster-ab")
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions.R")

# check that all of the dependent directories exist
if (!dir.exists("~/git/CFdb-searches"))
  stop("repository `CFdb-searches` does not exist")

# detect system
source("R/functions/detect_system.R")
source("R/functions/submit_job.R")
source("R/functions/write_sh.R")

# set up grid
grid = tidyr::crossing(
  prop_augmentation = c(0, 0.05, 0.1, 0.25, 0.33, 0.5)
) %>% 
  # other settings
  tidyr::crossing(
    classifier = 'RF',
    n_features = 1,
    feature_select = 'best_first',
    split_by = 'proteins',
    n_folds = 10,
    min_fractions = 4
  )

# create input files
input_dir = file.path(base_dir, 'transfer_learning', 'networks')
input_filenames = pmap_chr(grid, function(...) {
  current = tibble(...)
  current %>%
    map2(., names(.), ~ {
      paste0(.y, '=', .x)
    }) %>%
    paste0(collapse = '-') %>% 
    paste0('.rds')
})
input_files = file.path(input_dir, input_filenames)

# setup output files
output_dir = file.path(base_dir, 'transfer_learning', 'networks_50')
output_filenames = pmap_chr(grid, function(...) {
  current = tibble(...)
  current %>%
    dplyr::select(-(classifier:min_fractions)) %>% 
    map2(., names(.), ~ {
      paste0(.y, '=', .x)
    }) %>%
    paste0(collapse = '-') %>%
    paste0('network-', ., '.rds')
})
output_files = file.path(output_dir, output_filenames)

# subset grid
grid0 = grid %>% 
  mutate(input_file = input_files,
         output_file = output_files) %>% 
  filter(file.exists(input_file),
         !file.exists(output_file))

# write
grid_file = "sh/grids/filter-transfer-networks.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = TRUE)
write.table(grid0, grid_file, quote = FALSE, row.names = FALSE, sep = "\t")

# write the sh file dynamically
sh_file = '~/git/CFdb-analysis/sh/filter-transfer-networks.sh'
write_sh(job_name = 'filter-transfer-networks',
         sh_file = sh_file,
         grid_file = grid_file,
         inner_file = 'R/transfer_learning/inner-filter-transfer-networks.R',
         system = system,
         time = 8,
         mem = 64)

# finally, run the job on whatever system we're on
submit_job(grid0, sh_file, args$allocation, system)
