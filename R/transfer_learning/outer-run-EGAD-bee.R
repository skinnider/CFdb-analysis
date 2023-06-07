setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-run-EGAD-bee.R')
parser$add_argument('--allocation', type = 'character', 
                    default = "st-ljfoster-1")
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions.R")

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

# setup input files
input_dir = file.path(base_dir, 'transfer_learning', 'networks_50')
input_filenames = pmap_chr(grid, function(...) {
  current = tibble(...)
  current %>%
    dplyr::select(-(classifier:min_fractions)) %>% 
    map2(., names(.), ~ {
      paste0(.y, '=', .x)
    }) %>%
    paste0(collapse = '-') %>%
    paste0('network-', ., '.rds')
})
input_files = file.path(input_dir, input_filenames)

# setup output files
output_dir = file.path(base_dir, 'transfer_learning', 'EGAD')
output_files = file.path(output_dir, input_filenames)

# subset grid
grid0 = grid %>% 
  mutate(input_file = input_files,
         output_file = output_files) %>% 
  filter(file.exists(input_file),
         !file.exists(output_file))

# write
grid_file = "sh/grids/run-EGAD-bee.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = TRUE)
write.table(grid0, grid_file, quote = FALSE, row.names = FALSE, sep = "\t")

# write the sh file dynamically
sh_file = '~/git/CFdb-analysis/sh/run-EGAD-bee.sh'
write_sh(job_name = 'run-EGAD-bee',
         sh_file = sh_file,
         grid_file = grid_file,
         inner_file = 'R/transfer_learning/inner-run-EGAD-bee.R',
         system = system,
         time = 8,
         mem = 72)

# finally, run the job on whatever system we're on
submit_job(grid0, sh_file, args$allocation, system)
