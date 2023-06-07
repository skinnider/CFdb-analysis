# Perform interactome mapping with transfer learning by borrowing labelled pairs
# from other datasets.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-transfer-learn-networks.R')
parser$add_argument('--allocation', type = 'character', 
                    default = "rrg-ljfoster-ab")
args = parser$parse_args()

library(tidyverse)
library(magrittr)

source("R/functions/detect_system.R")
source("R/functions/write_sh.R")
source("R/functions/submit_job.R")

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

# create output files
output_dir = file.path(base_dir, 'transfer_learning', 'networks')
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
output_filenames = pmap_chr(grid, function(...) {
  current = tibble(...)
  current %>%
    map2(., names(.), ~ {
      paste0(.y, '=', .x)
    }) %>%
    paste0(collapse = '-') %>% 
    paste0('.rds')
})
output_files = file.path(output_dir, output_filenames)

# now, check for which parameters are already complete
grid0 = grid %>% 
  mutate(output_file = output_files) %>%
  filter(!file.exists(output_file))

# write
grid_file = "sh/grids/transfer-learn-networks.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = TRUE)
write.table(grid0, grid_file, quote = FALSE, row.names = FALSE, sep = "\t")

# write the sh file dynamically
sh_file = '~/git/CFdb-analysis/sh/transfer-learn-networks.sh'
write_sh(job_name = 'transfer-learn-networks',
         sh_file = sh_file,
         grid_file = grid_file,
         inner_file = 'R/transfer_learning/inner-transfer-learn-networks.R',
         system = system,
         time = 24,
         cpus = 1,
         mem = 120)

# finally, run the job on whatever system we're on
submit_job(grid0, sh_file, args$allocation, system)
