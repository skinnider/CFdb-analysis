setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-run-EGAD.R')
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

# read a list of species
expts = read.csv("~/git/CFdb-searches/data/experiments.csv")
# count the number of replicates for each species, in each version
n_repl1 = filter(expts, Version == 'V1') %>% 
  dplyr::count(Species, name = 'n_replicates') %>% 
  dplyr::rename(species = Species)
n_repl2 = expts %>% 
  dplyr::count(Species, name = 'n_replicates') %>% 
  dplyr::rename(species = Species)
n_repl = bind_rows(
  n_repl1 %>% mutate(version = 'V1'),
  n_repl2 %>% mutate(version = 'V2')
) %>% 
  # remove V2 that has not grown since V1
  group_by(species, n_replicates) %>% 
  arrange(version) %>% 
  dplyr::slice(1) %>% 
  ungroup() %>% 
  # remove bacteria
  filter(!species %in% c("Anabaena sp. 7120",
                         "Cyanothece ATCC 51142",
                         "Escherichia coli",
                         "Kuenenia stuttgartiensis",
                         "Salmonella typhimurium SL1344",
                         "Synechocystis sp. PCC 6803")) %>% 
  # escape spaces for argparse
  mutate(species = chartr(' ', '_', species))

# establish grid of analyses
grid = n_repl %>% 
  # the number of features depends on the number of replicates 
  mutate(n_features = ifelse(n_replicates <= 5, 4, 
                             ifelse(n_replicates <= 10, 2, 1))) %>% 
  # dplyr::select(-n_replicates) %>% 
  tidyr::crossing(
    classifier = 'RF',
    feature_select = 'best_first',
    split_by = 'proteins',
    n_folds = 10,
    min_fractions = 4
  ) %>% 
  # filter to 3+ replicates
  filter(n_replicates >= 3)

# set up input files
input_dir = file.path(base_dir, 'euk_interactomes', 'networks')
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
output_dir = file.path(base_dir, 'euk_interactomes', 'EGAD')
output_files = file.path(output_dir, input_filenames)

# subset grid
grid0 = grid %>% 
  mutate(input_file = input_files,
         output_file = output_files) %>% 
  filter(file.exists(input_file),
         !file.exists(output_file))

# write
grid_file = "sh/grids/run-EGAD.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = TRUE)
write.table(grid0, grid_file, quote = FALSE, row.names = FALSE, sep = "\t")

# write the sh file dynamically
sh_file = '~/git/CFdb-analysis/sh/run-EGAD.sh'
write_sh(job_name = 'run-EGAD',
         sh_file = sh_file,
         grid_file = grid_file,
         inner_file = 'R/euk_interactomes/inner-run-EGAD.R',
         system = system,
         time = 8,
         mem = 72)

# finally, run the job on whatever system we're on
submit_job(grid0, sh_file, args$allocation, system)
