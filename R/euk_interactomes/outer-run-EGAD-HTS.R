# Run EGAD on other large-scale protein interaction networks derived from 
# high-throughput screens
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-run-EGAD-HTS.R')
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

# get human networks
hs_files = list.files("data/high_throughput", full.names = TRUE, 
                      pattern = '\\.gz$')
# get mouse networks
mm_files = list.files("data/resources/interactomes", pattern = "Skinnider2021",
                      full.names = TRUE)

# construct list of experiments
grid = data.frame(input_file = c(hts_files1, hts_files2)) %>% 
  mutate(network_name = gsub("\\..*$", "", basename(input_file)) %>% 
           gsub("-.*$", "", .),
         species = ifelse(grepl('Skinnider2021', network_name),
                          'Homo sapiens', 'Mus musculus')) %>% 
  # escape spaces for argparse
  mutate(species = chartr(' ', '_', species))

# setup output files
output_dir = file.path(base_dir, 'euk_interactomes', 'EGAD')
output_filenames = basename(grid$input_file) %>% 
  gsub("\\..*$", "", .) %>% 
  paste0('.rds')
output_files = file.path(output_dir, output_filenames)

# subset grid
grid0 = grid %>% 
  mutate(output_file = output_files) %>% 
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
         mem = 120)

# finally, run the job on whatever system we're on
submit_job(grid0, sh_file, args$allocation, system)
