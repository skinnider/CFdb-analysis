# Write feature tables for labelled protein pairs in human and mouse CF-MS 
# experiments.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-write-feature-tables.R')
parser$add_argument('--allocation', type = 'character', 
                    default = "rrg-ljfoster-ab")
args = parser$parse_args()

library(tidyverse)
library(magrittr)

# check that all of the dependent directories exist
if (!dir.exists("~/git/CFdb-searches"))
  stop("repository `CFdb-searches` does not exist")

# what system are we on?
source("R/functions.R") ## contains metrics used to predict interactions
source("R/functions/detect_system.R")
source("R/functions/write_sh.R")
source("R/functions/submit_job.R")

# read experiments
expts = read.csv("~/git/CFdb-searches/data/experiments.csv") %>% 
  mutate(Replicate = fct_recode(Replicate,
                                'BN3to12_RAW' = 'BN3to12',
                                'BN4to16_RAW' = 'BN4to16',
                                'Ghosts_DDM_BiobasicSEC-2' = 'Ghosts_6061_DDM_BiobasicSEC-2',
                                'Hemolysate_old_9330_SEC' = 'Hemolysate_9330_old_SEC',
                                'whole_rbc_IP_lysis_6840' = 'whole_rbc_IP_lysis',
                                'Hemolysate_8994_IEX-3' = 'Hemolysate_8994_IEX') %>% 
           as.character())
# filter to human and mouse 
expts0 = filter(expts, Species %in% c('Homo sapiens', 'Mus musculus')) %>% 
  # escape spaces for argparse
  mutate(Species = chartr(' ', '_', Species)) %>% 
  # remove fractionation column
  dplyr::select(-Fractionation) %>% 
  # remove empty quantitations
  mutate(Quantitation = ifelse(Quantitation == 'TMT', Quantitation, 'other'))

# create input files
chrom_dirs = with(expts0, file.path("~/git/CFdb-searches/data/chromatograms_gene",
                                    Accession, Replicate, "default"))
chrom_files = paste0(ifelse(expts0$Quantitation == 'TMT', 
                            'Reporter_intensity', 'iBAQ'), '.rds') %>% 
  file.path(chrom_dirs, .)

# create output files
output_dir = file.path(base_dir, 'transfer_learning', 'feature_tables')
output_dirnames = with(expts0, file.path(Accession, Replicate, "default"))
output_files = file.path(output_dir, output_dirnames, 'features-top-4.rds')

# now, check for which parameters are already complete
grid0 = expts0 %>% 
  mutate(input_file = chrom_files,
         output_file = output_files) %>%
  filter(!file.exists(output_file))

# write
grid_file = "sh/grids/write-feature-tables.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = TRUE)
write.table(grid0, grid_file, quote = FALSE, row.names = FALSE, sep = "\t")

# write the sh file dynamically
sh_file = '~/git/CFdb-analysis/sh/write-feature-tables.sh'
write_sh(job_name = 'write-feature-tables',
         sh_file = sh_file,
         grid_file = grid_file,
         inner_file = 'R/transfer_learning/inner-write-feature-tables.R',
         system = system,
         time = 24,
         cpus = 1,
         mem = 96)

# finally, run the job on whatever system we're on
submit_job(grid0, sh_file, args$allocation, system)
