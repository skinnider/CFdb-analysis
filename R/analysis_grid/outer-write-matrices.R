# For a subset of metrics that take a very long time to compute, calculate
# correlation matrices and write them to a file.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-write-matrices.R')
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

# establish grid of analyses
opts = list(
  analysis = 'complexes',
  metric = metrics,
  scale = FALSE,
  transform = c('none', 'log'),
  min_pairs = 0,
  normalize = c('quantile', 'none'),
  missing = c('NA', 'zero', 'noise')
)
grid = do.call(expand.grid, c(opts, stringsAsFactors = F)) %>%
  # pick only one transformation
  filter(
    # log-transform
    (scale == FALSE & transform == 'log' & normalize == 'none') |
      # quantile normalize
      (scale == FALSE & transform == 'none' & normalize == 'quantile') |
      # scale
      (scale == TRUE & transform == 'none' & normalize == 'none') |
      # ... or none of them
      (scale == FALSE & transform == 'none' & normalize == 'none')) %>%
  # filter some combinations that cannot deal with NAs
  filter(!(missing == 'NA' & 
             metric %in% c('phi_s', 'rho_p', 'partial_cor', 'GENIE3',
                           'distance_cor', 'cosine', 'wccor'))) %>%
  # filter some combinations for which zeros and NAs are identical
  filter(!(missing == 'NA' &
             metric %in% c('dice', 'hamming', 'jaccard', 'zi_kendall',
                           'binomial', 'bayes_cor'))) 

## keep only the best feature for now
best_features = readRDS("data/analysis_grid/best_features.rds") %>%
  extract2('combined') %>%
  # keep only one metric
  group_by(metric) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  arrange(desc(median)) %>%
  # spread 'transformation' back out
  mutate(transform = ifelse(transformation == 'log-transform',
                            'log', 'none'),
         normalize = ifelse(transformation == 'quantile\nnormalization',
                            'quantile', 'none')) %>% 
  dplyr::select(-transformation, -median)
top_n = head(best_features, 4)
## filter grid accordingly
grid %<>% inner_join(top_n)

# combine with quantitation strategies (input files)
chrom_dir = '~/git/CFdb-searches/data/chromatograms_gene'
chrom_files = list.files(chrom_dir, pattern = '*.rds', recursive = T) %>%
  # ignore the metadata files
  extract(!grepl("metadata", .))
split = strsplit(chrom_files, '/')
accessions = map_chr(split, 1)
experiments = map_chr(split, 2)
searches = map_chr(split, 3)
quant_modes = gsub("\\..*$", "", basename(chrom_files))
inputs = data.frame(file = file.path(chrom_dir, chrom_files),
                    accession = accessions,
                    experiment = experiments,
                    search = searches,
                    quant_mode = quant_modes) %>%
  # limit to default search only
  filter(search == 'default') %>%
  # dplyr::select(-search) %>%
  # process quant modes in 'inner' script
  distinct(accession, experiment, search, quant_mode) %>%
  # merge these
  unite(input, accession, experiment, search, quant_mode, sep = '|')

# rep each analysis over each input
grid %<>%
  dplyr::slice(rep(1:n(), each = nrow(inputs))) %>%
  mutate(input = rep(inputs$input, nrow(grid))) %>%
  left_join(inputs, by = 'input') %>%
  separate(input, into = c("accession", "experiment", "search", "quant_mode"),
           sep = "\\|")

# clean up grid
grid %<>%
  dplyr::select(accession, experiment, search, quant_mode, analysis, 
                metric, scale, transform, min_pairs, normalize, missing)

# now, filter to the long metrics
grid %<>% 
  filter(metric %in% c("MI", "zi_kendall", "wccor", "GENIE3", "distance_cor"))
# we can also just apply min_pairs post-hoc
grid %<>%
  dplyr::select(-min_pairs) %>%
  distinct()

# write the raw array
grid_file = "sh/grids/write_matrices_raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "matrices")

# now, check for which parameters are already complete
grid0 = grid %>%
  mutate(output_dir = file.path(base_dir, "matrices", accession, 
                                experiment, search),
         output_filename = paste0(quant_mode, 
                                  '-metric=', metric,
                                  '-scale=', scale,
                                  '-transform=', transform,
                                  '-normalize=', normalize,
                                  '-missing=', missing,
                                  '.rds'),
         output_file = file.path(output_dir, output_filename),
         exists = file.exists(output_file),
         idx = row_number()) %>%
  filter(!exists) %>%
  dplyr::select(-output_filename, -exists, -idx) %>% 
  # set up input file
  mutate(input_file = file.path("~/git/CFdb-searches/data/chromatograms_gene", 
                                accession, experiment, search, 
                                paste0(quant_mode, ".rds"))) %>% 
  filter(file.exists(input_file)) %>% 
  # temporary: remove MS1
  filter(quant_mode != 'MS1_intensity')

# write
grid_file = "sh/grids/write_matrices.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = TRUE)
write.table(grid0, grid_file, quote = FALSE, row.names = FALSE, sep = "\t")

# write the sh file dynamically
sh_file = '~/git/CFdb-analysis/sh/write-matrices.sh'
write_sh(job_name = 'write-matrices',
         sh_file = sh_file,
         grid_file = grid_file,
         inner_file = 'R/analysis_grid/inner-write-matrices.R',
         system = system,
         time = 24,
         cpus = 1,
         mem = 64)

# finally, run the job on whatever system we're on
submit_job(grid0, sh_file, args$allocation, system)
