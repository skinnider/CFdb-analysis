setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
args = list(); source("R/functions/detect_system.R")

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
expts0 = filter(expts, Species %in% c('Homo sapiens', 'Mus musculus'))

# list feature table files
output_dir = file.path(base_dir, 'transfer_learning', 'feature_tables')
output_dirnames = with(expts0, file.path(Accession, Replicate, "default"))
output_files = file.path(output_dir, output_dirnames, 'features-top-4.rds')

# read them all
features = map(seq_along(output_files), ~ {
  row = expts0[.x, ] %>% 
    dplyr::select(Accession, Replicate)
  file = output_files[.x]  
  print(file)
  readRDS(file) %>% 
    cbind(row, .)
}) %>% 
  bind_rows()

# save
output_file = file.path(base_dir, "transfer_learning", "features-table.rds")
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
saveRDS(features, output_file)
