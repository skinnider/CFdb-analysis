# Summarize key statistics about CFdb.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(data.table)
source("R/theme.R")

# list all summary files
files = list.files('data/maxquant', pattern = 'summary.txt', 
                   full.names = T, recursive = T) %>%
  # ignore phosphorylation
  extract(!grepl("phosphorylation", .))

# loop over files
summary = map_dfr(files, ~ {
  tab = fread(.) %>% 
    filter(Experiment != '')
  
  ## number of input raw files
  n_files = n_distinct(tab$`Raw file`)
  
  ## number of fractions
  n_fractions = n_distinct(tab$Experiment)
  
  ## number of MS/MS spectra
  msms = sum(tab$`MS/MS`)
  
  ## number of sequenced peptides
  peptides = sum(tab$`MS/MS Identified`)
  
  data.frame(file = ., raw_files = n_files, fractions = n_fractions,
             msms = msms, peptides = peptides)
})

# flag accession, replicate, and search
summary %<>%
  mutate(file = gsub("^.*maxquant/", "", file)) %>%
  separate(file, c("accession", "replicate", "search", "x"), '/') %>%
  dplyr::select(-x)

# merge in versions 1 vs 2
expts = read.csv("~/git/CFdb-searches/data/experiments.csv") %>% 
  # drop channels
  mutate(Replicate = gsub("-heavy|-medium|-[1-6]$", "", Replicate) %>% 
           fct_recode('BN3to12_RAW' = 'BN3to12',
                      'BN4to16_RAW' = 'BN4to16',
                      'Ghosts_DDM_BiobasicSEC-2' = 'Ghosts_6061_DDM_BiobasicSEC',
                      'Hemolysate_old_9330_SEC' = 'Hemolysate_9330_old_SEC',
                      'whole_rbc_IP_lysis_6840' = 'whole_rbc_IP_lysis',
                      'Hemolysate_8994_IEX-3' = 'Hemolysate_8994_IEX') %>% 
           as.character()) %>% 
  # fix one species
  mutate(Species = fct_recode(Species, 'Salmonella enterica' = 
                                'Salmonella typhimurium SL1344') %>% 
           as.character())
colnames(expts) %<>% tolower()
expts %<>% left_join(summary, by = c('accession', 'replicate'))

# print some statistics
message("number of experiments searched: ", nrow(expts))
message("number of raw files searched: ",
        summary %>% filter(search == 'default') %>% pull(raw_files) %>% sum())
message("number of fractions searched: ",
        summary %>% filter(search == 'default') %>% pull(fractions) %>% sum())
message("number of MS/MS spectra analyzed: ",
        summary %>% filter(search == 'default') %>% pull(msms) %>% sum())
message("number of peptides sequenced: ",
        summary %>% filter(search == 'default') %>% pull(peptides) %>% sum())

# compare with 2021 meta-analysis
expts1 = filter(expts, version == 'V1')
summary1 = inner_join(summary, distinct(expts1, accession, replicate))
message("number of experiments searched: ", nrow(expts1))
message("number of raw files searched: ",
        summary1 %>% filter(search == 'default') %>% pull(raw_files) %>% sum())
message("number of fractions searched: ",
        summary1 %>% filter(search == 'default') %>% pull(fractions) %>% sum())
message("number of MS/MS spectra analyzed: ",
        summary1 %>% filter(search == 'default') %>% pull(msms) %>% sum())
message("number of peptides sequenced: ",
        summary1 %>% filter(search == 'default') %>% pull(peptides) %>% sum())
