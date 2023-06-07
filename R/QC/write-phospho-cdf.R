# Plot a CDF of phosphosites vs. # of fractions in each experiment and 
# on average.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(seqinr)

# read experiments
expts = read.csv("~/git/CFdb-searches/data/experiments.csv")

# list all phosphosite files
files = list.files('~/git/CFdb-searches/data/phosphorylation', 
                   pattern = '*.rds', full.names = T, recursive = T) %>%
  # ignore the metadata files
  extract(!grepl("phosphosite", .))

# extract accession/experiment/search/quantitation
split = strsplit(gsub("^.*CFdb-searches/", "", files), '/')
accessions = map_chr(split, 3)
experiments = map_chr(split, 4)
quant_modes = gsub("\\..*$", "", basename(files))

# read all FASTA files
fasta_files = list.files("~/git/CFdb-searches/data/fasta/filtered", 
                         pattern = '*.gz', full.names = T)
patts = with(expts, paste0(substr(Species, 1, 1), '.',
                           gsub("^.* ", "", Species))) %>%
  unique() %>%
  fct_recode('Cyanothece' = 'C.51142',
             'B.oleracea' = 'B.italica',
             'Anabaena.PCC7120' = 'A.7120',
             'Synechocystis' = 'S.6803',
             'S.typhimurium' = 'S.SL1344') %>%
  paste0(collapse = '|')
fasta_files %<>% extract(grepl(patts, .))
fastas = map(fasta_files, seqinr::read.fasta) %>%
  setNames(gsub("^.*-|\\.fasta.*$", "", basename(fasta_files)))
proteome_sizes = map_int(fastas, length)

# construct CDF for each file
cdf = data.frame()
for (file_idx in seq_along(files)) {
  file = files[file_idx]
  accession = accessions[file_idx]
  experiment = experiments[file_idx]
  quant_mode = quant_modes[file_idx]
  message('[', file_idx, '/', length(files), '] working on file: ', file)
  
  # read phosphosite metadata
  site_file = file.path(dirname(file), 'phosphosites.rds')
  sites = readRDS(site_file)
  
  # extract CDF for this file
  mat = readRDS(file)
  cdf0 = data.frame()
  for (n_fractions in seq(0, ncol(mat))) {
    mat0 = mat %>%
      replace(is.na(.), 0) %>%
      extract(rowSums(. > 0) >= n_fractions, , drop = F)
    n_phosphosites = nrow(mat0)
    # filter metadata
    sites0 = filter(sites, as.character(id) %in% rownames(mat0))
    
    # also calculate coverage as a % of proteome size
    if (nrow(mat0) > 0) {
      proteins = sites0$Protein %>% unique()
      protein_groups = sites0$Proteins %>% 
        strsplit(';') %>%
        unlist() %>%
        unique()
      n_proteins = length(proteins)
      n_protein_groups = length(protein_groups)
      experiment_clean = experiment %>% 
        fct_recode('Ghosts_6061_DDM_BiobasicSEC-2' = 'Ghosts_DDM_BiobasicSEC-2',
                   'Hemolysate_9330_old_SEC' = 'Hemolysate_old_9330_SEC',
                   'whole_rbc_IP_lysis' = 'whole_rbc_IP_lysis_6840',
                   'BN3to12' = 'BN3to12_RAW',
                   'BN4to16' = 'BN4to16_RAW') %>% 
        as.character()
      species = expts %>%
        filter(Accession == accession, Replicate == experiment_clean) %>%
        pull(Species) 
      key = paste0(substr(species, 1, 1), '.', gsub("^.* ", "", species)) %>%
        fct_recode('Cyanothece' = 'C.51142',
                   'B.oleracea' = 'B.italica',
                   'Anabaena.PCC7120' = 'A.7120',
                   'Synechocystis' = 'S.6803',
                   'S.typhimurium' = 'S.SL1344') %>%
        as.character()
      proteome_size = proteome_sizes[key]
      if (is.na(proteome_size)) 
        stop("couldn't get proteome size for species: ", species)
      coverage = length(protein_groups) / proteome_size
    } else {
      coverage = 0
    }
    
    # append
    cdf0 %<>% bind_rows(data.frame(n_fractions = n_fractions,
                                   n_phosphosites = n_phosphosites,
                                   n_proteins = n_proteins,
                                   n_protein_groups = n_protein_groups,
                                   coverage = coverage))
  }
  
  # add the rest of the metadata
  cdf0 %<>% mutate(accession = accession, experiment = experiment, 
                   quant_mode = quant_mode)
  
  # append to master CDF
  cdf %<>% bind_rows(cdf0)
}

# rearrange columns
cdf %<>% 
  dplyr::select(accession, experiment, quant_mode,
                n_fractions, n_phosphosites, n_proteins, n_protein_groups, 
                coverage) %>%
  # clean up quantitation mode
  mutate(quant_mode = chartr('_', ' ', quant_mode))

# save 
cdf_file = 'data/QC/phosphosites-CDF.rds'
saveRDS(cdf, cdf_file)
