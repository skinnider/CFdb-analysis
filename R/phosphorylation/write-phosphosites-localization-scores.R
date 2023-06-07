setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# read experiments
expts = read.csv("~/git/CFdb-searches/data/experiments.csv") %>% 
  # drop duplicates
  filter(!(grepl("-[2-6]$", Replicate) & Accession %in% c('PXD027704',
                                                          'PXD023211',
                                                          'MSV000082468')))

# iterate through
gene_cdf = pmap_dfr(expts, function(...) {
  current = tibble(...) 
  print(current)
  
  # set up chromatogram file
  chrom_dir = "~/git/CFdb-searches/data/phosphorylation"
  chrom_file = current %>%
    mutate(Replicate = fct_recode(Replicate, 
                                  'BN3to12_RAW' = 'BN3to12',
                                  'BN4to16_RAW' = 'BN4to16',
                                  'Ghosts_DDM_BiobasicSEC-2' = 'Ghosts_6061_DDM_BiobasicSEC-2',
                                  'Hemolysate_old_9330_SEC' = 'Hemolysate_9330_old_SEC',
                                  'whole_rbc_IP_lysis_6840' = 'whole_rbc_IP_lysis',
                                  'Hemolysate_8994_IEX-3' = 'Hemolysate_8994_IEX') %>% 
             as.character()) %$%
    file.path(chrom_dir, Accession, Replicate, 'intensity.rds')
  
  # read matrix
  mat = readRDS(chrom_file)
  
  # read metadata 
  meta_file = file.path(dirname(chrom_file), 'phosphosites.rds')
  meta = readRDS(meta_file)
  sites = with(meta, paste0(Protein, '-', `Amino acid`, '-', Position))
  
  # count number of fractions for each site
  fractions = rowSums(is.finite(mat) & mat > 0)
  # create data frame
  data.frame(gene = sites, n_fractions = fractions) %>% 
    mutate(localization_prob = meta$`Localization prob`,
           delta_score = meta$`Delta score`,
           intensity = meta$Intensity,
           ratio_mod_base = meta$`Ratio mod/base`) %>% 
    cbind(current, .) %>% 
    mutate(quant_mode = gsub("\\.rds$", "", basename(chrom_file)))
})

# rearrange columns
gene_cdf %<>% 
  dplyr::select(Accession, Replicate, Version, quant_mode, 
                Species, gene, n_fractions, everything()) %>% 
  dplyr::select(-Author, -Year, -Quantitation, -Fractionation) %>% 
  dplyr::rename(accession = Accession, 
                experiment = Replicate,
                version = Version, 
                species = Species) %>% 
  # clean up quantitation mode
  mutate(quant_mode = chartr('_', ' ', quant_mode))

# save 
cdf_file = 'data/QC/phosphosites-localization-scores.rds'
saveRDS(gene_cdf, cdf_file)
