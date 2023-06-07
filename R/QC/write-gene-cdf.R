# For each CF-MS experiment, write the total number of fractions 
# that each gene was detected in. 
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# read experiments
expts = read.csv("~/git/CFdb-searches/data/experiments.csv")

# iterate through
gene_cdf = pmap_dfr(expts, function(...) {
  current = tibble(...) 
  print(current)
  
  # set up chromatogram file
  chrom_dir = "~/git/CFdb-searches/data/chromatograms_gene"
  chrom_file = current %>%
    mutate(Replicate = fct_recode(Replicate, 
                                  'BN3to12_RAW' = 'BN3to12',
                                  'BN4to16_RAW' = 'BN4to16',
                                  'Ghosts_DDM_BiobasicSEC-2' = 'Ghosts_6061_DDM_BiobasicSEC-2',
                                  'Hemolysate_old_9330_SEC' = 'Hemolysate_9330_old_SEC',
                                  'whole_rbc_IP_lysis_6840' = 'whole_rbc_IP_lysis',
                                  'Hemolysate_8994_IEX-3' = 'Hemolysate_8994_IEX') %>% 
             as.character()) %$%
    file.path(chrom_dir, Accession, Replicate, 'default', 
              ifelse(Quantitation == 'TMT', 
                     'Reporter_intensity.rds', 
                     'MS1_intensity.rds'))
  
  # read matrix
  mat = readRDS(chrom_file)
  # count number of fractions for each gene
  fractions = rowSums(is.finite(mat) & mat > 0)
  # create data frame
  data.frame(gene = names(fractions), n_fractions = fractions) %>% 
    cbind(current, .) %>% 
    mutate(quant_mode = gsub("\\.rds$", "", basename(chrom_file)),
           search = 'default')
})

# rearrange columns
gene_cdf %<>% 
  dplyr::select(Accession, Replicate, Version, quant_mode, search,
                Species, gene, n_fractions) %>% 
  dplyr::rename(accession = Accession, 
                experiment = Replicate,
                version = Version, 
                species = Species) %>% 
  # clean up quantitation mode
  mutate(quant_mode = chartr('_', ' ', quant_mode))

# save 
cdf_file = 'data/QC/genes-CDF.rds'
saveRDS(gene_cdf, cdf_file)
