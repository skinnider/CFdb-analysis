# Create a chromatogram matrix containing the top-100 most frequently quantified
# human phosphosites across all experiments.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# get most frequently quantified sites
top_sites = readRDS("data/QC/phosphosites-localization-scores.rds") %>% 
  filter(n_fractions > 0) %>% 
  filter(species == 'Homo sapiens') %>% 
  group_by(gene) %>% 
  summarise(fractions = sum(n_fractions)) %>% 
  ungroup() %>% 
  arrange(desc(fractions)) %>% 
  head(100)

# extract tidy chromatogram df for these sites
expts = read.csv("~/git/CFdb-searches/data/experiments.csv") %>% 
  filter(Species == 'Homo sapiens') %>% 
  # drop duplicates
  filter(!(grepl("-[2-6]$", Replicate) & Accession %in% c('PXD027704',
                                                          'PXD023211',
                                                          'MSV000082468')))
chrom_df = pmap_dfr(expts, function(...) {
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
  
  # subset matrix
  keep = which(sites %in% top_sites$gene)
  mat0 = mat[keep, , drop = FALSE]
  sites0 = with(meta[keep, ], paste0(Protein, '-', `Amino acid`, '-', Position))
  if (nrow(mat0) == 0) return(data.frame())
  rownames(mat0) = sites0
  
  # convert to long
  df = reshape2::melt(mat0, varnames = c('gene', 'fraction'), 
                      value.name = 'intensity', as.is = TRUE) %>% 
    cbind(current, .) %>% 
    mutate(quant_mode = gsub("\\.rds$", "", basename(chrom_file)))
})

# map to gene names
fa = seqinr::read.fasta("~/git/CFdb-searches/data/fasta/filtered/UP000005640-H.sapiens.fasta.gz",
                        seqtype = 'AA', as.string = TRUE)
accns = map_chr(fa, ~ attr(., 'name')) %>% 
  strsplit('\\|') %>% 
  map_chr(2)
genes = map_chr(fa, ~ attr(., 'Annot')) %>% 
  map_chr(~ gsub("^.*GN=", "", .) %>% gsub(" .*$", "", .))
# merge
chrom_df %<>% 
  mutate(uniprot = strsplit(gene, '-') %>% 
           map_chr(~ head(., -2) %>% paste0(collapse = '-')),
         position = strsplit(gene, '-') %>% 
           map_chr(~ tail(., 1))) %>% 
  left_join(data.frame(uniprot = accns, gene_name = genes), by = 'uniprot')

# write
output_file = "data/phosphorylation/frequently-quantified-sites.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
saveRDS(chrom_df, output_file)
