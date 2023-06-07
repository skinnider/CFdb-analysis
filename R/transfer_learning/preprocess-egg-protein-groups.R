# Preprocess high-coverage bee proteomic data from McAfee et al., PLoS ONE
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(data.table)
args = list(); source("R/functions/detect_system.R")

# read protein groups
dat = read_delim(file.path(base_dir, "bee", "MSV000086862", "txt_eggs",
                           "proteinGroups_eggs.txt"))

# remove reverse, contaminant, etc.
dat0 = dat %>% 
  filter(is.na(Reverse), 
         is.na(`Only identified by site`),
         is.na(`Potential contaminant`))

# read idmapping
idmap = fread(file.path(base_dir, 'bee/idmapping', 'refseq.dat'),
              col.names = c('uniprot', 'db', 'id'), header = FALSE)
refseq = filter(idmap, db == 'RefSeq')

# filter to bee proteins
fa = seqinr::read.fasta("~/git/CFdb-searches/data/fasta/filtered/UP000005203-A.mellifera.fasta.gz",
                         seqtype = 'AA', as.string = TRUE)
accns = map_chr(fa, ~ attr(., 'name')) %>% strsplit('\\|') %>% map_chr(2)
refseq0 = filter(refseq, uniprot %in% accns)

# extract gene intensity matrix
mat = dat0 %>% 
  dplyr::select(starts_with('LFQ')) %>% 
  as.matrix()
genes = dat0$`Majority protein IDs` %>% 
  gsub(";.*$", "", .)
rownames(mat) = genes

# map refseq->uniprot->gene
map_df = data.frame(refseq = genes)
genes = map_chr(fa, ~ attr(., 'Annot')) %>% 
  map_chr(~ gsub("^.*GN=", "", .) %>% gsub(" .*$", "", .))
gn_map = data.frame(gene_name = genes, uniprot = accns)
refseq_map = refseq0 %>% 
  dplyr::select(uniprot, id) %>% 
  dplyr::rename(refseq = id)
map_df %<>%
  left_join(refseq_map, by = 'refseq') %>% 
  left_join(gn_map, by = 'uniprot') %>% 
  distinct(refseq, gene_name) %>% 
  drop_na()
n_distinct(map_df$refseq) == nrow(map_df) ## TRUE

# map the matrix to gene names
mat = mat[map_df$refseq, ] %>% 
  set_rownames(map_df$gene)

# log-transform
mat[mat == 0] = NA
mat %<>% log()
# keep proteins in 10+ samples
mat %<>% extract(rowSums(is.finite(.)) >= 10, )

# save
saveRDS(mat, "data/bee/protein-groups-egg.rds")
