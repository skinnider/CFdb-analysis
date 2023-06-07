# Extract one-to-one bee/fly orthologs from InParanoid, and map them to their
# gene names.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(XML)
library(seqinr)

# read InParanoid
data = xmlParse("data/bee/ortholog/orthologs_7460_7227.xml.gz")
list = xmlToList(data)

# get the two species
species1 = list[[1]]
species1Name = species1$.attrs['name']
species2 = list[[2]]
species2Name = species2$.attrs['name']

# extract ortholog group IDs
ortholog_groups = list$groups
ortholog_scores = map_chr(ortholog_groups, c("score", "value")) %>% as.numeric()
ortholog_refs = purrr::map(ortholog_groups, ~ .[names(.) == "geneRef"])
ortholog_ids = lapply(ortholog_refs, function(x) map(x, c(".attrs", "id")))

# extract ID to protein maps 
db1 = species1$database
db2 = species2$database
genes1 = db1$genes
genes2 = db2$genes
ids1 = unname(purrr::map_chr(genes1, "id"))
ids2 = unname(purrr::map_chr(genes2, "id"))
prots1 = unname(purrr::map_chr(genes1, "protId"))
prots2 = unname(purrr::map_chr(genes2, "protId"))
prots1 = setNames(prots1, ids1)
prots2 = setNames(prots2, ids2)

# map ortholog IDs to proteins
ortholog_prots1 = purrr::map(ortholog_ids, ~ na.omit(prots1[unlist(.)]))
ortholog_prots2 = purrr::map(ortholog_ids, ~ na.omit(prots2[unlist(.)]))

# map to data frame
orthologs = purrr::map_dfr(seq_len(length(ortholog_ids)), ~ 
                             tidyr::crossing(source = ortholog_prots1[[.]],
                                             target = ortholog_prots2[[.]]))

# filter to one-to-one orthologs
single_source = orthologs %>%
  group_by(source) %>%
  dplyr::summarise(n_source = n_distinct(target)) %>%
  ungroup() %>%
  dplyr::filter(n_source == 1) %>%
  pull(source)
single_target = orthologs %>%
  group_by(target) %>%
  dplyr::summarise(n_target = n_distinct(source)) %>%
  ungroup() %>%
  dplyr::filter(n_target == 1) %>%
  pull(target)
orthologs %<>% dplyr::filter(source %in% single_source & 
                               target %in% single_target) %>%
  dplyr::select(source, target) %>% 
  map_dfc(unname)

# map UniProt IDs to gene names
## human
fa1 = seqinr::read.fasta("~/git/CFdb-searches/data/fasta/filtered/UP000000803-D.melanogaster.fasta.gz",
                         seqtype = 'AA', as.string = TRUE)
accns1 = map_chr(fa1, ~ attr(., 'name')) %>% strsplit('\\|') %>% map_chr(2)
genes1 = map_chr(fa1, ~ attr(., 'Annot')) %>% 
  map_chr(~ gsub("^.*GN=", "", .) %>% gsub(" .*$", "", .))
map1 = setNames(genes1, accns1)
## bee
fa2 = seqinr::read.fasta("~/git/CFdb-searches/data/fasta/filtered/UP000005203-A.mellifera.fasta.gz",
                         seqtype = 'AA', as.string = TRUE)
accns2 = map_chr(fa2, ~ attr(., 'name')) %>% strsplit('\\|') %>% map_chr(2)
genes2 = map_chr(fa2, ~ attr(., 'Annot')) %>% 
  map_chr(~ gsub("^.*GN=", "", .) %>% gsub(" .*$", "", .))
map2 = setNames(genes2, accns2)
orthologs %<>%
  mutate(source_gene = map2[source],
         target_gene = map1[target]) %>% 
  drop_na(source_gene, target_gene)
  
# do one more round of one-to-one filtering 
# (few dozen UniProt IDs map to >1 gene)
orthologs %<>%
  group_by(source, target) %>%
  filter(n_distinct(source_gene) == 1,
         n_distinct(target_gene) == 1) %>%
  ungroup()

# write
write.table(orthologs, "data/bee/ortholog/InParanoid-bee-fly.txt",
            row.names = F, sep = "\t", quote = F)
