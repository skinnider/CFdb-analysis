setwd("~/git/CFdb-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)

# list eggNOG euk files
euk_files = list.files("data/resources/eggNOG/euk", full.names = TRUE, 
                       pattern = "gz")
proteomes = basename(euk_files) %>% gsub("-.*$", "", .)
ortho = map(euk_files, ~ 
              read.delim(.,
                         comment.char = '#', header = FALSE,
                         col.names = c('query_name', 'seed_eggNOG_ortholog',
                                       'seed_ortholog_evalue', 'seed_ortholog_score',
                                       'predicted_gene_name', 'GO_terms', 
                                       'KEGG_KOs', 'BiGG_reactions',
                                       'Annotation_tax_scope', 'OGs', 
                                       'bestOG|evalue|score', 'COG cat', 
                                       'eggNOG annot')) %>%
              dplyr::select(query_name, OGs) %>%
              mutate(uniprot = map_chr(strsplit(query_name, '\\|'), 2),
                     eggNOG = strsplit(OGs, ',')) %>%
              unnest(eggNOG) %>% 
              # remove isoforms
              mutate(uniprot = gsub("-.*$", "", uniprot)) %>% 
              distinct(uniprot, eggNOG)) %>% 
  setNames(proteomes) %>% 
  bind_rows(.id = 'proteome')

# read all experiments
expts = read.csv("~/git/CFdb-searches/data/experiments.csv")

# read species to taxid map
species_map = read.csv("data/GO/species-map.csv")

# read proteomes
proteome_files = list.files("~/git/CFdb-searches/data/fasta/filtered",
                            pattern = "gz", full.names = TRUE)
proteomes = gsub("-.*$", "", basename(proteome_files))
fastas = map(proteome_files, ~ seqinr::read.fasta(., seqtype = 'AA', 
                                                  as.string = TRUE)) %>% 
  setNames(proteomes)
# create gene name maps
gn_maps = map(fastas, ~ {
  fa = .x
  accns = map_chr(fa, ~ attr(., 'name')) %>% 
    strsplit('\\|') %>% 
    map_chr(2)
  has_gn = map_chr(fa, ~ attr(., 'Annot')) %>% grepl('GN=', .)
  genes = map_chr(fa, ~ attr(., 'Annot')) %>% 
    map_chr(~ gsub("^.*GN=", "", .) %>% gsub(" .*$", "", .))
  genes[!has_gn] = NA
  gn_map = data.frame(uniprot = accns, gene = genes)
  
  # deal with four species without gene names
  if (mean(!is.na(gn_map$gene)) < 0.5)
    gn_map %<>% mutate(gene = uniprot)
  
  drop_na(gn_map)
})

# read CORUM
complexes = read.delim("data/complex/CORUM/complexes_human.txt")

# now, create eggNOG maps
for (proteome in setdiff(species_map$proteome, "UP000005640")) {
  euk1 = filter(ortho, proteome == !!proteome) 
  euk2 = filter(ortho, proteome == "UP000005640")
  
  # map to genes
  map1 = gn_maps[[proteome]]
  map2 = gn_maps[["UP000005640"]]
  euk1 %<>% left_join(map1, by = 'uniprot') %>% 
    drop_na(gene)
  euk2 %<>% left_join(map2, by = 'uniprot') %>% 
    drop_na(gene)
  
  # create one-to-one map
  euk1 %<>% 
    dplyr::select(-proteome) %>% 
    dplyr::rename(target = gene) %>% 
    # filter to single-mapping
    group_by(eggNOG) %>% 
    filter(n_distinct(target) == 1) %>% 
    ungroup()
  euk2 %<>% 
    dplyr::select(-proteome, -uniprot) %>% 
    dplyr::rename(human = gene) %>% 
    # filter to single-mapping
    group_by(eggNOG) %>% 
    filter(n_distinct(human) == 1) %>% 
    ungroup()
  map = left_join(euk1, euk2, by = 'eggNOG') %>% 
    dplyr::select(eggNOG, human, target, uniprot) %>% 
    drop_na(human, target) %>% 
    distinct(human, target, uniprot) 
  
  # map CORUM
  mapped = complexes %>% 
    dplyr::rename(human = gene_name) %>% 
    left_join(map, by = 'human') %>% 
    distinct(complex, target) %>% 
    drop_na()
  # save
  output_file = paste0("data/euk_interactomes/CORUM/complexes_", 
                       proteome, ".txt")
  output_dir = dirname(output_file)
  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)
  write.csv(mapped, output_file, row.names = FALSE)
}
