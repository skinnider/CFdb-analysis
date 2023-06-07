setwd("~/git/CFdb-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)

# list eggNOG bact files
bact_files = list.files("data/resources/eggNOG/bact", full.names = TRUE, 
                        pattern = "gz")
proteomes = basename(bact_files) %>% gsub("-.*$", "", .)
ortho = map(bact_files, ~ 
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

# read species to taxid map
species_map = read.csv("data/GO/species-map.csv") %>% 
  # filter to bacterial
  filter(proteome %in% proteomes)

# read all experiments
expts = read.csv("~/git/CFdb-searches/data/experiments.csv") %>% 
  # filter to bacterial
  filter(Species %in% species_map$species)

# read proteomes
proteome_files = list.files("~/git/CFdb-searches/data/fasta/filtered", pattern = "gz", 
                            full.names = TRUE) %>% 
  extract(grepl(paste0(proteomes, collapse = '|'), .))
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

# read EcoCyc
complexes = read.delim("data/resources/EcoCyc/protcplxs.col.gz",
                       comment.char = '#')
complex_list = pmap(complexes, function(...) 
  dplyr::select(tibble(...), starts_with('GENE.NAME')) %>% 
    unlist() %>% 
    setdiff(''))
complex_names = complexes$NAME
complexes = setNames(complex_list, complex_names) %>% 
  extract(lengths(.) > 1)
complex_df = data.frame(complex = rep(names(complexes), lengths(complexes)),
                        gene_name = unlist(complexes))
# save E. coli
output_file = "data/prok_interactomes/EcoCyc/complexes_UP000000625.csv"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
write.csv(complex_df, output_file, row.names = FALSE)

# now, create eggNOG maps
for (proteome in setdiff(species_map$proteome, "UP000000625")) {
  message("working on proteome: ", proteome, " ...")
  bact1 = filter(ortho, proteome == !!proteome) 
  bact2 = filter(ortho, proteome == "UP000000625")
  
  # map to genes
  map1 = gn_maps[[proteome]]
  map2 = gn_maps[["UP000000625"]]
  bact1 %<>% left_join(map1, by = 'uniprot') %>% 
    drop_na(gene)
  bact2 %<>% left_join(map2, by = 'uniprot') %>% 
    drop_na(gene)
  
  # create one-to-one map
  bact1 %<>% 
    dplyr::select(-proteome) %>% 
    dplyr::rename(target = gene) %>% 
    # filter to single-mapping
    group_by(eggNOG) %>% 
    filter(n_distinct(target) == 1) %>% 
    ungroup()
  bact2 %<>% 
    dplyr::select(-proteome, -uniprot) %>% 
    dplyr::rename(ecoli = gene) %>% 
    # filter to single-mapping
    group_by(eggNOG) %>% 
    filter(n_distinct(ecoli) == 1) %>% 
    ungroup()
  map = left_join(bact1, bact2, by = 'eggNOG') %>% 
    dplyr::select(eggNOG, ecoli, target, uniprot) %>% 
    drop_na(ecoli, target) %>% 
    distinct(ecoli, target, uniprot) 
  
  # map CORUM
  mapped = complex_df %>% 
    dplyr::rename(ecoli = gene_name) %>% 
    left_join(map, by = 'ecoli') %>% 
    distinct(complex, target) %>% 
    drop_na()
  message("  mapped ", nrow(mapped), " of ", nrow(complex_df),
          " complex proteins")
  
  # save
  output_file = paste0("data/prok_interactomes/EcoCyc/complexes_", 
                       proteome, ".csv")
  write.csv(mapped, output_file, row.names = FALSE)
}
