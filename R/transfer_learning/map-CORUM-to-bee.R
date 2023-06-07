setwd("~/git/CFdb-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(flavin)
library(PrInCE)

# read CORUM
corum = read.delim("data/complex/CORUM/complexes_human.txt") %>%
  as_annotation_list('gene_name', 'complex')

# read one-to-one orthologs
ortho = read.delim("data/bee/ortholog/InParanoid-bee-human.txt")
ortho_map = with(ortho, setNames(source_gene, target_gene))

# map CORUM
corum_bee = map(corum, ~ ortho_map[.] %>% extract(!is.na(.))) %>% 
  extract(lengths(.) >= 2)

## print number of mapped complexes
length(corum)
length(corum_bee)
## print number of mapped interactions
map_dfr(corum_bee, ~ tidyr::crossing(protein_A = ., protein_B = .) %>% 
          filter(protein_A < protein_B)) %>% 
  distinct() %>% 
  nrow()

# how many chromatogram genes are in bee CORUM?
chrom = readRDS("data/bee/chromatograms.rds")
chrom_genes = map(chrom, rownames) %>% Reduce(union, .)
corum_genes = unlist(corum_bee) %>% unique()
length(chrom_genes)
length(corum_genes)
length(intersect(chrom_genes, corum_genes)) ## 487!

# write 
bee_df = map2_dfr(names(corum_bee), corum_bee,
                  ~ data.frame(complex = .x, gene_name = .y))
write.csv(bee_df, "data/bee/complexes_bee.csv", row.names = FALSE)
