setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-run-EGAD.R')
grid = read.delim("sh/grids/run-EGAD.txt")
for (param_name in colnames(grid))
  parser$add_argument(paste0('--', param_name),
                      type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(PrInCE)
library(flavin)
library(igraph)
library(ontologyIndex)
library(EGAD)

# un-escape species from argparse
args$species %<>% chartr('_', ' ', .)
# retrieve from species map
species_map = read.csv("data/GO/species-map.csv") %>% 
  filter(species == args$species)

# read the network
if (grepl("\\.rds$", args$input_file)) {
  network = readRDS(args$input_file)
} else if (grepl("\\.txt", args$input_file)) {
  network = read.delim(args$input_file)
} else if (grepl("\\.csv", args$input_file)) {
  network = read.csv(args$input_file)
}

# read gene sets from ...
### GO
ontology = get_ontology("data/GO/go-basic.obo.gz")
goa_file = paste0("data/GO/GOA/", species_map$file, ".gz")
goa = read_gaf(goa_file,
               filter.NOT = T,
               filter.evidence = c("ND", "IPI", "IEA"), 
               ontology = ontology)
# remove roots 
rootNames = c(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")
goa %<>% dplyr::filter(!GO_ID %in% rootNames)

# create annotations
ann1 = as_annotation_list(goa, 'DB_Object_Symbol', 'GO_ID') %>%
  setNames(paste0('GO|', names(.)))

# for human, also read disease genes
if (args$species == 'Homo sapiens') {
  ### disease genes
  dis_files = list.files("data/resources/disease_genes",
                         pattern = '*\\.txt\\.gz', full.names = TRUE)
  dis_names = gsub("\\.txt\\.gz$", "", basename(dis_files))
  dis = map(dis_files, read.delim) %>%
    setNames(dis_names) %>%
    bind_rows(.id = 'database') %>%
    unite(gene_set, database, source, disease, sep = '|', remove = FALSE)
  
  # compile list of gene sets
  ann2 = as_annotation_list(dis, 'gene', 'gene_set') %>%
    setNames(paste0('disease_gene|', names(.)))
  gene_sets = c(ann1, ann2)
} else {
  gene_sets = ann1
}

# save the original size of every gene set
gene_set_sizes = data.frame(gene_set = names(gene_sets),
                            n_genes = lengths(gene_sets))
# standardize column names
colnames(network)[c(1, 2)] = c('gene_A', 'gene_B')

# alphabetize interactors
sorted = t(apply(network[, 1:2], 1, sort))
network$gene_A = sorted[, 1]
network$gene_B = sorted[, 2]
network %<>% distinct()
# tibbles break EGAD
network %<>% as.data.frame()

# extract relevant gene sets
g = graph_from_data_frame(network, directed = FALSE)
nodes = names(V(g))
ann0 = map(gene_sets, ~ intersect(., nodes)) %>%
  # doesn't make sense to look for anything less than 3 nodes
  extract(lengths(.) >= 3)

## EGAD AUC
message("calculating EGAD AUCs ...")
# make EGAD network
genelist = make_genelist(as.data.frame(network))
gene_network = make_gene_network(network, genelist)
# make annotations
ann_long = map(ann0, ~ data.frame(gene = .)) %>%
  bind_rows(.id = 'gene_set') %>%
  dplyr::select(gene, gene_set)
annotations = make_annotations(ann_long, genelist, names(ann0))
# run GBA with no threshold
gba = run_GBA(gene_network, annotations, min = 3, max = 1e5)
# set up output data frame
aurocs = gba[[1]][, "auc"]
# get number of proteins to which each term was annotated
all_terms = colSums(annotations)
n_terms = all_terms[names(aurocs)]
# get result
EGAD = data.frame(network = basename(args$input_file),
                  term = names(aurocs), auroc = aurocs, 
                  n_proteins = n_terms, 
                  pct_proteins = n_terms / length(genelist))

# save gene set sizes, too
attr(EGAD, 'gene_set_sizes') = gene_set_sizes

# save output
output_dir = dirname(args$output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
saveRDS(EGAD, args$output_file)
