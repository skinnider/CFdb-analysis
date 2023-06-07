setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-run-EGAD-bee.R')
grid = read.delim("sh/grids/run-EGAD-bee.txt")
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

# read the network
if (grepl("\\.rds$", args$input_file)) {
  network = readRDS(args$input_file)
} else if (grepl("\\.txt", args$input_file)) {
  network = read.delim(args$input_file)
} else if (grepl("\\.csv", args$input_file)) {
  network = read.csv(args$input_file)
}

# read gene sets from GO
ontology = get_ontology("data/GO/go-basic.obo.gz")
gaf.colnames = c("DB", "DB_Object_ID", "DB_Object_Symbol", 
                 "Qualifier", "GO_ID", "DB_Reference", "Evidence_Code", 
                 "With_From", "Aspect", "DB_Object_Name", "DB_Object_Synonym", 
                 "DB_Object_Type", "Taxon", "Date", "Assigned_By", 
                 "Annotation_Extension", "Gene_Product_Form_ID")
goa = readr::read_tsv("data/GO/GOA/4859278.A_mellifera.goa.gz", comment = "!", 
                      col_names = gaf.colnames)
## filter.NOT
goa = goa[!grepl("NOT", goa$Qualifier), ]
## skip filter.evidence for bee
## propagate annotations
goa$ancestors = ontology$ancestors[goa$GO_ID]
goa = goa[lengths(goa$ancestors) > 0, ]
goa = tidyr::unnest(goa, ancestors)
goa[["GO_ID"]] = goa[["ancestors"]]
goa = goa[, -ncol(goa)]
# remove roots 
rootNames = c(BP = "GO:0008150", CC = "GO:0005575", MF = "GO:0003674")
goa %<>% dplyr::filter(!GO_ID %in% rootNames)

# map UniProt IDs to genes
fa = seqinr::read.fasta("~/git/CFdb-searches/data/fasta/filtered/UP000005203-A.mellifera.fasta.gz",
                        seqtype = 'AA', as.string = TRUE)
accns = map_chr(fa, ~ attr(., 'name')) %>% strsplit('\\|') %>% map_chr(2)
genes = map_chr(fa, ~ attr(., 'Annot')) %>% 
  map_chr(~ gsub("^.*GN=", "", .) %>% gsub(" .*$", "", .))
gn_map = setNames(genes, accns)
goa %<>% mutate(gene = gn_map[DB_Object_ID]) %>% 
  drop_na(gene)

# create annotations
ann1 = as_annotation_list(goa, 'gene', 'GO_ID') %>%
  setNames(paste0('GO|', names(.)))
gene_sets = ann1

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
