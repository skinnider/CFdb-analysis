setwd("~/git/CFdb-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(seqinr)

## bee
fa = seqinr::read.fasta("~/git/CFdb-searches/data/fasta/filtered/UP000005203-A.mellifera.fasta.gz",
                         seqtype = 'AA', as.string = TRUE)
accns = map_chr(fa, ~ attr(., 'name')) %>% strsplit('\\|') %>% map_chr(2)
genes = map_chr(fa, ~ attr(., 'Annot')) %>% 
  map_chr(~ gsub("^.*GN=", "", .) %>% gsub(" .*$", "", .))
gn_map = setNames(genes, accns)

# read first set of replicates
dat1 = read_tsv("data/bee/raw/rep1to10_pg_matrix.tsv.gz")
attr(dat1, 'spec') = NULL
# extract intensities
mat = dat1 %>% 
  extract(, grepl('_', colnames(.))) %>% 
  as.matrix()
# handle infinite or NaN values
mat[!is.finite(mat)] = NA
mat[is.nan(mat)] = NA
# map to gene names
proteins = strsplit(dat1$Protein.Group, ';')
matches = map(proteins, ~ gn_map[.] %>% extract(!is.na(.)))
n_uniq = map_int(matches, n_distinct)
genes = map_chr(matches, ~ unique(.) %>% c(., NA) %>% head(1), .default = NA)
rownames(mat) = genes
mat %<>% extract(!is.na(genes), )

# keep unique rows only
genes0 = unique(na.omit(genes))
gene_mat = matrix(NA, nrow = length(genes0), ncol = ncol(mat),
                  dimnames = list(genes0, colnames(mat)))
n_fractions = rowSums(!is.na(mat) & is.finite(mat) & mat != 0)
out_map = data.frame()
for (gene in genes0) {
  idxs = which(rownames(mat) == gene)
  # pick the best protein for this replicate
  n_fractions0 = n_fractions[idxs]
  best = names(which(n_fractions0 == max(n_fractions0))) %>%
    dplyr::first()
  gene_mat[gene, ] = mat[best, ]
  out_map %<>% bind_rows(data.frame(gene = gene, protein_group = best))
}

# split matrices
reps = gsub("_.*$", "", colnames(gene_mat)) %>% paste0('R', .)
fractions = gsub("^.*_", "", colnames(gene_mat)) %>% 
  str_pad(width = 2, side = 'left', pad = '0') %>% 
  paste0("F", .)
colnames(gene_mat) = fractions
mats = map(unique(reps), ~ gene_mat[, reps == .x])
# filter proteins that are absent from each matrix
mats %<>% map(~ extract(.x, rowSums(is.finite(.x)) > 0, ))

# read additional replicate
dat2 = read_tsv("data/bee/raw/report.pg_matrix_rep1_DIA.tsv.gz")
attr(dat2, 'spec') <- NULL
dat2_mat = dat2 %>% 
  extract(, grepl("\\.d", colnames(.))) %>% 
  as.matrix()
# handle infinite or NaN values
dat2_mat[!is.finite(dat2_mat)] = NA
dat2_mat[is.nan(dat2_mat)] = NA
# map to gene names
proteins = strsplit(dat2$Protein.Group, ';')
matches = map(proteins, ~ gn_map[.] %>% extract(!is.na(.)))
n_uniq = map_int(matches, n_distinct)
genes = map_chr(matches, ~ unique(.) %>% c(., NA) %>% head(1), .default = NA)
rownames(dat2_mat) = genes
dat2_mat %<>% extract(!is.na(genes), )

# clean up column names
colnames(dat2_mat) %<>% gsub("^.*Slot1-", "", .) %>% 
  gsub("_.*$", "", .) %>% 
  paste0('F', .)

# keep unique rows only
genes0 = unique(na.omit(genes))
gene_mat = matrix(NA, nrow = length(genes0), ncol = ncol(dat2_mat),
                  dimnames = list(genes0, colnames(dat2_mat)))
n_fractions = rowSums(!is.na(dat2_mat) & is.finite(dat2_mat) & dat2_mat != 0)
out_map = data.frame()
for (gene in genes0) {
  idxs = which(rownames(dat2_mat) == gene)
  # pick the best protein for this replicate
  n_fractions0 = n_fractions[idxs]
  best = names(which(n_fractions0 == max(n_fractions0))) %>%
    dplyr::first()
  # matrixStats::rowSds(mat[proteins, ], na.rm = T)
  gene_mat[gene, ] = dat2_mat[best, ]
  # save the mapping
  out_map %<>% bind_rows(data.frame(gene = gene, protein_group = best))
}

# save all matrices
out = mats
names(out) = unique(reps)
out %<>% c(list(R11 = gene_mat))
saveRDS(out, "data/bee/chromatograms.rds")
