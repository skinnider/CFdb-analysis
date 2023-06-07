# Merge features (gene-level distance correlation matrices) across all human
# experiments for input to UMAP.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
args = list(); source("R/functions/detect_system.R")

# read single best feature (distance correlation) parameters
## keep only the best feature for now
best_features = readRDS("data/analysis_grid/best_features.rds") %>%
  extract2('combined') %>%
  # keep only one metric
  group_by(metric) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  arrange(desc(median)) %>%
  # spread 'transformation' back out
  mutate(transform = ifelse(transformation == 'log-transform',
                            'log', 'none'),
         normalize = ifelse(transformation == 'quantile\nnormalization',
                            'quantile', 'none'))
params = head(best_features, 1)

# set up input files
chrom_dir = '~/git/CFdb-searches/data/chromatograms_gene'
chrom_files = list.files(chrom_dir, pattern = '*.rds', recursive = T) %>%
  # ignore the metadata files
  extract(!grepl("metadata", .))
split = strsplit(chrom_files, '/')
accessions = map_chr(split, 1)
experiments = map_chr(split, 2)
searches = map_chr(split, 3)
quant_modes = gsub("\\..*$", "", basename(chrom_files))
inputs = data.frame(file = file.path(chrom_dir, chrom_files),
                    accession = accessions,
                    experiment = experiments,
                    search = searches,
                    quant_mode = quant_modes) %>%
  # limit to default search only
  filter(search == 'default') %>%
  dplyr::select(-search)

# read experiments
expts = read.csv("~/git/CFdb-searches/data/experiments.csv") %>% 
  mutate(Replicate = fct_recode(Replicate,
                                'BN3to12_RAW' = 'BN3to12',
                                'BN4to16_RAW' = 'BN4to16',
                                'Ghosts_DDM_BiobasicSEC-2' = 'Ghosts_6061_DDM_BiobasicSEC-2',
                                'Hemolysate_old_9330_SEC' = 'Hemolysate_9330_old_SEC',
                                'whole_rbc_IP_lysis_6840' = 'whole_rbc_IP_lysis',
                                'Hemolysate_8994_IEX-3' = 'Hemolysate_8994_IEX')) %>% 
  dplyr::rename(experiment = Replicate, accession = Accession)
# filter to human
expts0 = filter(expts, Species == 'Homo sapiens')
inputs %<>% inner_join(expts0)

# set up correlation matrix filepaths
inputs %<>%
  # set up filepaths
  mutate(mat_dir = file.path(base_dir, "matrices", accession, 
                             experiment, "default"),
         mat_filename = paste0(quant_mode,
                               '-metric=', params$metric,
                               '-scale=', FALSE,
                               '-transform=', params$transform,
                               '-normalize=', params$normalize,
                               '-missing=', params$missing,
                               '.rds'),
         mat_file = file.path(mat_dir, mat_filename))

# read matrices
files = inputs$mat_file
mats = map(files, readRDS)

# get all genes
genes = map(mats, rownames) %>% Reduce(union, .)

# calculate feature matrices for 100 genes at a time
loop_size = 100
loop_len = ceiling(length(genes) / loop_size)
chunks = split(genes, ceiling(seq_along(genes) / loop_size))
feat_mats = map(seq_len(loop_len), ~ {
  message(".. working on chunk ", .x, " of ", loop_len, " ...")
  genes0 = chunks[[.x]]
  
  # match dimensions of each matrix
  new_mats = list()
  for (mat_idx in seq_along(mats)) {
    mat = mats[[mat_idx]]
    mat0 = mat[rownames(mat) %in% genes0, ]
    new_mat = matrix(NA, nrow = length(genes0), ncol = length(genes),
                     dimnames = list(genes0, genes))
    new_mat[rownames(mat0), colnames(mat0)] = mat0
    new_mats[[mat_idx]] = new_mat
  }
  
  # average correlations over all matrices
  array = simplify2array(new_mats)
  feat_mat = apply(array, c(1, 2), mean, na.rm = TRUE)
  return(feat_mat)
})

# rbind
feat_mat = do.call(rbind, feat_mats)

# remove any rows with all NAs
keep = rowSums(is.finite(feat_mat)) > 0
feat_mat %<>% extract(keep, keep)

# write matrix
output_file = file.path(base_dir, "UMAP/features/human.rds")
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
saveRDS(feat_mat, output_file)
