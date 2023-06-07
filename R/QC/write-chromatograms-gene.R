# Convert MaxQuant chromatograms at the protein group level to gene-level
# matrices, to enable comparisons across experiments.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# read experiments
expts = read.csv("~/git/CFdb-searches/data/experiments.csv")

# set up chromatogram files
chrom_dir = "~/git/CFdb-searches/data/chromatograms"
chrom_files = expts %>% 
  mutate(Replicate = Replicate %>% 
           fct_recode('BN3to12_RAW' = 'BN3to12',
                      'BN4to16_RAW' = 'BN4to16',
                      'Ghosts_DDM_BiobasicSEC-2' = 'Ghosts_6061_DDM_BiobasicSEC-2',
                      'Hemolysate_old_9330_SEC' = 'Hemolysate_9330_old_SEC',
                      'whole_rbc_IP_lysis_6840' = 'whole_rbc_IP_lysis',
                      'Hemolysate_8994_IEX-3' = 'Hemolysate_8994_IEX') %>% 
           as.character()) %$%
  file.path(chrom_dir, Accession, Replicate, 'default', 
            ifelse(Quantitation == 'TMT', 
                   'Reporter_intensity.rds', 
                   'iBAQ.rds'))
accessions = expts$Accession

# set up data frame to check gene name parsing
gn_parse = data.frame(file = chrom_files, genes = NA, protein_groups = NA,
                      example = NA)

# iterate through each in turn
for (file_idx in seq_along(chrom_files)) {
  chrom_file = chrom_files[file_idx]
  message("[", file_idx, "/", length(chrom_files), "] ", chrom_file, " ...")
  
  # read chromatogram matrix
  mat = readRDS(chrom_file)
  # handle infinite or NaN values
  mat[!is.finite(mat)] = NA
  mat[is.nan(mat)] = NA
  
  # read metadata
  meta_file = gsub("iBAQ|MS1_intensity|Reporter_intensity", "metadata", 
                   chrom_file)
  meta = readRDS(meta_file)
  
  # for Plasmodium datasets, filter out mouse proteins
  accession = accessions[file_idx]
  if (accession == 'PXD009039') {
    keep = !grepl('Mus musculus', meta$`Fasta headers`)
    mat %<>% extract(keep, )
    meta %<>% extract(keep, )
    meta$`Gene names` = NULL
  }
  
  # add gene names, if not already present
  if (!'Gene names' %in% colnames(meta)) {
    gene_names = ifelse(grepl('GN', meta$`Fasta headers`),
                        meta$`Fasta headers` %>% 
                          gsub("^.*GN=", "", .) %>% 
                          gsub(" .*$", "", .),
                        NA)
    
    # catch species without gene names (e.g. S. purpuratus)
    ratio = n_distinct(gene_names) / nrow(meta)
    if (ratio < 0.79) { ## Podospora anserina is 0.799 and is OK
      gene_names = meta$`Majority protein IDs` %>% 
        gsub(";.*$", "", .) %>% 
        gsub(" .*$", "", .) %>% 
        gsub("-.*$", "", .)
    }
    meta$`Gene names` = gene_names
    gn_parse$genes[file_idx] = n_distinct(gene_names)
    gn_parse$protein_groups[file_idx] = nrow(meta)
    gn_parse$example[file_idx] = head(gene_names, 10) %>% paste0(collapse = ' ')
    message(".. parsed ", n_distinct(gene_names), " gene names manually from ",
            nrow(meta), " protein groups")
  }
  
  # keep one row per gene
  gene_map = meta %>%
    dplyr::select(`Majority protein IDs`, `Gene names`) %>%
    set_colnames(c("protein_group", "gene")) %>%
    mutate(gene = strsplit(gene, ';')) %>%
    unnest(gene) %>%
    drop_na()
  genes = unique(gene_map$gene)
  gene_mat = matrix(NA, nrow = length(genes), ncol = ncol(mat),
                    dimnames = list(genes, colnames(mat)))
  n_fractions = rowSums(!is.na(mat) & is.finite(mat) & mat != 0)
  out_map = data.frame()
  for (gene in genes) {
    protein_groups = gene_map$protein_group[gene_map$gene == gene]
    # pick the best protein for this replicate
    n_fractions0 = n_fractions[protein_groups]
    best = names(which(n_fractions0 == max(n_fractions0))) %>%
      dplyr::first()
    # matrixStats::rowSds(mat[proteins, ], na.rm = T)
    gene_mat[gene, ] = mat[best, ]
    # save the mapping
    out_map %<>% bind_rows(data.frame(gene = gene, protein_group = best))
  }
  
  # set attribute on the gene matrix
  attr(gene_mat, 'gene_map') = out_map
  
  # save
  output_file = gsub("chromatograms", "chromatograms_gene", chrom_file)
  output_dir = dirname(output_file)
  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)
  saveRDS(gene_mat, output_file)
}
