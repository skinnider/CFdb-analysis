# Calculate the number of CORUM complex proteins quantified in each 
# human and mouse dataset.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(seqinr)
source("R/theme.R")

# read CORUM 
corum_hs = read.delim("data/complex/CORUM/complexes_human.txt")
corum_mm = read.delim("data/complex/CORUM/complexes_mouse.txt")

# filter human complexes to bona fide human genes
fa = read.fasta("~/git/CFdb-searches/data/fasta/filtered/UP000005640-H.sapiens.fasta.gz",
  as.string = T, seqtype = 'AA')
human_genes = map_chr(fa, getAnnot) %>%
  gsub("^.*GN=", "", .) %>%
  gsub(" .*$", "", .) %>%
  unique() %>%
  # remove proteins without gene names
  extract(!startsWith(., ">"))
corum_hs %<>% filter(gene_name %in% human_genes)

# create list of complexes
complexes = list(
  `Homo sapiens` = list(CORUM = corum_hs),
  `Mus musculus` = list(CORUM = corum_mm)
)

# read gene CDF
cdf = readRDS("data/QC/genes-CDF.rds") %>% 
  filter(species != 'Arabidopsis thaliana')

# iterate through species
detection = tidyr::crossing(version = c('Version 1', 'Version 2')) %>% 
  pmap_dfr(function(...) {
    current = tibble(...)
    cdf0 = filter(cdf, species == 'Homo sapiens') %>% 
      filter(n_fractions > 0)
    complexes0 = complexes$`Homo sapiens`$CORUM
    if (current$version == 'Version 1') {
      cdf0 %<>% filter(version == 'V1')
    }
    n_fractions = cdf0 %>% 
      group_by(gene) %>% 
      summarise(n_fractions = sum(n_fractions)) %>% 
      ungroup()
    complexes0 %<>% 
      dplyr::rename(gene = gene_name) %>% 
      left_join(n_fractions, by = 'gene') %>% 
      replace_na(list(n_fractions = 0)) %>% 
      distinct(gene, n_fractions) %>% 
      mutate(detected = n_fractions > 0) %>% 
      cbind(current, .)
  }) %>% 
  group_by(version, gene) %>% 
  summarise(detected = any(detected),
            n_fractions = sum(n_fractions)) %>% 
  ungroup()
# print 
detection %>% group_by(version) %>% summarise(mean(detected))

# plot CDF
detection_v1 = detection %>% 
  filter(version == 'Version 1') %>% 
  filter(n_fractions > 0) %>% 
  arrange(desc(n_fractions)) %>% 
  mutate(gene = factor(gene, levels = unique(.$gene)),
         gene_idx = as.integer(gene))
detection_v2 = detection %>% 
  filter(version == 'Version 2') %>% 
  filter(n_fractions > 0) %>% 
  arrange(desc(n_fractions)) %>% 
  mutate(gene = factor(gene, levels = unique(.$gene)),
         gene_idx = as.integer(gene))
detection_version = bind_rows(
  mutate(detection_v1, version = 'Version 1'),
  mutate(detection_v2, version = 'Version 2')
)
version_pal = c('grey88', pals::stepped()[3])
p1 = detection_version %>% 
  ggplot(aes(x = n_fractions, y = gene_idx / n_distinct(detection$gene),
             color = version)) +
  geom_path(size = 0.3) + 
  scale_x_log10('# of fractions', labels = fancy_scientific) +
  scale_y_continuous('% of CORUM', labels = ~ . * 100) +
  annotation_logticks(sides = 'b', short = unit(0.06, 'lines'),
                      mid = unit(0.12, 'lines'), long = unit(0.18, 'lines'),
                      size = 0.2) +
  scale_color_manual('', values = version_pal,
                     breaks = c('Version 1', 'Version 2')) +
  boxed_theme() +
  theme(aspect.ratio = 1,
        legend.key.width = unit(0.4, 'lines'),
        legend.key.height = unit(0.2, 'lines'))
p1
# inset pie chart
n1 = nrow(detection_v1)
n2 = nrow(detection_v2) - n1
n0 = n_distinct(detection$gene) - n2 - n1
title = format(n2, big.mark = ',') %>% 
  paste('additional\nproteins detected')
p1_inset = data.frame(n = c(n1, n2, n0),
                      fill = c('Version 1', 'Version 2', 'Undetected')) %>% 
  mutate(title = title) %>% 
  ggplot(aes(x = '1', y = n, fill = fill)) +
  facet_grid(~ title) +
  geom_col(size = 0.15, color = 'grey20') + 
  coord_polar(theta = 'y') + 
  scale_fill_manual('', values = c('Version 1' = version_pal[1],
                                   'Version 2' = version_pal[2],
                                   'Undetected' = 'white')) +
  # guides(fill = guide_legend())
  boxed_theme() + 
  theme(plot.margin = margin(rep(0, 4)),
        plot.background = element_blank(),
        legend.position = 'none',
        legend.key.size = unit(0.35, 'lines'),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
p1_inset
# superimpose 
p1_comb = p1 + 
  inset_element(p1_inset, left = 0, bottom = 0, right = 0.5, top = 0.5,
                align_to = 'panel')
p1_comb
ggsave(paste0("fig/CORUM/cdf-version.pdf"), p1_comb, 
       width = 5.25, height = 5.4, units = "cm", useDingbats = FALSE)

# now, do GO enrichment on these proteins ####
map = read_tsv("data/identifier/HUMAN_9606_idmapping.dat.gz",
               col_names = c("uniprot", "db", "id"))
gn = filter(map, db == 'Gene_Name') %>% dplyr::select(-db)
id = filter(map, db == 'GeneID') %>% dplyr::select(-db)
gn2id = left_join(gn, id, by = 'uniprot') %>%
  dplyr::select(starts_with('id')) %>%
  set_colnames(c('gene', 'id'))
entrez = detection %>%
  left_join(gn2id, by = 'gene') %>%
  drop_na(id) %>%
  distinct(version, id, detected)

# GO analysis:
library(GOstats)
# - hits = genes detected in version 2 only
# - universe = genes detected in version 2
universe = filter(entrez, detected) %$% unique(id)
v1 = filter(entrez, detected, version == 'Version 1') %$% unique(id)
v2 = filter(entrez, detected, version == 'Version 2') %$% unique(id)
hits = setdiff(v2, v1)
# run GOstats enrichment analysis
enrs = list()
for (ontology in c("BP", "MF", "CC")) {
  message(".. analyzing enrichment of ", ontology, " terms ...")
  for (cond in c(T, F)) {
    message("... performing test with conditional parameter = ", cond)
    params = new("GOHyperGParams", 
                 geneIds = hits, universeGeneIds = universe,
                 annotation = "org.Hs.eg.db", ontology = ontology,
                 pvalueCutoff = 1, conditional = cond, 
                 testDirection = 'over')
    hgOver = hyperGTest(params)
    enr = summary(hgOver) %>%
      mutate(root = ontology, 
             conditional = cond)
    enrs[[length(enrs) + 1]] = enr
  }
}
# save
saveRDS(enrs, "data/complex/GO-never-detected.rds")

# manually filter redundant terms and plot final
enr = enrs %>%
  bind_rows() %>% 
  filter(!conditional) %>% 
  filter(Pvalue < 0.05)
enr %>% 
  arrange(Pvalue) %>% 
  pull(Term) %>% 
  head(100)
terms = c("integral component of plasma membrane",
          "ion channel activity",
          "transmembrane signaling receptor activity",
          "G protein-coupled receptor activity",
          "transporter complex",
          "synaptic membrane",
          "ligand-gated channel activity",
          "potassium ion transport",
          "GABA receptor complex",
          "blood vessel morphogenesis",
          "glycosyltransferase activity",
          "sphingosine-1-phosphate receptor activity",
          "bioactive lipid receptor activity")
pal = get_pal("early meadow")[c(1, 3)]
p2 = enr %>% 
  filter(!conditional) %>%
  filter(Term %in% terms) %>% 
  ggplot(aes(x = reorder(Term, -Pvalue), y = -log10(Pvalue))) + 
  # facet_grid(~ "Detected") +
  geom_col(width = 0.7, alpha = 0.7, size = 0.3, aes(color = '1', fill = '1')) +
  geom_hline(aes(yintercept = -log10(0.05)), color = 'red', size = 0.3) +
  geom_vline(aes(xintercept = -Inf), color = 'grey50', size = 0.5) +
  scale_fill_manual(values = pal[2], guide = F) +
  scale_color_manual(values = pal[2], guide = F) +
  scale_x_reordered() +
  scale_y_continuous(expression(-log[10](P)), limits = c(0, 20), 
                     expand = c(0, 0)) +
  coord_flip() +
  boxed_theme() +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 1.2)
p2
ggsave("fig/CORUM/GO-enrichment.pdf", p2, width = 10, height = 5, 
       units = "cm", useDingbats = FALSE)
