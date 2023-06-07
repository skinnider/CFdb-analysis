# Plot the number of understudied proteins with phosphosites in CF-MS data.
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(igraph)
library(readxl)
source("R/theme.R")

# read gene2pubmed
g2p = read.delim("data/resources/uncharacterized/gene2pubmed.gz") %>% 
  filter(X.tax_id == 9606) %>% 
  dplyr::count(GeneID)
# map to gene names
map = read_tsv("data/identifier/HUMAN_9606_idmapping.dat.gz",
               col_names = c("uniprot", "db", "id"))
gn = filter(map, db == 'Gene_Name') %>% dplyr::select(-db)
id = filter(map, db == 'GeneID') %>% dplyr::select(-db)
gn2id = left_join(gn, id, by = 'uniprot') %>%
  dplyr::select(starts_with('id')) %>%
  set_colnames(c('gene', 'id'))
# execute mapping
g2p %<>%
  dplyr::rename(id = GeneID, n_pubs = n) %>% 
  mutate(id = as.character(id)) %>% 
  left_join(gn2id, by = 'id') %>% 
  drop_na(gene) %>% 
  distinct(gene, n_pubs)

# read phosphosites
phos = readRDS("data/QC/phosphosites-localization-scores.rds") %>% 
  filter(species == 'Homo sapiens')
# count the number of times each phosphosite was quantified
n_frac = phos %>% 
  group_by(gene) %>% 
  summarise(n = sum(n_fractions)) %>% 
  ungroup()
n_expt = phos %>% 
  filter(n_fractions > 0) %>% 
  dplyr::count(gene)
# extract 'frequently quantified' sites
freq = n_frac %>% 
  filter(n >= 5) %>% 
  pull(gene)
# what phosphoproteins are those p-sites on?
freq_prots = unique(gsub("-[STY]-.*$", "", freq)) %>% 
  unique()

# map to gene names
fa = seqinr::read.fasta("~/git/CFdb-searches/data/fasta/filtered/UP000005640-H.sapiens.fasta.gz",
                        seqtype = 'AA', as.string = TRUE)
accns = map_chr(fa, ~ attr(., 'name')) %>% 
  strsplit('\\|') %>% 
  map_chr(2)
genes = map_chr(fa, ~ attr(., 'Annot')) %>% 
  map_chr(~ gsub("^.*GN=", "", .) %>% gsub(" .*$", "", .)) %>% 
  setNames(accns)
freq_genes = genes[freq_prots] %>% na.omit() %>% unique()

# we also need a list of every phosphoprotein ever quantified
all_prots = phos$gene %>% 
  gsub("-[STY]-.*$", "", .) %>% 
  unique()
all_genes = genes[all_prots] %>% na.omit() %>% unique()

# bar chart: % of understudied proteins with a phosphosite vs interaction
net = readRDS("data/euk_interactomes/networks/network-species=Homo_sapiens-n_replicates=166-version=V2-n_features=1.rds")
g = net %>% 
  dplyr::select(protein_A, protein_B) %>% 
  igraph::graph_from_data_frame()
degrees = degree(g)
has_interaction = names(degrees) %>% 
  intersect(understudied$gene)
has_phosphorylation = freq_genes %>% 
  intersect(understudied$gene)
understudied = g2p %>%
  filter(n_pubs <= 10) %>% 
  distinct(gene) %>% 
  mutate(has_interaction = gene %in% has_interaction,
         has_phospho = gene %in% has_phosphorylation,
         both = has_interaction & has_phospho,
         color = ifelse(both,
                        'Both',
                        ifelse(has_interaction, 'Interaction',
                               ifelse(has_phospho, 'Regulatory phosphosite',
                                      'Neither'))) %>% 
           fct_relevel('Both', 
                       'Interaction',
                       'Regulatory phosphosite', 
                       'Neither'))
pal = c('Neither' = 'grey88',
        'Interaction' = pals::stepped2()[4],
        'Regulatory phosphosite' = pals::stepped2()[15],
        'Both' = 'black'
)
counts = dplyr::count(understudied, color)
labels = with(counts, paste0(color, ' (n = ', format(n, big.mark = ',') %>% 
                               trimws(), ')')) %>% 
  setNames(counts$color)
p = understudied %>% 
  ggplot(aes(x = '1', fill = color)) +
  geom_bar(size = 0.15, color = 'grey20') + 
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 2.5e3),
               color = 'grey50', size = 0.4) +
  scale_y_continuous('Understudied proteins', expand = c(0, 200),
                     limits = c(0, 2800),
                     breaks = seq(0, 2500, 500)) +
  scale_fill_manual('', values = pal, breaks = names(pal), labels = labels) +
  guides(fill = guide_legend(ncol = 1)) +
  coord_flip() +
  boxed_theme() +
  theme(panel.border = element_blank(),
        aspect.ratio = 0.2,
        legend.key.size = unit(0.3, 'lines'),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
p
ggsave("fig/understudied/interaction-vs-phospho-bar-chart.pdf", p9,
       width = 3.75, height = 3, units = "cm", useDingbats = FALSE)
