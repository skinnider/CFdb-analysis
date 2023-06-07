# Are pY sites detected by the Beltrao group meta-analysis enriched for 
# protein complexes? 
setwd("~/git/CFdb-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(readxl)
source("R/theme.R")

# read Ochoa dataset
phos = read_excel("data/resources/Ochoa/41587_2019_344_MOESM4_ESM.xlsx")

# get universe
universe = unique(phos$uniprot)
# filter to pY
py = filter(phos, residue == 'Y') %$% unique(uniprot)

# read CORUM proteins
corum = read.delim("data/complex/CORUM/complexes_human.txt")

# map to uniprot accessions
map = read_tsv("data/identifier/HUMAN_9606_idmapping.dat.gz",
               col_names = c('uniprot', 'db', 'id'))
complex_proteins = map %>% 
  filter(db == 'Gene_Name', id %in% corum$gene_name) %$%
  unique(uniprot)

# test for enrichment
chisq.test(universe %in% complex_proteins, universe %in% py)

# plot
df = data.frame(protein = universe) %>% 
  mutate(pY = protein %in% py,
         complex = protein %in% complex_proteins)
freq = df %>% 
  group_by(pY) %>% 
  summarise(mean = mean(protein %in% complex_proteins),
            sd = sqrt(mean * (1 - mean) / n()))
pal = pals::stepped3()[c(1, 5) + 2] %>% c(pals::kelly()[1])
p = freq %>% 
  ggplot(aes(x = pY, fill = pY, y = mean)) +
  geom_col(width = 0.6, color = 'black', size = 0.25) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, size = 0.25) +
  scale_y_continuous('% of proteins', expand = expansion(c(0, 0.2)),
                     labels = ~ . * 100) +
  scale_x_discrete(limits = c(TRUE, FALSE), labels = c('pY', 'pS or pT only')) +
  scale_fill_manual(values = c('TRUE' = 'grey50', 'FALSE' = 'grey88')) +
  boxed_theme() +
  theme(aspect.ratio = 1.5,
        legend.position = 'none',
        axis.title.x = element_blank())
p
ggsave("fig/phosphorylation/pY-enrichment-complexes-Beltrao.pdf", p,
       width = 4.5, height = 4, units = "cm", useDingbats = FALSE)
