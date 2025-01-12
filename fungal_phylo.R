library(ggtree)
library(tidyverse)

##ONYGENALES##
Ony2 = read.newick("Ony_outgroup.trim.contree")
Ony2 = ape::root(Ony2, outgroup = "NR_121481.1") #Aspergillus fumigatus
plt_Ony2 = ggtree(Ony2, ladderize = T)

plt_Ony2$data = plt_Ony2$data %>%
  mutate(bootstrap = ifelse(isTip, 11, as.numeric(label))) #lowest bootstrap value
plt_Ony2 = plt_Ony2 +
  aes(linewidth = bootstrap) +
  scale_linewidth_continuous(range = c(0.25, 1)) +
  geom_treescale() +
  geom_tiplab(size = 3)

ggsave(filename = "crane_Ony2.pdf",
       plot = plt_Ony2,
       width = 11,
       height = 13,
       unit = "in", 
       device = "pdf")

##ASPERGILLACEAE##
Asp2 = read.newick("Aspergillus_outgroup_v4.trim.contree")
Asp2 = ape::root(Asp2, outgroup = "NR_157446.1") #Coccidioides immitis
plt_Asp2 = ggtree(Asp2, ladderize = T)

plt_Asp2$data = plt_Asp2$data %>%
  mutate(bootstrap = ifelse(isTip, 11, as.numeric(label))) 
plt_Asp2 = plt_Asp2 +
  aes(linewidth = bootstrap) +
  scale_linewidth_continuous(range = c(0.25, 1)) +
  geom_treescale() +
  geom_tiplab(size = 3)

ggsave(filename = "crane_Asp2.pdf",
       plot = plt_Asp2,
       width = 11,
       height = 13,
       unit = "in", 
       device = "pdf")

##TREMELLOMYCETES/CRYPTOCOCCUS-LIKE##
Trem = read.newick("Cryptococcus_outgroup_v2.trim.contree")
Trem = ape::root(Trem, outgroup = "MH858320.1") #Mycosarcoma (Ustilago) maydis
plt_Trem = ggtree(Trem, ladderize = T)

plt_Trem$data = plt_Trem$data %>%
  mutate(bootstrap = ifelse(isTip, 11, as.numeric(label))) 
plt_Trem = plt_Trem +
  aes(linewidth = bootstrap) +
  scale_linewidth_continuous(range = c(0.25, 1)) +
  geom_treescale() +
  geom_tiplab(size = 3)

ggsave(filename = "crane_Cryptococcus-like.pdf",
       plot = plt_Trem,
       width = 11,
       height = 13,
       unit = "in", 
       device = "pdf")

##SACCHAROMYCETALES##
Sacch2 = read.newick("Sacch_outgroup_v2.trim.contree")
Sacch2 = ape::root(Sacch2, outgroup = "MH860307.1") #Neurospora crassa CBS 709.71
plt_Sacch2 = ggtree(Sacch2, ladderize = T)

plt_Sacch2$data = plt_Sacch2$data %>%
  mutate(bootstrap = ifelse(isTip, 11, as.numeric(label))) 
plt_Sacch2 = plt_Sacch2 +
  aes(linewidth = bootstrap) +
  scale_linewidth_continuous(range = c(0.25, 1)) +
  geom_treescale() +
  geom_tiplab(size = 3)

ggsave(filename = "crane_Sacch2.pdf",
       plot = plt_Sacch2,
       width = 11,
       height = 13,
       unit = "in", 
       device = "pdf")
