#if you need to reinstall phyloseq
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
#BiocManager::install("phyloseq")
#devtools::install_github('schuyler-smith/phylosmith')
#remotes::install_github("taowenmicro/ggClusterNet")
#devtools::install_github("jfq3/QsRutils", build_vignettes = TRUE)
#BiocManager::install("microbiome")

# Load Packages, Themes, etc... ####
pacman::p_load(phyloseq,
BiodiversityR,
tidyverse,
plyr,
geosphere,
phylosmith,
ggClusterNet,
microbiome,
multtest,
data.table, 
indicspecies,
bipartite,
simplecolors, 
ggnetwork,
igraph,
metagMisc,
QsRutils,
plotly,
FSA,
microbiome, 
DirichletMultinomial,
magrittr,
reshape2,
ggalluvial, 
wesanderson, 
tmap,
alluvial)

source("~/Desktop/ggplot_themes/ggplot_themes.R")
source("~/Desktop/R_color_palettes/Gadek_custom_colors.R")

theme_set(theme_arial_clean())

#Set constants ####
merged = F #to run merged analysis where all species are merged or separate
mergeby <- "sci_name"

#Read in physeq object processed by Phyloseq_wrangler.R file####
load("data/physeq_all.Rdata")
load("data/physeq_rarefied.Rdata")
load("data/physeq_merged_spp.Rdata")
load("data/physeq_rarefied_merged_spp.Rdata")

#read in processed sample data  and taxonomy file
all_meta <- read_csv("data/Bird_predictors.csv")
taxonomy <- read.table("data/zotus_tax_v5.txt", header = TRUE, sep = '\t', row.names = 1)
#remove white space
taxonomy <-taxonomy %>% 
  mutate(across(where(is.character), str_trim))

# Estimate sample coverage ####
phyloseq_coverage(physeq, correct_singletons = FALSE, add_attr = T)%>%
  dplyr::summarise(mean = mean(SampleCoverage, na.rm=T),
                   sd = sd(SampleCoverage, na.rm=T),
                   min = min(SampleCoverage, na.rm=T),
                   max = max(SampleCoverage, na.rm=T))

# calculate Good's coverage (only uses singletons)
#summary(goods(otu_table(physeq)))

# calculate extrapolated species richness and plot against observed
estimate_richness(physeq, measures=c('Observed', 'Chao1'))%>%
ggplot( aes(x=Observed, y=Chao1)) +
  geom_point(color="gray55")+
  geom_abline(slope=1, intercept = 0)+
  theme(aspect.ratio = 1)

# check whether extrapolated species richness appears to be sensitive to sequencing depth
cbind(estimate_richness(physeq, measures=c('Observed', 'Chao1')), readcounts=sample_sums(physeq))%>%
  ggplot( aes(x=readcounts, y=Chao1)) +
  geom_point(color="gray55")+
  geom_smooth(method="lm")

#linear model to get significance stats, probably violating assumptions here
summary(lm(Chao1 ~ readcounts, data = cbind(estimate_richness(physeq, measures=c('Observed', 'Chao1')), readcounts=sample_sums(physeq))))

#extrapolated species richness appears to be sensitive to readcounts, but that's okay because using rarefied dataset to conduct community analyses

#Get summary statistics on samples sizes ####
all_meta%>%
  group_by(Order1)%>%
  dplyr::summarise(n=n())

all_meta%>%
  group_by(Order1, sci_name)%>%
  dplyr::summarise(n=n())%>%
  #View()%>%
  write.csv(., "data/taxon_table.csv")

all_meta%>%
  group_by(sci_name)%>%
  #View()%>%
  dplyr::summarise(n=n())

# Confirm physeq object is correct NOTE this is before rarefaction
ntaxa(physeq) #526 zOTUs in samples
phyloseq::nsamples(physeq) #183 seq but 167 once <2000 reads were removed and 160 if we remove suspicious cranes were removed
sample_names(physeq)
taxa_names(physeq)
rank_names(physeq)
sample_variables(physeq)

##Animal symbiont summary stats####

anm.sym<-read_csv("data/animal_symb_FUNGuild_v2.csv")

anm.sym%>%
  dplyr::select(-c(1))%>%
  filter(zotu %in% rownames(tax_table(physeq)))%>%
  group_by(animal.trophic)%>%
  dplyr::summarise(n =n()/526)


anm.s<-psmelt(physeq)%>%
  as.data.frame()%>%
  dplyr::select(OTU,  sci_name, Sample, family, Abundance)%>%
  #dplyr::rename(sci_name = Sample)%>%
  left_join(anm.sym%>%dplyr::rename(OTU=zotu))%>%
  filter(Abundance >0)%>%
  group_by(Sample)%>%
  dplyr::mutate(n_total = n())%>%
  group_by(Sample, animal.trophic)%>%
  mutate(n_trophic = n(),
         perc_trophic = n_trophic/n_total)


t.test(anm.s%>%filter(animal.trophic=="y")%>%pull(Abundance), anm.s%>%filter(animal.trophic=="n")%>%pull(Abundance))

###Cranes

psmelt(physeq)%>%
  as.data.frame()%>%
  filter(bird == "crane",
         Abundance >0)%>%
  dplyr::select(OTU,  sci_name, bird)%>%
  left_join(anm.sym%>%dplyr::rename(OTU=zotu))%>%
  distinct()%>%
  group_by(sci_name)%>%
  mutate(n_total=n())%>%
  group_by(sci_name, animal.trophic)%>%
  mutate(n_trophic =n(),
         perc_trophic =n_trophic/n_total)%>%
  dplyr::select(sci_name, perc_trophic)%>%
  distinct()




#Core microbiome####

# Calculate compositional version of the data
# (relative abundances)
pseq.rel <- microbiome::transform(physeq, "compositional")

#calculate core microbiome fungi that wer present in >50% of samples and greater than 0.01 relative abundance

pseq.core <- microbiome::core(pseq.rel, detection = 0.5, prevalence = 0.01) #make phyloseq object
core.taxa <- taxa_names(pseq.core)
core.taxa

# Composition plot
d<-psmelt(physeq)%>%
  #filter(OTU %in% core.taxa)%>%
  mutate(all_reads=sum(Abundance))%>%
  group_by(family, all_reads)%>%
  dplyr::summarise(reads.family = sum(Abundance),
                   perc_family = reads.family/all_reads)%>%
  filter(family !="")%>%
  distinct()%>%
  arrange(perc_family)

ggplot(d, aes(x=all_reads, y=perc_family, fill=family))+
  geom_bar(color="black", stat="identity", position="stack")+
  labs(x=NULL, y="Percentage of core reads")+
  theme(aspect.ratio = 1.5,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Assess host-fungi specifics through summary stats####
#Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
f.summs = filter_taxa(physeq, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

f.summs.mat<- 
  physeq@otu_table%>%
  as.data.frame()%>%
  mutate_if(is.character, str_trim)%>%
  rownames_to_column(var="guid")%>%
  pivot_longer(cols =  contains("otu"), names_to = c("OTU"), values_to = "Abundance") %>% #make it long form
  filter(Abundance > 0) %>%
  left_join(taxonomy%>%rownames_to_column(var="OTU"), by = "OTU") %>%
  mutate(order = replace_na(order, "Unknown"),
         order = if_else(order=="", "Unknown", order),
         order = if_else(order==" ", "Unknown", order)) %>%
  group_by(guid, order) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup() %>%
  pivot_wider(names_from = c(order), values_from = Abundance, values_fill = list(Abundance = 0))%>%
  left_join(all_meta%>%dplyr::select(guid, sci_name, Family1), by="guid")
  

physeq@otu_table%>%
 as.data.frame()%>%
  rownames_to_column(var="guid")%>%
  pivot_longer(2:length(.))%>%
  left_join(taxonomy%>%rownames_to_column(var="name"))%>%
  mutate_if(is.character, str_trim)%>%
  left_join(all_meta%>%dplyr::select(Family1, sci_name, guid))%>%
  group_by(class)%>%
  mutate(sumorder = sum(value))%>%
  ungroup()%>%
  filter(value >0)%>%
  group_by(Family1,sci_name,  class)%>%
  dplyr::summarise(n=n(),
    sum =sum(value))%>%
  filter(Family1 =="Gruidae",
         !class %in% c("Atractiellomycetes", "Exobasidiomycetes", "Ustilaginomycotina_cls_Incertae_sedis", ""))%>%
  mutate(class = factor(class, levels = c("Dothideomycetes", 
                                          "Eurotiomycetes",
                                          "Sordariomycetes",
                                          "Saccharomycetes", 
                                          "Pezizomycetes", 
                                          "Leotiomycetes",
                                          "Orbiliomycetes",
                                          "Arthoniomycetes",
                                          "Lecanoromycetes",
                                          "Tremellomycetes",
                                          "Malasseziomycetes",
                                          "Agaricomycetes",
                                          "Microbotryomycetes",
                                          "Agaricostilbomycetes",
                                          "Wallemiomycetes",
                                          "Ustilaginomycetes",
                                          "Tritirachiomycetes",
                                          "Mucoromycetes",
                                          "Mortierellomycetes",
                                          "Chytridiomycetes")))%>%
  ggplot(aes(x= class, y = log10(sum), fill = sci_name))+
  geom_col(position = position_dodge())+
  scale_fill_manual(values=c("grey33", "grey65"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = c(0.8, 0.9))+
  scale_fill_manual(values =c("#803B31", "#B2AAA2"))+
  labs(y = "zOTUs", x = NULL)

f.summs.mat%>%mutate(sum = rowSums(across(where(is.numeric)), na.rm=TRUE))%>%
  dplyr::select(guid, Family1, sum)%>%
  group_by(Family1)%>%
  mutate(n = n())%>%
  ungroup()%>%
  group_by(Family1, n)%>%
  dplyr::summarise(sum = sum(sum))%>%
  arrange(.,sum)%>%
  #filter(Family1 !="Gruidae")%>%
  ggplot2::ggplot(.,aes(y= sum, x = reorder(Family1, -sum), label = n))+
  geom_col()+
  geom_text(vjust = -0.5)+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), limits = c(0, 400000))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = NULL, y = "Reads")

#Rarified zOTUs ####

ntaxa(physeq.r) #508

#Find taxa unique to cranes and other birds####

orders<-taxonomy%>%
  filter(order != "")%>%
  select(order)%>%
  distinct()%>%
  pull()

N <- length(orders)

##Unique to cranes####

unique.to.cranes <- data.frame(order   = character(N),
                               unique = character(N),
                               number_zotu = numeric(N),
                               number_reads = numeric(N))

for(i in orders){
  
s<-physeq.r@otu_table%>%
  as.data.frame()%>%
  rownames_to_column("guid")%>%
  mutate(guid = gsub("^([0-9]+)", "S\\1", guid))%>%
  pivot_longer(cols=2:509, names_to = "zOTU", values_to = "reads")%>%
  left_join(., all_meta%>%dplyr::select(guid, sci_name, bird))%>%
  left_join(., taxonomy%>%rownames_to_column("zOTU"))%>%
  filter(order !="",
         reads >0, 
         order == i)

crane_not <-unique(s$bird)

if(length(crane_not)==1 && crane_not =="crane"){
  unique <- "unique"
}else{
  unique <- "not unique"
}

number_zotu <- s%>%
  select(zOTU)%>%
  distinct()%>%
  pull()%>%
  length()

number_zotu <- s%>%
  select(zOTU)%>%
  distinct()%>%
  pull()%>%
  length()

number_reads <-s%>%
  ungroup()%>%
  select(reads)%>%
  dplyr::summarise(sum = sum(reads))%>%
  pull()

unique.to.cranes[i, "order"] <- i
unique.to.cranes[i, "unique"] <- unique
unique.to.cranes[i, "number_zotu"] <- number_zotu
unique.to.cranes[i, "number_reads"] <- number_reads

unique.to.cranes <- unique.to.cranes%>%
  filter(number_reads>0)

}

##Unique to other birds####

unique.to.others <- data.frame(order   = character(N),
                               unique = character(N),
                               number_zotu = numeric(N),
                               number_reads = numeric(N))

for(i in orders){
  
  s<-physeq.r@otu_table%>%
    as.data.frame()%>%
    rownames_to_column("guid")%>%
    mutate(guid = gsub("^([0-9]+)", "S\\1", guid))%>%
    pivot_longer(cols=2:509, names_to = "zOTU", values_to = "reads")%>%
    left_join(., all_meta%>%dplyr::select(guid, sci_name, bird))%>%
    left_join(., taxonomy%>%rownames_to_column("zOTU"))%>%
    filter(order !="",
           reads >0, 
           order == i)
  
  crane_not <-unique(s$bird)
  
  if(length(crane_not)==1 && crane_not =="other"){
    unique <- "unique"
  }else{
    unique <- "not unique"
  }
  
  number_zotu <- s%>%
    select(zOTU)%>%
    distinct()%>%
    pull()%>%
    length()
  
  number_zotu <- s%>%
    select(zOTU)%>%
    distinct()%>%
    pull()%>%
    length()
  
  number_reads <-s%>%
    ungroup()%>%
    select(reads)%>%
    dplyr::summarise(sum = sum(reads))%>%
    pull()
  
  unique.to.others[i, "order"] <- i
  unique.to.others[i, "unique"] <- unique
  unique.to.others[i, "number_zotu"] <- number_zotu
  unique.to.others[i, "number_reads"] <- number_reads
  
  unique.to.others <- unique.to.others%>%
    filter(number_reads>0)
  
}

#Distance Decay####
library(patchwork)
library(vegan)

#vegan object
veg <- read.table("data/all_lungs_v2.txt", header = TRUE, sep = '\t') %>%
  mutate(guid = paste0("S", guid))%>%
  left_join(as.data.frame(physeq.r@otu_table)%>%rownames_to_column("guid"), by = "guid") %>%
  drop_na()
veg_meta <-veg[,1:13] #metadata columns
veg_otus <-veg[,14:521] #otus
veg_rel<-decostand(veg_otus, "total")

#non-cranes
veg.nc <- read.table("data/all_lungs_v2.txt", header = TRUE, sep = '\t') %>%
  mutate(guid = paste0("S", guid))%>%
  filter(bird == "other")%>%
  left_join(as.data.frame(physeq.r@otu_table)%>%rownames_to_column("guid"), by = "guid") %>%
  drop_na()
veg_meta.nc <-veg.nc[,1:13] #metadata columns
veg_otus.nc <-veg.nc[,14:521] #otus
veg_rel.nc<-decostand(veg_otus.nc, "total")

# Spatial Diversity
## Mantel test
### All birds
latlon <-data.frame(veg_meta$long, veg_meta$lat) #extract geographic coordinates
spat<-distm(latlon) # create geographic distance matrix
dismat<-vegan::vegdist(vegan::decostand(veg_otus, "total"), "bray")
mc<-mantel.correlog(dismat, as.dist(spat, upper=FALSE, diag=FALSE),cutoff = FALSE)  #mantel correlogram
plot(mc)
mantel(dismat, as.dist(spat, upper=FALSE, diag=FALSE))
#r=0.03357  p=0.156 not significant
mc #describe stats
###samples collected within 58 km tend to be more similar to each other

## Distance decay
plot(as.dist(spat,upper=FALSE, diag=FALSE), (1-dismat), xlab = "Distance (m)", ylab = "Bray-Curtis Similarity")
mod<-lm(as.vector((1-dismat))~as.vector(as.dist(spat,upper=FALSE, diag=FALSE)))
mod #y-int =6.264e-02 m = -2.981e-08  (very little slope)
abline(mod, col="#b21f80", lwd=3)
###very low decline in similarity of two samples as the distance between them increases 

patchwork::mc + mod + plot_layout(nrow = 1)


# Spatial Diversity
## Mantel test
### Non-crane birds
latlon <-data.frame(veg_meta.nc$long, veg_meta.nc$lat) #extract geographic coordinates
spat<-distm(latlon) # create geographic distance matrix
dismat<-vegan::vegdist(vegan::decostand(veg_otus.nc, "total"), "bray")
mc<-mantel.correlog(dismat, as.dist(spat, upper=FALSE, diag=FALSE),cutoff = FALSE)  #mantel correlogram
plot(mc)
mantel(dismat, as.dist(spat, upper=FALSE, diag=FALSE))
#r=0.03357  p=0.156 not significant
mc #describe stats
###samples collected within 58 km tend to be more similar to each other

## Distance decay
plot(as.dist(spat,upper=FALSE, diag=FALSE), (1-dismat), xlab = "Distance (m)", ylab = "Bray-Curtis Similarity")
mod<-lm(as.vector((1-dismat))~as.vector(as.dist(spat,upper=FALSE, diag=FALSE)))
mod #y-int =6.264e-02 m = -2.981e-08  (very little slope)
abline(mod, col="#b21f80", lwd=3)
###very low decline in similarity of two samples as the distance between them increases 

patchwork::mc + mod + plot_layout(nrow = 1)

# Alpha Diversity #### 
#NOTE according to function documentation this must be used with untrimmed data
er<-estimate_richness(physeq.r, measures=c("Observed", "InvSimpson", "Shannon", "Chao1", "ACE"))%>%
  rownames_to_column(var = "guid")%>%
  left_join(all_meta)

##Summary stats####
#mean and sd
er%>%
  dplyr::summarise(mean.r = mean(Observed, na.rm=T),
                   sd.r = sd(Observed, na.rm=T),
                   mean.sh = mean(Shannon, na.rm=T),
                   sd.sh = sd(Shannon, na.rm=T),
                   mean.is = mean(InvSimpson, na.rm=T),
                   sd.is = sd(InvSimpson, na.rm=T),
                   mean.ch = mean(Chao1, na.rm=T),
                   sd.ch = sd(Chao1, na.rm=T))

#ranges
er%>%
  dplyr::summarise(min.r = min(Observed, na.rm=T),
                   max.r = max(Observed, na.rm=T),
                   min.sh = min(Shannon, na.rm=T),
                   max.sh = max(Shannon, na.rm=T),
                   min.is = min(InvSimpson, na.rm=T),
                   max.is = max(InvSimpson, na.rm=T),
                   min.ch = min(Chao1, na.rm=T),
                   max.ch = max(Chao1, na.rm=T))

##Box plots####
p <-er%>%
  ggplot(aes(x=Order1, y= Chao1))+
  geom_boxplot(aes(), fill="grey")+
  theme(axis.text.x = element_text(angle=45, vjust=1.1, hjust = 1, size=10),
        aspect.ratio = 0.7)+
  labs(x=NULL)
ggsave(p, file = "figures/Chao1_Order.pdf", width=5, height=3, device = cairo_pdf)

p <-er%>%
  ggplot(aes(x=Order1, y= Shannon))+
  geom_boxplot(aes(), fill="grey")+
  theme(axis.text.x = element_text(angle=45, vjust=1.1, hjust = 1, size=10),
        aspect.ratio = 0.7)+
  labs(x=NULL)
p
ggsave(p, file = "figures/Shannnon_Order.pdf", width=5, height=3, device = cairo_pdf)

p <-er%>%
  ggplot(aes(x=Order1, y= Observed))+
  geom_boxplot(aes(), fill="grey")+
  theme(axis.text.x = element_text(angle=45, vjust=1.1, hjust = 1, size=10),
        aspect.ratio = 0.7)+
  labs(x=NULL)
p
ggsave(p, file = "figures/Observed_Order.pdf", width=5, height=3, device = cairo_pdf)

p <-er%>%
  ggplot(aes(x=Order1, y= InvSimpson))+
  geom_boxplot(aes(), fill="grey")+
  theme(axis.text.x = element_text(angle=45, vjust=1.1, hjust = 1, size=10),
        aspect.ratio = 0.7)+
  labs(x=NULL)
p
ggsave(p, file = "figures/InvSimp_Order.pdf", width=5, height=3, device = cairo_pdf)

#Beta Diversity####
#Permanova/ANOSIM  ####
### All birds####
# For all ANOSIM and PERMANOVA we will use rarefied datasets

metadata <- as(sample_data(physeq.r), "data.frame")
metadata <- 
set.seed(1982)
vegan::adonis2(phyloseq::distance(physeq.r, method="bray") ~  + Order1,
               by = "margin", data = metadata, permutations = 10000)
#          Df SumOfSqs      R2     F   Pr(>F)   
# Order1    11    5.635 0.07582 1.156 0.005099 **
# Residual 155   68.688 0.92418                  
# Total    166   74.323 1.00000   

# Calculate multivariate dispersions
#check for homogeneity of variances between groups with betadisp
set.seed(1982)
permutest(betadisper(phyloseq::distance(physeq.r, method="bray"), physeq.r@sam_data$Order1),  permutations = 10000)

# betadispersion significant!! -> variances different USE ANOSIM for Order

set.seed(1982)
anosim(phyloseq::distance(physeq.r, method="bray"),  physeq.r@sam_data$Order1, permutations = 10000)
# ANOSIM statistic R: 0.1648 
# Significance: 0.0015998 

set.seed(1982)
vegan::adonis2(phyloseq::distance(physeq.r, method="bray") ~  + Migratory.status,
       by = "margin", data = metadata, permutations = 10000)
# Df SumOfSqs      R2     F Pr(>F)   
# Migratory.status   2    1.450 0.01951 1.632 0.0011 **
# Residual         164   72.873 0.98049                
# Total            166   74.323 1.00000                     

# Calculate multivariate dispersions
#check for homogeneity of variances between groups with betadisp
set.seed(1982)
permutest(betadisper(phyloseq::distance(physeq.r, method="bray"), physeq.r@sam_data$Migratory.status),  permutations = 10000)

# betadispersion significant!! -> variances different USE ANOSIM for Migration

set.seed(1982)
anosim(phyloseq::distance(physeq.r, method="bray"),  physeq.r@sam_data$Migratory.status, permutations = 10000)
#ANOSIM statistic R: 0.1753 
#Significance: 0.0015998 

#looks good with ANOSIM so ok

set.seed(1982)
vegan::adonis2(phyloseq::distance(physeq.r, method="bray") ~  + Trophic.Level,
               by = "margin", data = metadata, permutations = 10000)

#                 Df SumOfSqs      R2      F    Pr(>F)    
# Trophic.Level   2    1.716 0.02308 1.9377 9.999e-05 ***
# Residual      164   72.608 0.97692                     
# Total         166   74.323 1.00000  

# Trophic level highly significant

## Calculate multivariate dispersion
set.seed(1982)
permutest(betadisper(phyloseq::distance(physeq.r, method="bray"), physeq.r@sam_data$Trophic.Level), permutations = 10000)

# betadispersion significant!! -> variances different USE ANOSIM for Trophic.Level

set.seed(1982)
anosim(phyloseq::distance(physeq.r, method="bray"),  physeq.r@sam_data$Trophic.Level, permutations = 10000)

# ANOSIM statistic R: 0.1129 
# Significance: 0.018798

#Year
set.seed(1982)
vegan::adonis2(phyloseq::distance(physeq.r, method="bray") ~  + year,
               by = "margin", data = metadata, permutations = 10000)

#             Df SumOfSqs      R2      F Pr(>F)    
# year       1    0.997 0.01342 2.2444  5e-04 ***
# Residual 165   73.326 0.98658                  
# Total    166   74.323 1.00000    

# year highly significant

## Calculate multivariate dispersion
set.seed(1982)
permutest(betadisper(phyloseq::distance(physeq.r, method="bray"), physeq.r@sam_data$year), permutations = 10000)

# betadispersion significant!! -> variances different USE ANOSIM for Trophic.Level

set.seed(1982)
anosim(phyloseq::distance(physeq.r, method="bray"),  physeq.r@sam_data$year, permutations = 10000)

#ANOSIM statistic R: 0.1745 
#Significance: 9.999e-05 
##Highly significant for year!!


### Just Cranes ####
metadata <- as(sample_data(subset_samples(physeq.r, bird == "crane")), "data.frame")

#subspecies
sus.nk <-c(281742, 281737, 281745, 281794, 281809, 281811, 281821)

metadata <- as(sample_data(subset_samples(physeq.r, c(bird == "crane", !NK %in% sus.nk))), "data.frame")

phseq.crane.sub <- subset_samples(physeq.r, c(bird == "crane"))

set.seed(1982)
vegan::adonis2(phyloseq::distance(phseq.crane.sub, method="bray") ~  + sci_name,
               by = "margin", data = metadata, permutations = 10000)

#             Df SumOfSqs      R2      F Pr(>F)
# sci_name   1    0.544 0.01015 1.2307 0.1316
# Residual 120   53.049 0.98985              
# Total    121   53.593 1.00000    

## Calculate multivariate dispersions
permutest(betadisper(phyloseq::distance(phseq.crane.sub, method="bray"), phseq.crane.sub@sam_data$sci_name),  permutations = 10000)

#year
metadata <- as(sample_data(subset_samples(physeq.r, c(bird == "crane", sex %in% c("male", "female")))), "data.frame")%>%
  filter(!is.na(sex))

set.seed(1982)
vegan::adonis2(phyloseq::distance(phseq.crane.sex, method="bray") ~  + factor(year),
               by = "margin", data = metadata, permutations = 10000)

#                 Df SumOfSqs      R2      F Pr(>F)    
# factor(year)   1    1.063 0.02195 2.4465  2e-04 ***
# Residual     109   47.380 0.97805                  
# Total        110   48.444 1.00000  

## Calculate multivariate dispersions
permutest(betadisper(phyloseq::distance(phseq.crane.sex, method="bray"), phseq.crane.sex@sam_data$year),  permutations = 10000)
#             Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups      1 0.009626 0.0096257 4.0814  10000 0.0469 *
# Residuals 109 0.257069 0.0023584    

#not significant use adonis
set.seed(1982)
anosim(phyloseq::distance(phseq.crane.sex, method="bray"),  subset_samples(phseq.crane.sex)@sam_data$year, permutations = 10000)

# ANOSIM statistic R: 0.1058 
# Significance: 9.999e-05 

### Non-cranes ####

metadata <- as(sample_data(subset_samples(physeq.r, bird == "other")), "data.frame")

set.seed(1982)
vegan::adonis2(phyloseq::distance(subset_samples(physeq.r, bird == "other"), method="bray") ~  + Order1,
               by = "margin", data = metadata, permutations = 10000)

#           Df SumOfSqs     R2      F Pr(>F)
# Order1   10   4.4658 0.2283 1.0059 0.4338
# Residual 34  15.0948 0.7717              
# Total    44  19.5606 1.0000                  

## Calculate multivariate dispersions
permutest(betadisper(phyloseq::distance(subset_samples(physeq.r, bird == "other"), method="bray"), subset_samples(physeq.r, bird == "other")@sam_data$Order1),  permutations = 10000)
#           Df  Sum Sq  Mean Sq      F N.Perm    Pr(>F)    
# Groups    10 1.45186 0.145186 27.465  10000 9.999e-05 ***
# Residuals 34 0.17973 0.005286   

#Very Significant use ANOSIM

set.seed(1982)
anosim(phyloseq::distance(subset_samples(physeq.r, bird == "other"), method="bray"),  subset_samples(physeq.r, bird == "other")@sam_data$Order1, permutations = 10000)

#ANOSIM statistic R: -0.03049 
#Significance: 0.63834

# betadispersion not significant!!
set.seed(1982)
vegan::adonis2(phyloseq::distance(subset_samples(physeq.r, bird == "other"), method="bray") ~  + Trophic.Level,
               by = "margin", data = metadata, permutations = 10000)

#               Df SumOfSqs      R2      F Pr(>F)  
# Trophic.Level  2   1.2082 0.06177 1.3825  0.027 *
# Residual      42  18.3523 0.93823                
# Total         44  19.5606 1.00000        

## Calculate multivariate dispersions
permutest(betadisper(phyloseq::distance(subset_samples(physeq.r, bird == "other"), method="bray"), subset_samples(physeq.r, bird == "other")@sam_data$Trophic.Level),  permutations = 10000)

# betadispersion significant!! -> variances different USE ANOSIM for Migration

set.seed(1982)
anosim(phyloseq::distance(subset_samples(physeq.r, bird == "other"), method="bray"),  subset_samples(physeq.r, bird == "other")@sam_data$Trophic.Level, permutations = 10000)
#ANOSIM statistic R: 0.2392 
#Significance: 0.00039996 

#Migration
set.seed(1982)
vegan::adonis2(phyloseq::distance(subset_samples(physeq.r, bird == "other"), method="bray") ~  + Migratory.status,
               by = "margin", data = metadata, permutations = 10000)

#                   Df SumOfSqs     R2      F Pr(>F)  
# Migratory.status  2   1.1755 0.0601 1.3427 0.0416 *
# Residual         42  18.3850 0.9399                
# Total            44  19.5606 1.0000             

## Calculate multivariate dispersions
permutest(betadisper(phyloseq::distance(subset_samples(physeq.r, bird == "other"), method="bray"), subset_samples(physeq.r, bird == "other")@sam_data$Migratory.status),  permutations = 10000)

# Not significant use PERMANOVA


#year
set.seed(1982)
vegan::adonis2(phyloseq::distance(subset_samples(physeq.r, bird == "other"), method="bray") ~  + year,
               by = "margin", data = metadata, permutations = 10000)

#           Df SumOfSqs      R2      F Pr(>F)
# year      1   0.4307 0.02202 0.9682  0.507
# Residual 43  19.1298 0.97798              
# Total    44  19.5606 1.00000                 

## Calculate multivariate dispersions
permutest(betadisper(phyloseq::distance(subset_samples(physeq.r, bird == "other"), method="bray"), subset_samples(physeq.r, bird == "other")@sam_data$year),  permutations = 10000)
#           Df  Sum Sq  Mean Sq      F N.Perm    Pr(>F)    
# Groups    10 1.41168 0.141168 22.042  10000 9.999e-05 ***
# Residuals 34 0.21775 0.006405   

#Very Significant use ANOSIM

set.seed(1982)
anosim(phyloseq::distance(subset_samples(physeq.r, bird == "other"), method="bray"),  subset_samples(physeq.r, bird == "other")@sam_data$year, permutations = 10000)

#ANOSIM statistic R: -0.06973 
#Significance: 0.90541 

### Pathogens only####

metadata <- as(sample_data(anm.sym.physeq.r), "data.frame")

set.seed(1982)
vegan::adonis2(phyloseq::distance(anm.sym.physeq.r, method="bray") ~  + Order1,
               by = "margin", data = metadata, permutations = 10000)

#           Df SumOfSqs      R2      F Pr(>F)  
# Order1    11    5.292 0.07677 1.1717 0.0249 *
# Residual 155   63.642 0.92323                
# Total    166   68.934 1.00000            

## Calculate multivariate dispersions
permutest(betadisper(phyloseq::distance(anm.sym.physeq.r, method="bray"), anm.sym.physeq.r@sam_data$Order1),  permutations = 10000)

#           Df Sum Sq  Mean Sq  F     N.Perm    Pr(>F)    
# Groups     11 1.8906 0.171872 11.616  10000 9.999e-05 ***
# Residuals 155 2.2934 0.014796  

# betadispersion significant!! -> variances different USE ANOSIM for Order

set.seed(1982)
anosim(phyloseq::distance(anm.sym.physeq.r, method="bray"), anm.sym.physeq.r@sam_data$Order1, permutations = 10000)
#ANOSIM statistic R: 0.1112 
#Significance: 0.020698 

#Trophic Level
set.seed(1982)
vegan::adonis2(phyloseq::distance(anm.sym.physeq.r, method="bray") ~  + Trophic.Level,
               by = "margin", data = metadata, permutations = 10000)

#                   Df SumOfSqs R2   F       Pr(>F)   
# Trophic.Level   2    1.582 0.02296 1.9266  0.002 **
# Residual      164   67.352 0.97704                 
# Total         166   68.934 1.00000    

## Calculate multivariate dispersions
permutest(betadisper(phyloseq::distance(anm.sym.physeq.r, method="bray"), anm.sym.physeq.r@sam_data$Trophic.Level),  permutations = 10000)

# betadispersion significant!! -> variances different USE ANOSIM for Migration

set.seed(1982)
anosim(phyloseq::distance(anm.sym.physeq.r, method="bray"), anm.sym.physeq.r@sam_data$Trophic.Level, permutations = 10000)
#ANOSIM statistic R: 0.0911 
#Significance: 0.046395 

set.seed(1982)
vegan::adonis2(phyloseq::distance(anm.sym.physeq.r, method="bray") ~  + Migratory.status,
               by = "margin", data = metadata, permutations = 10000)

#                   Df SumOfSqs      R2      F Pr(>F)   
# Migratory.status   2    1.432 0.02078 1.7401 0.0031 **
# Residual         164   67.502 0.97922                 
# Total            166   68.934 1.00000       

## Calculate multivariate dispersions
permutest(betadisper(phyloseq::distance(anm.sym.physeq.r, method="bray"), anm.sym.physeq.r@sam_data$Migratory.status),  permutations = 10000)

# betadispersion significant Use ANOSIM

set.seed(1982)
anosim(phyloseq::distance(anm.sym.physeq.r, method="bray"), anm.sym.physeq.r@sam_data$Migratory.status, permutations = 10000)
# ANOSIM statistic R: 0.06585 
# Significance: 0.14249 

#Year

set.seed(1982)
vegan::adonis2(phyloseq::distance(anm.sym.physeq.r, method="bray") ~  + year,
               by = "margin", data = metadata, permutations = 10000)

#           Df SumOfSqs      R2      F Pr(>F)   
# year       1    0.964 0.01398 2.3391 0.0016 **
# Residual 165   67.971 0.98602                 
# Total    166   68.934 1.00000         

## Calculate multivariate dispersions
permutest(betadisper(phyloseq::distance(anm.sym.physeq.r, method="bray"), anm.sym.physeq.r@sam_data$year),  permutations = 10000)

# betadispersion significant!! -> variances different USE ANOSIM for Order

set.seed(1982)
anosim(phyloseq::distance(anm.sym.physeq.r, method="bray"), anm.sym.physeq.r@sam_data$year, permutations = 10000)

#ANOSIM statistic R: 0.1243 
#Significance: 9.999e-05 

#Year non-cranes
metadata <- as(sample_data(subset_samples(physeq.r, bird == "other")), "data.frame")
set.seed(1982)
vegan::adonis2(phyloseq::distance(subset_samples(anm.sym.physeq.r, bird == "other"), method="bray") ~  + year,
               by = "margin", data = metadata, permutations = 10000)

#           Df SumOfSqs      R2     F Pr(>F)
# year      1    0.298 0.01665 0.728 0.7915
# Residual 43   17.602 0.98335             
# Total    44   17.900 1.00000            

## Calculate multivariate dispersions
permutest(betadisper(phyloseq::distance(subset_samples(anm.sym.physeq.r, bird == "other"), method="bray"), subset_samples(anm.sym.physeq.r, bird == "other")@sam_data$year),  permutations = 10000)

# significant use PERMANOVA
set.seed(1982)
anosim(phyloseq::distance(subset_samples(anm.sym.physeq.r, bird == "other"), method="bray"),  subset_samples(anm.sym.physeq.r, bird == "other")@sam_data$Trophic.Level, permutations = 10000)


#NMDS ####
##All birds####
psprn.ord <- ordinate(physeq.r, "NMDS", "bray")

all_meta_ordin <- as_tibble(all_meta)

nmds_points <- 
  as_tibble(psprn.ord$points,  rownames = "guid") %>%
  left_join(y = all_meta_ordin, by = "guid")


## Plot NMDS ####
{
  c25 <- c(
    "skyblue2", 
    # lt pink
    "maroon",
    "darkorange4",
    "palegreen3",
    "darkorange2",
    "purple3",
    # lt purple
    "grey66",
    "gold2",
    "burlywood2",
    # lt orange
    "steelblue4",
    "pink2",
    "bisque3",
    "yellow4"
    
    
  )
  
  hull <- nmds_points %>%
    group_by(Family1)%>%
    dplyr::mutate(n = dplyr::n())%>%
    filter(#!Family1 == "Gruidae",
      n >2)%>%
    ungroup()%>%
    group_by(Family1) %>%
    dplyr::slice(chull(MDS1, MDS2))
  
  p<-nmds_points %>%
    group_by(Family1)%>%
    dplyr::mutate(n = dplyr::n())%>%
    filter(#!Family1 == "Gruidae",
      n >2)%>%
    dplyr::select(Family1, MDS1, MDS2, n)%>%
    mutate(Family1=fct_relevel(Family1,c("Corvidae","Fringillidae", "Strigidae", "Caprimulgidae", "Accipitridae", "Turdidae", "Tytonidae", "Gruidae")))%>%
    
    ggplot(aes(x = MDS1, y = MDS2, color = Family1, fill=Family1)) + 
    geom_point(size = 1)+
    # geom_polygon(data = hull,
    #              alpha = 0.2,
    #              color = "transparent") +
    stat_ellipse(linetype=2)+
    scale_color_manual(values = c25) +
    scale_fill_manual(values = c25)+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme(aspect.ratio=1)
  
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_with_crane_family.pdf", sep = ""),
    width = 5,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
}

p


{
  
  
  p<-nmds_points %>%
    group_by(Order1)%>%
    dplyr::mutate(n = dplyr::n())%>%
    #filter(!Family1 == "Gruidae")%>%
    ggplot(aes(x = MDS1, y = MDS2, color = Trophic.Level, fill = Trophic.Level)) + 
    geom_point(size = 0.5)+
    # geom_polygon(data = hull,
    #              alpha = 0.2,
    #              color = "transparent") +
    stat_ellipse(linetype=2)+
    scale_color_manual(values=c("red4","#018571", "#CD9600"))+
    scale_fill_manual(values=c("red4","#018571", "#CD9600"))+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme(aspect.ratio=1)
  
  
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_with_crane_trophic_level_ellispe.pdf", sep = ""),
    width = 5,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
}

p

{
  
  cmig <- c(
    
    # lt purple
    "#077551",
    "#362392",
    "#AF700B",
    "#AF370B"
    
  )
 
  p<-nmds_points %>%
    group_by(Order1)%>%
    dplyr::mutate(n = dplyr::n(),
                 Migration = factor(Migration, levels = c(3,2,1)) )%>%
    #filter(!Family1 == "Gruidae")%>%
    ggplot(aes(x = MDS1, y = MDS2, color = as.factor(Migration), fill = as.factor(Migration))) + 
    geom_point(size = 0.5)+
    # geom_polygon(data = hull,
    #              alpha = 0.2,
    #              color = "transparent") +
    stat_ellipse(linetype=2)+
    scale_color_manual(values = cmig) +
    scale_fill_manual(values = cmig)+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme(aspect.ratio=1)
  
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_with_crane_migration_ellipse.pdf", sep = ""),
    width = 5,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
}

p



{

  
  p<-nmds_points %>%
    group_by(year)%>%
    dplyr::mutate(n = dplyr::n())%>%
    #filter(!Family1 == "Gruidae")%>%
    ggplot(aes(x = MDS1, y = MDS2, color = year, fill = year)) + 
    geom_point(size = 0.5)+
    # geom_polygon(data = hull,
    #              alpha = 0.2,
    #              color = "transparent") +
    stat_ellipse(linetype=2)+
    scale_color_npg() +
    scale_fill_npg()+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme(aspect.ratio=1)
  
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_with_crane_sample_year.pdf", sep = ""),
    width = 5,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
}

p


### Just Cranes####
physeq.crane<- subset_samples(physeq, bird == "crane")
psprn.ord.crane <- ordinate(physeq.crane, "NMDS", "bray")
physeq.crane@sam_data$Family1 <- as.factor(physeq.crane@sam_data$Family1)
physeq.crane@sam_data$Trophic.Level <- as.factor(physeq.crane@sam_data$Trophic.Level)
physeq.crane@sam_data$Trophic.Niche <- as.factor(physeq.crane@sam_data$Trophic.Niche)
physeq.crane@sam_data$Migration <- as.factor(physeq.crane@sam_data$Migration)
physeq.crane@sam_data$Primary.Lifestyle <- as.factor(physeq.crane@sam_data$Primary.Lifestyle)

all_meta_ordin <- as_tibble(all_meta)

crane.nmds <- 
  as_tibble(psprn.ord.crane$points,  rownames = "guid") %>%
  left_join(y = all_meta_ordin, by = "guid")

{
  hull <- crane.nmds %>%
    group_by(sci_name) %>%
    slice(chull(MDS1, MDS2))
  
  p <- crane.nmds %>%
    ggplot(aes(
      x = MDS1,
      y = MDS2,
      color = sci_name,
      fill = sci_name
    )) +
    geom_point(size = 2) +
    stat_ellipse(linetype = 2)+
    # geom_polygon(data = hull,
    #              alpha = 0.3,
    #              color = "transparent") +
    scale_color_manual(values = c("grey53", "brown")) +
    scale_fill_manual(values = c("grey53", "brown"))+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme(aspect.ratio=1,
          legend.position = "none")
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_cranes_ellipse.pdf", sep = ""),
    width = 5,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
  }
p


#Lesser by sex

{
  
  p <- nmds_points %>%
    filter(Family1 == "Gruidae",
           sex%in% c("male", "female")) %>%
    ggplot(aes(
      x = MDS1,
      y = MDS2,
      color = sex,
      fill = sex
    )) +
    geom_point(size = 2) +
    stat_ellipse(linetype = 2)+
    # geom_polygon(data = hull,
    #              alpha = 0.3,
    #              color = "transparent") +
    scale_color_manual(values = c("yellow3", "steelblue")) +
    scale_fill_manual(values = c("yellow3", "steelblue"))+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme(aspect.ratio=1,
          legend.position = "none")
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_cranes_all_sex_ellipse.pdf", sep = ""),
    width = 6,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
}
p

#Lesser by sex

{
  
  p <- nmds_points %>%
    filter(Family1 == "Gruidae",
           sci_name ==  "Antigone canadensis canadensis",
           sex%in% c("male", "female")) %>%
    ggplot(aes(
      x = MDS1,
      y = MDS2,
      color = sex,
      fill = sex
    )) +
    geom_point(size = 2) +
    stat_ellipse(linetype = 2)+
    # geom_polygon(data = hull,
    #              alpha = 0.3,
    #              color = "transparent") +
    scale_color_manual(values = c("yellow3", "steelblue")) +
    scale_fill_manual(values = c("yellow3", "steelblue"))+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme(aspect.ratio=1)
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_cranes_lesser_sex.pdf", sep = ""),
    width = 6,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
  }
p

#Greater by sex
{
  
  p <- nmds_points %>%
    filter(Family1 == "Gruidae",
           sci_name ==  "Antigone canadensis tabida",
           sex%in% c("male", "female")) %>%
    ggplot(aes(
      x = MDS1,
      y = MDS2,
      color = sex,
      fill = sex
    )) +
    geom_point(size = 2) +
    stat_ellipse(linetype = 2)+
    # geom_polygon(data = hull,
    #              alpha = 0.3,
    #              color = "transparent") +
    scale_color_manual(values = c("yellow3", "steelblue")) +
    scale_fill_manual(values = c("yellow3", "steelblue"))+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme_update(aspect.ratio=1)
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_cranes_greater_sex.pdf", sep = ""),
    width = 6,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
  }
p


#Lesser by age

{
  
  p <- nmds_points %>%
    filter(Family1 == "Gruidae",
           age%in% c("adult", "juvenile")) %>%
    ggplot(aes(
      x = MDS1,
      y = MDS2,
      color = age,
      fill = age
    )) +
    geom_point(size = 2) +
    stat_ellipse(linetype = 2)+
    # geom_polygon(data = hull,
    #              alpha = 0.3,
    #              color = "transparent") +
    scale_color_manual(values = c("grey", "steelblue3")) +
    scale_fill_manual(values = c("grey", "steelblue3"))+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme(aspect.ratio=1,
          legend.position = "none")
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_cranes_all_age_ellipse.pdf", sep = ""),
    width = 6,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
}
p



#Lesser by age

{
  p <- nmds_points %>%
    filter(Family1 == "Gruidae",
           sci_name ==  "Antigone canadensis canadensis",
           age%in% c("adult", "juvenile")) %>%
    ggplot(aes(
      x = MDS1,
      y = MDS2,
      color = age,
      fill = age
    )) +
    geom_point(size = 2) +
    stat_ellipse(linetype = 2)+
    # geom_polygon(data = hull,
    #              alpha = 0.3,
    #              color = "transparent") +
    scale_color_manual(values = c("grey", "steelblue3")) +
    scale_fill_manual(values = c("grey", "steelblue3"))+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme(aspect.ratio=1)
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_cranes_lesser_age.pdf", sep = ""),
    width = 6,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
}
p

#Greater by age
{
  p <- nmds_points %>%
    filter(Family1 == "Gruidae",
           sci_name ==  "Antigone canadensis tabida",
           age%in% c("adult", "juvenile")) %>%
    ggplot(aes(
      x = MDS1,
      y = MDS2,
      color = age,
      fill = age
    )) +
    geom_point(size = 2) +
    stat_ellipse(linetype = 2)+
   scale_color_manual(values = c("grey", "steelblue3")) +
   scale_fill_manual(values = c("grey", "steelblue3"))+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme_update(aspect.ratio=1)
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_cranes_greater_age.pdf", sep = ""),
    width = 6,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
}
p


# Birds other than cranes
p<- subset_samples(physeq.r, bird == "other")
psprn.ord.non.crane <- ordinate(physeq.non.crane, "NMDS", "bray")
physeq.non.crane@sam_data$Family1 <- as.factor(physeq.non.crane@sam_data$Family1)
physeq.non.crane@sam_data$Trophic.Level <- as.factor(physeq.non.crane@sam_data$Trophic.Level)
physeq.non.crane@sam_data$Trophic.Niche <- as.factor(physeq.non.crane@sam_data$Trophic.Niche)
physeq.non.crane@sam_data$Migration <- as.factor(physeq.non.crane@sam_data$Migration)
physeq.non.crane@sam_data$Primary.Lifestyle <- as.factor(physeq.non.crane@sam_data$Primary.Lifestyle)

all_meta_ordin <- as_tibble(all_meta)

bird.nmds <- 
  as_tibble(psprn.ord.non.crane$points,  rownames = "guid") %>%
  left_join(y = all_meta_ordin, by = "guid")

# all birds order
{
  p <- bird.nmds %>%
    group_by(Order1) %>%
    dplyr::mutate(n = dplyr::n())%>%
    filter(n >=3)%>%
    ggplot(aes(
      x = MDS1,
      y = MDS2,
      color = Order1,
      fill = Order1
    )) +
    geom_point(size = 2) +
    stat_ellipse(linetype = 2)+
    # geom_polygon(data = hull,
    #              alpha = 0.3,
    #              color = "transparent") +
    #scale_color_manual(values = c("grey53", "brown")) +
   # scale_fill_manual(values = c("grey53", "brown"))+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme(aspect.ratio=1)
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_all_birds_order.pdf", sep = ""),
    width = 6,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
  }
p

#Birds by sex

{
  p <- bird.nmds %>%
    filter(
      sex%in% c("male", "female")) %>%
    ggplot(aes(
      x = MDS1,
      y = MDS2,
      color = sex,
      fill = sex
    )) +
    geom_point(size = 2) +
    stat_ellipse(linetype = 2)+
    # geom_polygon(data = hull,
    #              alpha = 0.3,
    #              color = "transparent") +
    scale_color_manual(values = c("yellow3", "steelblue")) +
    scale_fill_manual(values = c("yellow3", "steelblue"))+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme(aspect.ratio=1)
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_birds_non_crane_sex.pdf", sep = ""),
    width = 6,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
}
p

#all birds trophic level
{
  p <- bird.nmds %>%
    ggplot(aes(
      x = MDS1,
      y = MDS2,
      color = Trophic.Level,
      fill = Trophic.Level
    )) +
    geom_point(size = 2) +
    stat_ellipse(linetype = 2)+
    # geom_polygon(data = hull,
    #              alpha = 0.3,
    #              color = "transparent") +
    scale_color_manual(values=c("red4","darkgreen", "gold3"))+
    scale_fill_manual(values=c("red4","darkgreen", "gold3"))+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme(aspect.ratio=1)
  
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_birds_non_crane_trophic_level.pdf", sep = ""),
    width = 6,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
}
p

# all birds migration
{
  p <- bird.nmds %>%
    ggplot(aes(
      x = MDS1,
      y = MDS2,
      color = as.factor(Migration),
      fill = as.factor(Migration)
    )) +
    geom_point(size = 2) +
    stat_ellipse(linetype = 2)+
    # geom_polygon(data = hull,
    #              alpha = 0.3,
    #              color = "transparent") +
    scale_color_manual(values=cmig)+
    scale_fill_manual(values=cmig)+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme(aspect.ratio=1)
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_birds_non_crane_migration.pdf", sep = ""),
    width = 6,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
}
p

# all birds primary lifestyle
{
  p <- bird.nmds %>%
    ggplot(aes(
      x = MDS1,
      y = MDS2,
      color = as.factor(Primary.Lifestyle),
      fill = as.factor(Primary.Lifestyle)
    )) +
    geom_point(size = 2) +
    #stat_ellipse(linetype = 2)+
    geom_polygon(data = hull,
                 alpha = 0.3,
                 color = "transparent") +
    scale_color_manual(values=clife)+
    scale_fill_manual(values=clife)+
    labs(x="NMDS 1", y = "NMDS 2")+
    theme(aspect.ratio=1)
  ggsave(
    p,
    filename = paste(getwd(), "/figures/NMDS_birds_non_crane_lifestyle.pdf", sep = ""),
    width = 6,
    height = 3,
    units = "in",
    device = cairo_pdf
  )
}
p

## Network ####

## Co-occurrence ####

#with core mycobiome at family level
pseq.core.family <- subset_taxa(pseq.core, family != "Other")



core_subset <- subset_taxa(physeq, rownames(tax_table(physeq)) %in% core.taxa)


pseq.core.family <- conglomerate_taxa(core_subset, "family")

coot.f <-co_occurrence(pseq.core.family, treatment = NULL, method ="spearman", rho = 0, p = 1, cores = 1)
coot.f
write.csv(coot.f%>%dplyr::select(-Treatment), file="tables/co-occurrence_core_family.csv")

co_fam_net <-co_occurrence_network(pseq.core.family, treatment = NULL, 
                      classification = 'family', co_occurrence_table = coot.f)

co_fam_net

pseq.anm.sym.rel <- microbiome::transform(anm.sym.physeq.r, "compositional")


#calculate core microbiome fungi that wer present in >50% of samples and greater than 0.01 relative abundance

pseq.anm.sym.core <- core(pseq.anm.sym.rel, detection = 0.5, prevalence = 0.01) #make phyloseq object

pseq.anm.sym.core <- subset_taxa(pseq.anm.sym.core, family != "Other")

pseq.anm.sym.core <- conglomerate_taxa(pseq.anm.sym.core, "family")

coot.as <-co_occurrence(pseq.anm.sym.core, treatment = NULL, method ="spearman", rho = 0, p = 1, cores = 1)
coot.as
write.csv(coot.as%>%dplyr::select(-Treatment), file="tables/co-occurrence_animal_symbiont.csv")

co_anm_sym_net <-co_occurrence_network(pseq.anm.sym.core, treatment = NULL, 
                            classification = 'family', co_occurrence_table = coot.as)

co_anm_sym_net


## Alluvial plot ####

#get vector of cryptococcus-like fungi
cl.v <- read_csv("data/animal_symb_FUNGuild_v2.csv")%>%
  dplyr::select(-c(1) )%>%
  filter(crypt.like == "y")%>%
  dplyr::select(zotu)%>%
  distinct()%>%
  pull(zotu)

alluv<-psmelt(physeq.MS)%>%
  as.data.frame()%>%
  dplyr::select(OTU,  Sample, order, Abundance)%>%
  dplyr::rename(sci_name = Sample)%>%
  left_join( taxonomy%>%rownames_to_column(var="OTU")%>%dplyr::select(OTU, order, family))%>%
  filter(sci_name %in% c("Antigone canadensis canadensis", "Antigone canadensis tabida"),
         Abundance >0)%>%
  group_by(sci_name, order,family, OTU)%>%
  mutate(n = n(),
         order = if_else(family %in%c("Aspergillaceae"), "Aspergillaceae",  order),
         order = if_else(OTU %in% cl.v, "Cryptococcus-like", order))%>%
  ungroup()%>%
  group_by(sci_name, order, OTU)%>%
  summarise(Abundance = sum(Abundance),
            num_otu = sum(n))%>%
  filter(sci_name %in% c("Antigone canadensis canadensis", "Antigone canadensis tabida"),
         order %in% c("Onygenales","Saccharomycetales", "Aspergillaceae", "Cryptococcus-like")|OTU %in% cl.v)%>%
  mutate(order = factor(order, levels = c("Aspergillaceae","Saccharomycetales","Cryptococcus-like", "Onygenales")))


# alluv<-psmelt(physeq.MS)%>%
#   as.data.frame()%>%
#   dplyr::select(OTU,  Sample, order, Abundance)%>%
#   dplyr::rename(sci_name = Sample)%>%
#   left_join( taxonomy%>%rownames_to_column(var="OTU")%>%dplyr::select(OTU, family))%>%
#   filter(sci_name %in% c("Antigone canadensis canadensis", "Antigone canadensis tabida"),
#          Abundance >0)%>%
#   group_by(sci_name,family, OTU)%>%
#   mutate(n = n())%>%
#   ungroup()%>%
#   group_by(sci_name, family, OTU)%>%
#   summarise(Abundance = sum(Abundance),
#             num_otu = sum(n))%>%
#   mutate(family = if_else(family %in%c("Saccharomycetaceae", "Debaryomycetaceae", "Metschnikowiaceae"), "Saccharomycetales",  family))%>%
#   filter(sci_name %in% c("Antigone canadensis canadensis", "Antigone canadensis tabida"),
#          family %in% c("Ajellomycetaceae","Saccharomycetales", "Aspergillaceae" )|OTU %in% cl.v)%>%
#   mutate(family = if_else(OTU %in% cl.v, "Cryptococcus-like", family),
#          family = factor(family, levels = c("Aspergillaceae","Saccharomycetales","Cryptococcus-like", "Ajellomycetaceae")))

alluv%>%
  ungroup()%>%
  dplyr::select(order, OTU)%>%
  group_by(order)%>%
  distinct()%>%
  dplyr::summarise(num_otu = n())

lode_ord<-alluv%>%
  arrange(order)%>%
  ungroup()%>%
  dplyr::select(OTU)%>%
  distinct(.)%>%
  pull(OTU)

# is_alluvia_form(as.data.frame(.), axes = 1:3, silent = TRUE)
#attempt to reorder 
l_ord <- alluv%>%
  mutate(family = as.character(order))%>%
  arrange(OTU)%>%
  arrange(order)%>%
  dplyr::select(OTU)%>%
  pull()

alluv <- alluv%>%
  mutate(OTU = factor(OTU, levels = lode_ord))

alluvial(alluv[,1:3], freq=alluv$num_otu, 
         col = ifelse(alluv$sci_name == "Antigone canadensis canadensis",  "#A51209", "#FFDB58"),
         alpha = 0.8,
         blocks="TRUE",
         ordering = list(
           NULL,
           NULL,
           NULL
         )
)

#Figure out distribution of select pathogenic taxa in crane samples i.e. what percentage of lessers and greaters had apsergillosus
p<-psmelt(physeq)%>%
  as.data.frame()%>%
  mutate(order = if_else(family %in%c("Aspergillaceae"), "Aspergillaceae",  order),
         order = if_else(OTU %in% cl.v, "Cryptococcus-like", order))%>%
  dplyr::select(OTU,  sci_name, order, Abundance)%>%
  left_join( taxonomy%>%rownames_to_column(var="OTU")%>%dplyr::select(OTU, order, family))%>%
  filter(sci_name %in% c("Antigone canadensis canadensis", "Antigone canadensis tabida"),
         Abundance >0)%>%
  group_by(sci_name)%>%
  mutate(sum_zotu_per_crane_subsp = sum(Abundance))%>%
  ungroup()%>%
  group_by(sci_name, order, OTU)%>%
  dplyr::mutate(otu_reads_per_subspeccies = sum(Abundance),
                perc_prev = otu_reads_per_subspeccies/sum_zotu_per_crane_subsp)%>%
  filter(order %in% c("Onygenales","Saccharomycetales", "Aspergillaceae", "Cryptococcus-like")|OTU %in% cl.v)%>%
  mutate(order = factor(order, levels = c("Aspergillaceae","Saccharomycetales","Cryptococcus-like", "Onygenales")))%>%
  ungroup()%>%
  dplyr::select(OTU, sci_name, perc_prev, order)%>%
  distinct()%>%
  mutate(OTU = factor(OTU, levels = lode_ord))%>%
  ggplot(., aes(x = perc_prev, y=OTU,  fill = sci_name))+
  geom_col(position = position_dodge())+
  scale_fill_manual(values = c("#A51209", "#FFDB58"))+
  scale_x_log10()+
  theme(legend.position = "none")
  
p
ggsave(p, filename="figures/Prevelance_otus.pdf", width=3, height=7, device = cairo_pdf)

#Animal symbiont composition plot
anm.sym <- read_csv("data/animal_symb_FUNGuild_v2.csv")%>%
  dplyr::select(-c(1))

p<-psmelt(physeq)%>%
  as.data.frame()%>%
  dplyr::select(OTU,  sci_name, Sample, family, Abundance)%>%
  #dplyr::rename(sci_name = Sample)%>%
  left_join(anm.sym%>%dplyr::rename(OTU=zotu))%>%
  filter(sci_name %in% c("Antigone canadensis canadensis", "Antigone canadensis tabida"),
         Abundance >0)

#lessers
p1<-p %>%
  filter(sci_name == "Antigone canadensis canadensis") %>%
  group_by(Sample) %>%
  mutate(total_zotu = n()) %>%
  ungroup() %>%
  group_by(Sample, animal.trophic) %>%
  mutate(
    anm_cat_count = n(),
    perc_anm_symb = anm_cat_count / total_zotu,
    
    animal.trophic = factor(animal.trophic, levels = c("y", "n", "u"))) %>%
  ungroup() %>%
  dplyr::select(Sample, perc_anm_symb, animal.trophic) %>%
  distinct() %>%
  # Arrange by descending perc_anm_symb within each animal.trophic
  arrange(animal.trophic, desc(-perc_anm_symb)) %>%
  # Set Sample as a factor with ordered levels
  mutate(Sample = factor(Sample, levels = unique(Sample))) %>%
  ggplot(aes(Sample, perc_anm_symb, fill = animal.trophic)) +
  geom_col(width=1) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c("darkturquoise", "dodgerblue", "grey53")) +
  labs(
    x = "",
    y = "Proportion",
    fill = "Animal Trophic"
  )

p2<-p %>%
  filter(sci_name == "Antigone canadensis tabida") %>%
  group_by(Sample) %>%
  mutate(total_zotu = n()) %>%
  ungroup() %>%
  group_by(Sample, animal.trophic) %>%
  mutate(
    anm_cat_count = n(),
    perc_anm_symb = anm_cat_count / total_zotu,
  
    animal.trophic = factor(animal.trophic, levels = c("y", "n", "u"))) %>%
  ungroup() %>%
  dplyr::select(Sample, perc_anm_symb, animal.trophic) %>%
  distinct() %>%
  # Arrange by descending perc_anm_symb within each animal.trophic
  arrange(animal.trophic, desc(-perc_anm_symb)) %>%
  # Set Sample as a factor with ordered levels
  mutate(Sample = factor(Sample, levels = unique(Sample))) %>%
  ggplot(aes(Sample, perc_anm_symb, fill = animal.trophic)) +
  geom_col(width=0.7) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c("darkturquoise", "dodgerblue", "grey53")) +
  labs(
    x = "",
    y = "Proportion",
    fill = "Animal Trophic"
  )

g <- gridExtra::grid.arrange(p1, p2, ncol=1)

ggsave(g, filename= "figures/animal_symbiont_cranes_2_panel.pdf", width=8, height=5, device = cairo_pdf)

# Plot
ggplot(data, aes(Sample, perc_anm_symb, fill = animal.trophic)) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "Proportion of Animal.Trophic Categories per Sample",
    x = "Sample",
    y = "Proportion",
    fill = "Animal Trophic"
  )


p<-psmelt(physeq)%>%
  as.data.frame()%>%
  mutate(order = if_else(family %in%c("Aspergillaceae"), "Aspergillaceae",  order),
         order = if_else(OTU %in% cl.v, "Cryptococcus-like", order))%>%
  dplyr::select(OTU,  sci_name, order, Abundance)%>%
  left_join( taxonomy%>%rownames_to_column(var="OTU")%>%dplyr::select(OTU, order, family))%>%
  filter(sci_name %in% c("Antigone canadensis canadensis", "Antigone canadensis tabida"),
         Abundance >0)%>%
  group_by(sci_name, order)%>%
  dplyr::summarise(sum_abun = sum(Abundance))  
    
    
## Percent prevelance pathogen orders####
psmelt(physeq)%>%
  as.data.frame()%>%
  mutate(order = if_else(family %in%c("Aspergillaceae"), "Aspergillaceae",  order),
         order = if_else(OTU %in% cl.v, "Cryptococcus-like", order))%>%
  dplyr::select(OTU,  sci_name, order, Abundance)%>%
  left_join( taxonomy%>%rownames_to_column(var="OTU")%>%dplyr::select(OTU, order, family))%>%
  filter(sci_name %in% c("Antigone canadensis canadensis", "Antigone canadensis tabida"),
         Abundance >0)%>%
  group_by(sci_name)%>%
  mutate(sum_reads_per_crane_subsp = sum(Abundance))%>%
  ungroup()%>%
  group_by(sci_name, order)%>%
  dplyr::mutate(order_reads_per_subspeccies = sum(Abundance),
                perc_prev = order_reads_per_subspeccies/sum_reads_per_crane_subsp)%>%
  filter(order %in% c("Onygenales","Saccharomycetales", "Aspergillaceae", "Cryptococcus-like")|OTU %in% cl.v)%>%
  mutate(order = factor(order, levels = c("Aspergillaceae","Saccharomycetales","Cryptococcus-like", "Onygenales")))%>%
  ungroup()%>%
  dplyr::select(sci_name, perc_prev, order)%>%
  distinct())
