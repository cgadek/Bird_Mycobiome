# Phyloseq Wrangler#### 
#combine and process zOTU, taxonomy, host traits, and host phylogeny files for analysis

# Load Packages, Themes, etc... ####
pacman::p_load(phyloseq,
               tidyverse,
               ape,
               plyr,
               phylosmith,
               ggClusterNet,
               microViz,
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
               reshape2)

source("functions/Merging_functions.R")

# Read in data /minor cleaning####
## Transposed otu table
otu <- read.table("data/zotutab_v4.txt", header = TRUE, sep = '\t', row.names=1)
motus<-as.matrix(otu)
class(motus) <- "numeric"

#load in FunGuild table to make animal symbiont subset

anm.sym <- read_csv("data/animal_symb_FUNGuild_v2.csv")%>%
  dplyr::select(-c(1))

taxonomy <- read.table("data/zotus_tax_v4.txt", header = TRUE, sep = '\t', row.names = 1)

#get rid of taxa names and zotu with spaces
taxonomy <- taxonomy%>%
  mutate( genus = gsub(" ", "", genus))%>%
  rownames_to_column("zotu")%>%
  mutate(zotu = gsub(" ", "", zotu))%>%
  column_to_rownames("zotu")

mtax<-as.matrix(taxonomy)

otu.otu <- as.data.frame(colnames(otu))
otu.tax <- as.data.frame(rownames(taxonomy))

not.shared <- otu.otu%>%
  filter(!`colnames(otu)` %in% otu.tax$`rownames(taxonomy)`)
#NOTE there is single zotu that does not match the otu file. Paris identified as Basidiomycete that did not have hit in updated UNITE so we are proceeding without zOTU 157

#all_meta <- read.table("data/all_lungs_v3.txt", header = TRUE, sep = '\t') #version three has all the avonet and arctos data integrated

all_meta <- read.csv("data/Bird_predictors.csv")

#filter out non-amplicon birds
all_meta <- all_meta%>%
  filter(!is.na(guid))%>%
  mutate(guid.row = guid)%>%
  column_to_rownames("guid.row")

# Merge into a phyloseq object ####
otutab<-otu_table(motus,taxa_are_rows=FALSE)
sample_names(otutab) <- paste0( "S", sample_names(otutab)) #need to have a character in front
taxtab<-tax_table(mtax)
samtab<-sample_data(all_meta)
physeq<-merge_phyloseq(otutab,taxtab,samtab)

#Rarefy dataset####
#this is standard practice for microbiome research
#NOTE according to function documentation this must be used with untrimmed data; but Michael and internet recommend rarefication
set.seed(1982)
physeq.r <- rarefy_even_depth(physeq, rngseed=1, sample.size= min(sample_sums(physeq)), replace=F)

# we chose rarified based on minimum sample size of data
# 180TUs were removed because they are no longer 
#present in any sample after random subsampling

#subset by animal symbionts
anm.sym.physeq <-tax_select(physeq, tax_list = anm.sym%>%filter(animal.trophic =="y")%>%dplyr::select(cat)%>%unique()%>%pull(cat))

#subset by non-animal symbionts
not.anm.sym.physeq <-tax_select(physeq, tax_list = anm.sym%>%filter(animal.trophic =="n")%>%
                                  dplyr::select(cat)%>%
                                  unique()%>%
                                  filter(!cat %in% c("Hypocreales", "Eurotiales", "Mucorales", "Saccharomycetales", "-"))%>%
                                  pull(cat))

#rarefy these
anm.sym.physeq.r <- rarefy_even_depth(anm.sym.physeq, rngseed=1, sample.size= min(sample_sums(anm.sym.physeq)), replace=F)

not.anm.sym.physeq.r <- rarefy_even_depth(not.anm.sym.physeq, rngseed=1, sample.size= min(sample_sums(not.anm.sym.physeq)), replace=F)

# merge by species
mergeby <- "sci_name"

#Merge by samples this sums all counts at samples
physeq.MS <- merge_samples(physeq, mergeby)
physeq.MS.r <- merge_samples(physeq.r, mergeby)

anm.sym.physeq.MS <- merge_samples(anm.sym.physeq, mergeby)
anm.sym.physeq.MS.r <- merge_samples(anm.sym.physeq.r, mergeby)

not.anm.sym.physeq.MS <- merge_samples(not.anm.sym.physeq, mergeby)
not.anm.sym.physeq.MS.r <- merge_samples(not.anm.sym.physeq.r, mergeby)


# To simplify the GP data
#ps = subset_taxa(GlobalPatterns, Kingdom == "Bacteria")
#ps = prune_taxa(taxa_sums(ps) > 1000, ps)

total = median(sample_sums(physeq))
standf = function(x, t=total) round(t * (x / sum(x)))
normalized_ps = transform_sample_counts(physeq, standf)

#save as .Rdata files
save(physeq, file ="data/physeq_all.Rdata") #not rarefied all individuals except for suspicious cranes
save(physeq.r, file ="data/physeq_rarefied.Rdata") #rarefied physeq object
save(physeq.MS, file ="data/physeq_merged_spp.Rdata") #physeq merged by species not rarefied
save(physeq.MS.r, file ="data/physeq_rarefied_merged_spp.Rdata") #physeq merged by species and rarefied
save(anm.sym.physeq, file ="data/physeq_animal_symbiont.Rdata") #physeq animal symbiont not rarefied
save(not.anm.sym.physeq, file ="data/physeq_non_animal_symbiont.Rdata") #physeq non-animal symbiont not rarefied
save(anm.sym.physeq.r, file ="data/physeq_animal_symbiont.Rdata") #physeq animal symbiont rarefied
save(not.anm.sym.physeq.r, file ="data/physeq_non_animal_symbiont.Rdata") #physeq non-animal symbiont rarefied


# Create data frame for GDM####
#rarefied community data, merged(summed) by species (includes subsepcies of cranes), now obtaining relative/compositional abundance per sample
physeq.MS.r%>%
  transform_sample_counts(., function(OTU) OTU/sum(OTU))%>%
  psmelt()%>%
  dplyr::select(OTU, Sample, Abundance)%>%
  write_csv(file="data/gdm_compositional.csv")

# host phylo tree ####
#make max clade cred from Linck et al. 2022
tree <- ape::read.tree(file =paste("~/Dropbox/Research/Marsh Birds/Marsh_birds/data/birds_mcc.tre", sep=""))

tips<-as.data.frame(tree$tip.label)

#Some tip labels have old names and we need to change to match our data and AVONET
tree$tip.label<-gsub("Anas_strepera", "Mareca_strepera", tree$tip.label)
tree$tip.label<-gsub("Spinus_pinus", "Carduelis_pinus", tree$tip.label)
tree$tip.label<-gsub("Carduelis_psaltria", "Spinus_psaltria",  tree$tip.label)
tree$tip.label<-gsub("Haemorhous_mexicanus", "Carpodacus_mexicanus",  tree$tip.label)
tree$tip.label<-gsub("Vermivora_celata", "Oreothlypis_celata", tree$tip.label)
tree$tip.label<-gsub("Geococcyx_californianus", "Geococcyx_californicus", tree$tip.label)
tree$tip.label<-gsub("Grus_canadensis", "Antigone_canadensis tabida", tree$tip.label)

tree <- keep.tip(tree, all_meta%>%
                   filter(sci_name != "Antigone canadensis canadensis")%>%
                   dplyr::select(sci_name)%>%
                   distinct()%>%
                   mutate(sci_name = gsub(" ", "_", sci_name),
                          sci_name = gsub("densis_tabida", "densis tabida", sci_name))%>%
                   pull(sci_name))

plot(tree, cex=0.5)
nodelabels(cex=0.2)


tree<-phytools::bind.tip(tree, "Antigone_canadensis canadensis", edge.length=NULL, where=31, position=1)

d<-cophenetic.phylo(tree)

plot(tree, cex=0.2)
nodelabels()

plot.phylo(tree, type="tidy",edge.width = 2, cex=0.5)
write.tree(tree, "data/bird_myco_tree_processed.tre")

