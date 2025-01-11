# Jacccard/Bray-Curtis heatmap ####

pacman::p_load(phyloseq,
               tidyverse,
               plyr,
               ape,
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
               microViz, 
               vegan)

#Read in physeq object processed by Phyloseq_wrangler.R file####
load("data/physeq_all.Rdata")
load("data/physeq_rarefied.Rdata")
load("data/physeq_merged_spp.Rdata")
load("data/physeq_rarefied_merged_spp.Rdata")

tree <- read.tree("data/bird_myco_spp_level.tre")

tips <- gsub("_", " ", tree$tip.label)

#read in processed sample data  and taxonomy file
all_meta <- read_csv("data/all_meta.csv")
taxonomy <- read.table("data/zotus_tax_v3.txt", header = TRUE, sep = '\t', row.names = 1)
samtab<-sample_data(all_meta)

# Lesser
{d.data<-physeq.r@otu_table%>%
  as.data.frame()%>%
  rownames_to_column(var="guid")%>%
  left_join(samtab%>%as.tibble()%>%dplyr::select(guid, bird, sci_name))%>%
  filter(bird=="crane" & sci_name =="Antigone canadensis canadensis")%>%
  dplyr::select(-bird, -sci_name)%>%
  column_to_rownames(var="guid")

bray <-avgdist(d.data, dmethod = "bray", sample=400)  %>%
  as.matrix()%>%
  as.tibble(rownames = "A")%>%
  pivot_longer(-A, names_to = "B", values_to = "distances")
  
bray%>%
  ggplot(aes(x=A, y = B, fill=distances))+
  geom_tile()}

# Greater
{d.data<-physeq.r@otu_table%>%
    as.data.frame()%>%
    rownames_to_column(var="guid")%>%
    left_join(samtab%>%as.tibble()%>%dplyr::select(guid, bird, sci_name))%>%
    filter(bird=="crane" & sci_name =="Antigone canadensis tabida")%>%
    dplyr::select(-bird, -sci_name)%>%
    column_to_rownames(var="guid")
  
  bray <-avgdist(d.data, dmethod = "bray", sample=400)  %>%
    as.matrix()%>%
    as.tibble(rownames = "A")%>%
    pivot_longer(-A, names_to = "B", values_to = "distances")
  
  bray%>%
    ggplot(aes(x=A, y = B, fill=distances))+
    geom_tile()
}


{d.data<-physeq.r@otu_table%>%
    as.data.frame()%>%
    rownames_to_column(var="guid")%>%
    left_join(samtab%>%as.tibble()%>%dplyr::select(guid, bird, sci_name))%>%
    #filter(bird=="crane" & sci_name =="Antigone canadensis tabida")%>%
    dplyr::select(-bird, -sci_name)%>%
    column_to_rownames(var="guid")
  
  bray <-avgdist(d.data, dmethod = "bray", sample=400)  %>%
    as.matrix()%>%
    as.tibble(rownames = "A")%>%
    pivot_longer(-A, names_to = "B", values_to = "distances")
  
  bray%>%
    ggplot(aes(x=A, y = B, fill=distances))+
    geom_tile()
}


{d.data<-physeq.MS.r@otu_table%>%
    as.data.frame()%>%
    rownames_to_column(var="guid")%>%
    left_join(samtab%>%as.tibble()%>%dplyr::select(guid, bird, sci_name))%>%
    #filter(bird=="crane" & sci_name =="Antigone canadensis tabida")%>%
    dplyr::select(-bird, -sci_name)%>%
    column_to_rownames(var="guid")
  
  bray <-avgdist(d.data, dmethod = "bray", sample=400)  %>%
    as.matrix()%>%
    as.tibble(rownames = "A")%>%
    pivot_longer(-A, names_to = "B", values_to = "distances")
  
  bray%>%
    ggplot(aes(x=A, y = B, fill=distances))+
    geom_tile()
}


# All species merged

#make epmty vector for distance matrices
d.v <- vector(mode="list", 1000)
set.seed(1982)

{
  for(i in 1:length(d.v)){
  #repeatedly rarefy 
  physeq.m.r.spp <- rarefy_even_depth(physeq.MS, replace=T)
  
  d.v[[i]]<-physeq.m.r.spp@otu_table%>%
    as.data.frame()%>%
    rownames_to_column(var="guid")%>%
    left_join(samtab%>%as.tibble()%>%dplyr::select(guid, bird, sci_name))%>%
    dplyr::select(-bird, -sci_name)%>%
    column_to_rownames(var="guid")%>%
    avgdist(., dmethod = "bray", sample=500)%>%
    as.matrix()
  
  }
  
  save(d.v, file="data/repeatedly_rarefied_beta_matrices.Rdata")
  
}

load(file="data/repeatedly_rarefied_beta_matrices.Rdata")
pt<-Reduce(`+`, d.v) / length(d.v)
#heatmap( Reduce(`+`, d.v) / length(d.v),scale="column")

pt <- pt[tips,tips]
#plot heat map in order of bird phylogeny
gplots::heatmap.2(pt, Rowv = as.dendrogram(tree), Colv = as.dendrogram(tree), trace="none", margins = c(9, 9), col = hcl.colors(25, palette = "cividis") )

#get dendrogram in beta diversity order
gplots::heatmap.2(pt, trace="none", margins = c(9, 9), col = hcl.colors(25, palette = "cividis") )

#Now heatmap of phylo distances
phy.d <- cophenetic.phylo(tree)
gplots::heatmap.2(phy.d, Rowv = as.dendrogram(tree), Colv = as.dendrogram(tree), trace="none", margins = c(9, 9), col = hcl.colors(25, palette = "mako") )

physeq.MS <-merge_samples(physeq, "sci_name")

ps.o<-tax_glom(physeq.MS, "order", NArm = TRUE)

sample_data(ps.o)$sci_name <- as.character(sample_data(ps.o)$sci_name)

tip <- as.data.frame(tree$tip.label)%>%
  mutate(tip = gsub("_", " ",.[[1]]))%>%
  pull(tip)

sample_order <- sample_names(physeq.MS)
sample_order <- sample_order[order(match(sample_order,tip))]


ps.o <- ps.o%>%
  ps_reorder(sample_order)

om <-ps.o@otu_table%>%
  as.data.frame()%>%
  rownames_to_column(var="sci_name")%>%
  pivot_longer(cols =  contains("otu"), names_to = c("OTU"), values_to = "Abundance") %>% #make it long form
  left_join(taxonomy%>%rownames_to_column(var="OTU"), by = "OTU") %>%
  mutate(order = replace_na(order, "Unknown"),
         order = if_else(order=="", "Unknown", family),
         order = if_else(order ==" ", "Unknown", family)) %>%
  group_by(sci_name, order) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup() %>%
  filter(!order =="")%>%
  pivot_wider(names_from = c(order), values_from = Abundance, values_fill = list(Abundance = 0))


om<-om[order(factor(om$sci_name, levels=unique(sample_order))),]

om<-om%>%
  column_to_rownames(var="sci_name")%>%
  as.matrix()

heatmap(as.matrix(om), scale="column")

om <-  log10(om+0.00001)

bray_dist <- function(x) {
  vegdist(x, method = "bray", na.rm = T)
}

ComplexHeatmap::Heatmap(om, row_names_side = "left", column_names_side = "bottom", 
        row_dend_side = "left", rect_gp = grid::gpar(col = "grey"), 
        row_names_gp = grid::gpar(cex=0.5, fontface = "bold"), 
        column_names_gp = grid::gpar(cex=0.5, fontface = "bold"), 
        row_dend_width = unit(4, "cm"), column_dend_height = unit(3, "cm"), 
        column_title = NULL, column_names_rot = 35, row_title = NULL,
        cluster_rows = F, clustering_distance_columns = bray_dist, 
        clustering_method_columns = "ward.D2", 
        row_km = 3, column_km = 4)

gplots::heatmap.2(om, trace="none", margins = c(9, 9), col = hcl.colors(25, palette = "cividis") )

p<- plot_heatmap(ps.o, "NMDS", "bray", "sci_name", "order")
p

 ps.o@sam_data$sci_name
 
physeq.m.spp
plot_heatmap(subset_taxa(physeq.m.r.spp, order == "Eurotiales"), method = "NMDS", distance="bray")

gplots::heatmap.2(as.matrix(om),  main="Host UniFracs", trace="none")
