# GDM ####

library(gdm)
library(tidyverse)
library(vegan)
library(phyloseq)
require(parallel)
require(doParallel)
require(tidyverse)
require(foreach)
require(ape)
library(dutchmasters)


source("~/Desktop/ggplot_themes/ggplot_themes.R")
theme_set(theme_arial_clean())

#color palette
n <- 10

h1 <- hcl.colors(n, palette = "Dynamic")
h2 <- hcl.colors(n, palette = "Earth")
h3 <- hcl.colors(n, palette = "Berlin")
h4 <- hcl.colors(n, palette = "Fall")
h5 <- hcl.colors(n, palette = "Sunset")

library(unikn)
seecol(list(h1, h2, h3, h4, h5), 
       col_brd = "white", lwd_brd = 4, 
       title = "Example palettes from hcl.colors(n = 10)", 
       pal_names = c("Dynamic", "Earth", "Berlin", "Fall",  "Sunset"))

contr.col <- c(h2[6], h2[2], h1[3], h4[10], h2[9])

extrafont::loadfonts()

#make max clade cred from Linck et al. 2022
tree <- read.tree(file =paste("data/bird_myco_tree_processed.tre", sep=""))

#load merged by species phyloseq object
load("data/physeq_merged_spp.Rdata")

trait_tab <-read_csv("data/Bird_predictors.csv")%>%
  #filter(guid %in% rand.guid$guid)%>%
  #mutate(Mass = if_else(is.na(weight), Mass, weight))%>%
  dplyr::select(sci_name, Range.Size, Mass, migratory_distance, Hand.Wing.Index, lat, long)%>%
  na.omit()%>%
  group_by(sci_name)%>%
  dplyr::mutate(Mass = mean(Mass, na.rm=T),
         sci_name = gsub(" ", "_", sci_name))%>%
  slice_sample(n=1)

# Do it again with rarefied compositional data like Jenn Rudgers did, but note that we summed otus during merge and Jen took mean
trait_tab$migratory_distance <-  trait_tab$migratory_distance /max(trait_tab$migratory_distance)
trait_tab$Hand.Wing.Index <- trait_tab$Hand.Wing.Index /max(trait_tab$Hand.Wing.Index)
trait_tab$Mass <- trait_tab$Mass /max(trait_tab$Mass)
trait_tab$Range.Size <- trait_tab$Range.Size /max(trait_tab$Range.Size)

#get phylogenetic distances
phy.d <-cophenetic.phylo(tree) 
phy.d<-phy.d / max(phy.d)
phy.d<-phy.d%>%
  as.data.frame()%>%
  rownames_to_column(var="sci_name")


# Repeatedly rareify and generate 1000 gdms
iters <- 1000
dev.exp <- vector(mode='list', length=iters)
splines <- vector(mode='list', length=iters)
varImp <- vector(mode='list', length=iters)

gdm.list <- list(dev.exp, splines, varImp)

#set seed
set.seed(3892594)

job::job({
  for (i in 1:iters) {
    # This method repeatedly rarefies
    
    physeq.m.r.spp <- rarefy_even_depth(physeq.MS,
                                        sample.size = min(sample_sums(physeq.MS)),
                                        replace = F)
    
    #long format species table
    sppTab.l <- psmelt(physeq.m.r.spp) %>%
      dplyr::select(Sample, Abundance, OTU) %>%
      dplyr::rename(sci_name = Sample) %>%
      mutate(sci_name = gsub(" ", "_", sci_name)) %>%
      left_join(trait_tab %>% dplyr::select(sci_name, lat, long))
    
    # # format species data
    # d.data<-physeq.m.r.spp@otu_table%>%
    #   as.data.frame()%>%
    #   rownames_to_column(var="guid")%>%
    #   mutate(guid = gsub(" ", "_", guid))%>%
    #   left_join(samtab%>%as.tibble()%>%dplyr::select(guid, sci_name))%>%
    #   dplyr::select(-sci_name)%>%
    #   column_to_rownames(var="guid")%>%
    #   avgdist(., dmethod = "bray", sample=400)  %>%
    #   as.matrix()
    #
    # spptab <- d.data%>%
    #   as.data.frame()%>%
    #   rownames_to_column(var="sci_name")
    
    trait_tab <- as.data.frame(trait_tab)
    
    comp.gdm.temp <- formatsitepair(
      bioData = sppTab.l,
      abundance = T,
      bioFormat = 2,
      #species site pair
      XColumn = "long",
      YColumn = "lat",
      abundColumn = "Abundance",
      predData = trait_tab,
      sppColumn = "OTU",
      distPreds = list(as.matrix(phy.d)),
      siteColumn = "sci_name"
    )
    
    gdm.list[[1]][[i]] <- gdm(comp.gdm.temp, geo = T)
    
    # extract splines
    gdm.list[[2]][[i]] <- data.frame(isplineExtract(gdm.list[[1]][[i]])) %>%
      dplyr::select(!contains("y.")) %>%
      pivot_longer(
        cols = c(1:length(.)),
        names_to = c("predictor"),
        values_to = "x"
      ) %>%
      bind_cols(
        .,
        data.frame(isplineExtract(gdm.list[[1]][[i]])) %>%
          dplyr::select(contains("y.")) %>%
          pivot_longer(
            cols = c(1:length(.)),
            names_to = c("predictor"),
            values_to = "y"
          ) %>%
          mutate(predictor = gsub("y.", "x.", predictor))
      ) %>%
      mutate(iteration = paste(i)) %>%
      dplyr::rename(predictor = predictor...1)
    
    gdm.list[[3]][[i]] <- gdm.varImp(
      comp.gdm.temp,
      geo = F,
      nPerm = 500,
      parallel = F
    )
    
    print(paste0("Finished with rarefied gdm:", i))
  }
  
  save(gdm.list, file = "data/gdm_Geo.RData")
  #load(file = "data/gdm_noGeo.RData")
  
  t <- do.call("rbind", gdm.list[[2]]) #bind splines for plotting
  
  #add full model predictions
  
  p <- ggplot(t %>% filter(iteration > 0), aes(x, y, group = iteration)) +
    geom_line(
      data = t %>% dplyr::filter(predictor == "x.Hand.Wing.Index"),
      aes(x, y),
      alpha = 0.02,
      color = "#D8DFC3"
    ) +
    geom_line(
      data = t %>% dplyr::filter(predictor == "x.migratory_distance"),
      aes(x, y),
      alpha = 0.02,
      color = "#C7522B"
    ) +
    geom_line(
      data = t %>% dplyr::filter(predictor == "x.matrix_1"),
      aes(x, y),
      alpha = 0.02,
      color = "#9DB469"
    ) +
    geom_line(
      data = t %>% dplyr::filter(predictor == "x.Range.Size"),
      aes(x, y),
      alpha = 0.02,
      color = "#509F9E"
    ) +
    geom_line(
      data = t %>% dplyr::filter(predictor == "x.Mass"),
      aes(x, y),
      alpha = 0.02,
      color = "#B2874C"
    ) +
    #geom_line(data=t%>%filter(iteration==0),aes( x, y), linewidth=1, linetype="dashed")+
    #facet_wrap(.~predictor, scales="free")+
    theme(aspect.ratio = 1) +
    labs(x = "Predictor dissimilarity", y = "Partial ecological distance")
  
  p
  
  ggsave(
    p,
    file = paste0("figures/gdm_dissimilarity_geo", iters, ".pdf"),
    width = 8,
    height = 5
  )
  
  
  t <- vector(length = iters) # create vector for deviance explained
  
  for (i in 1:iters) {
    t[i] <- gdm.list[[1]][[i]][["explained"]]
    
  }
})

#plot deviance explained
p<-as.data.frame(t)%>%
  ggplot(aes(t))+
  geom_density(linewidth=1)+
  labs(x= "Mean deviance explained")

p 
ggsave(p, file=paste0("figures/gdm_mean_deviance_explained_n", iters, ".pdf"), width=4, height=4)
# explained deviance low generally

as.data.frame(t)%>%
  dplyr::summarise(mean = mean(t),
                   sd = sd(t))

# mean        sd
# 1 9.023566 0.2514831

#Model coefficients
modcoef <-list() 

for(i in 1:iters){
  
  preds<-gdm.list[[1]][[i]][["predictors"]]
  modcoef[[i]]<-as.data.frame(gdm.list[[1]][[i]][["coefficients"]])%>% mutate(predictors = rep(preds, each=3))%>%dplyr::rename(coef=1)
  modcoef[[i]]$iter <- i
  
}

modcoef <- do.call("rbind", modcoef)

modcoef%>%
  group_by(predictors)%>%
  dplyr::summarise(mean = mean(coef, na.rm=T),
                   sd = sd(coef, na.rm=T))

# find vraible importance for resampled models
varimp<-list() # create vector for deviance explained

for(i in 1:iters){
  
  varimp[[i]]<-as.data.frame(gdm.list[[3]][[i]][["Predictor Importance"]])%>%rownames_to_column(var="predictor")
  varimp[[i]]$iter <- i
  
}

varimp <- do.call("rbind", varimp)

varimp%>%
  group_by(predictor)%>%
  dplyr::summarise(mean = mean(`All predictors`),
                   sd = sd(`All predictors`))


p<-varimp%>%
  group_by(predictor)%>%
  dplyr::summarise(mean = mean(`All predictors`),
                   sd = sd(`All predictors`))%>%
  mutate(predictor = fct_reorder(predictor, -mean))%>%
  ggplot(., aes(predictor, mean, fill = predictor))+
  geom_col(colour ="black")+
  geom_errorbar(aes(predictor, mean, ymin = mean-sd,
                                ymax = mean+sd), width=0.3)+
  scale_fill_manual(values = contr.col)+
  labs(y = "Mean variable importance", x=NULL)+
  theme(legend.position = "none", 
        aspect.ratio = 1)
p

ggsave(p, file="figures/variable_importance_boxplot.pdf", width = 5, height = 5)


# Get p values for all model
modass<-list() # create vector for deviance explained

for(i in 1:iters){
  
  modass[[i]]<-as.data.frame(gdm.list[[3]][[i]][["Model assessment"]])%>%rownames_to_column(var="predictor")
  modass[[i]]$iter <- i
  
}

modass <- do.call("rbind", modass)

p<-ggplot(modass%>%filter(predictor == "Model p-value"), aes(x=`All predictors`))+
  geom_density(adjust=3)+
  geom_vline(xintercept = 0.05, linetype="dashed")
p

ggsave(p, file="figures/model_p_value_density.pdf", width = 5, height = 5)

#Count variables in each model
varimp%>%
  group_by(predictor)%>%
  dplyr::summarise(n = n())



# Get p values for all predictors
predp<-list() # create vector for deviance explained

for(i in 1:iters){
  
  predp[[i]]<-as.data.frame(gdm.list[[3]][[i]][["Predictor p-values"]])%>%rownames_to_column(var="predictor")
  predp[[i]]$iter <- i
  
}

predp <- do.call("rbind", predp)

p<-ggplot(predp, aes(x=`All predictors`, y = predictor, fill=predictor))+
  geom_boxplot()+
  geom_vline(xintercept = 0.05, linetype="dashed")+
  theme(legend.position = "none")
p

ggsave(p, file="figures/predictor_p_value_boxplot.pdf", width = 5, height = 5)



