### create tree for mycobiome results
library(ape)
library(tidyverse)
library(phytools)
library(phylobase)
library(adephylo)
library(phylosignal)
library(picante)
library(ade4)
library(AER)
library(brms)
library(geiger)
library(phytools)
library(nlme)
library(visreg)
library(MASS)

set.seed(1982)

#import extra fonts or else there will be issues with custom theme
extrafont::font_import() #must select yes
y
extrafont::loadfonts()

# read in bird data from don
#df <- read_csv("data/natvig_bird_lung_loan_cleaned_county_centroids.csv")
#df <- read_csv("data/phylo_df_merged.csv") #this has community metrics merged by species
df <-read_csv("data/phylo_df_unmerged.csv") #this one takes median of individual bird community metrics and manually adjusted for crane subspecies i.e estimated tabida latitude, lon, mass, etc..
#df <-read_csv("data/phylo_df_unmerged_rarefied.csv")

phylo.lm.df <- df%>%
  mutate(sci_name = str_replace(sci_name, " ", "_"))

#read in phylo tree
tree <- read.tree("data/bird_myco_tree_processed.tre")


#Migration in AVONET data has three categories 
# "1 = Sedentary. 
# 2 = Partially migratory, i.e. minority of population migrates long distances, or most of population undergoes short-distance migration, nomadic movements, distinct altitudinal migration, etc.
# 3 = Migratory, i.e. majority of population undertakes long-distance migration"


# Phylogenetic signal
# create phylo 4d object
#start with Shannon and sci_name as rownames
###All otus
er <-estimate_richness(physeq.r, measures=c("Observed", "InvSimpson", "Shannon", "Chao1", "ACE"))%>%
  rownames_to_column(var = "guid")%>%
  left_join(all_meta)%>%
  group_by(sci_name)%>%
  dplyr::summarise(Observed = mean(Observed, na.rm=T),
                   Shannon = mean(Shannon, na.rm=T),
                   Chao1 = mean(Chao1, na.rm=T),
                   InvSimpson = mean(InvSimpson, na.rm=T))


dat <- er%>%
  dplyr::select(sci_name, Observed, Shannon,InvSimpson, Chao1)%>%
  mutate(sci_name = gsub(' ', "_", sci_name))%>%
  column_to_rownames(., var="sci_name")

p4d <- phylo4d(tree, dat[tree$tip.label,])
barplot(p4d, tree.ladderize =T, tree.ratio = 0.33, grid.vertical = F, show.box = T)


#Test phylosignal
## Moran's I
moran.test <- abouheif.moran(p4d,method="Abouheif", nrepet = 10000, alter = "two-sided")

moran.test
plot(moran.test)
# Adjustment method for multiple comparisons:   none 
# Permutation number:   10000 
#       Test         Obs    Std.Obs     Alter    Pvalue
# 1   Observed -0.04686905 -0.1512166 two-sided 0.8856114
# 2    Shannon -0.14560085 -0.9659180 two-sided 0.3448655
# 3 InvSimpson -0.06740234 -0.4482736 two-sided 0.6084392
# 4      Chao1  0.01800978  0.3833496 two-sided 0.7148285

##Abouheifs Cmean
abhouief.test <- abouheif.moran(p4d,method="oriAbouheif", nrepet = 10000, alter = "two-sided")

abhouief.test
# Adjustment method for multiple comparisons:   none 
# Permutation number:   10000 
#       Test         Obs    Std.Obs     Alter    Pvalue
# 1   Observed -0.03593955 -0.3071405 two-sided 0.7586241
# 2    Shannon -0.10762427 -0.9242510 two-sided 0.3639636
# 3 InvSimpson -0.03287402 -0.3736161 two-sided 0.7109289
# 4      Chao1  0.02745115  0.2478070 two-sided 0.8068193

##Pagel
###observed
trait <- dat[, 1]
names(trait) <- rownames(dat)

pagel.o <- phylosig(tree,
                    trait[tree$tip.label],
                    method = "lambda",
                    test = TRUE,
                    nsim = 10000)

pagel.o

#Phylogenetic signal lambda : 7.40706e-05 
#logL(lambda) : -140.526 
#LR(lambda=0) : -0.000801728 
#P-value (based on LR test) : 1 

###Shannon
trait <- dat[, 2]
names(trait) <- rownames(dat)

pagel.s <- phylosig(tree,
                    trait[tree$tip.label],
                    method = "lambda",
                    test = TRUE,
                    nsim = 10000)

pagel.s

# Phylogenetic signal lambda : 7.40706e-05 
# logL(lambda) : -38.2551 
# LR(lambda=0) : -0.000857602 
# P-value (based on LR test) : 1 

###Inverse Simpson
trait <- dat[, 3]
names(trait) <- rownames(dat)

pagel.is <- phylosig(tree,
                     trait[tree$tip.label],
                     method = "lambda",
                     test = TRUE,
                     nsim = 10000)

pagel.is

# Phylogenetic signal lambda : 7.40706e-05 
# logL(lambda) : -91.8281 
# LR(lambda=0) : -0.000613954 
# P-value (based on LR test) : 1 

###Chao1
trait <- dat[, 4]
names(trait) <- rownames(dat)

pagel.c1 <- phylosig(tree,
                     trait[tree$tip.label],
                     method = "lambda",
                     test = TRUE,
                     nsim = 10000)

pagel.c1

# Phylogenetic signal lambda : 7.40706e-05 
# logL(lambda) : -155.837 
# LR(lambda=0) : -0.00060123 
# P-value (based on LR test) : 1 


###animal symbionts
er.a.r <-estimate_richness(anm.sym.physeq.r, measures=c("Observed", "InvSimpson", "Shannon", "Chao1", "ACE"))%>%
  rownames_to_column(var = "guid")%>%
  left_join(all_meta)%>%
  group_by(sci_name)%>%
  dplyr::summarise(Observed = mean(Observed, na.rm=T),
                   Shannon = mean(Shannon, na.rm=T),
                   Chao1 = mean(Chao1, na.rm=T),
                   InvSimpson = mean(InvSimpson, na.rm=T))



dat <- er.a.r%>%
  dplyr::select(sci_name, Observed, Shannon,InvSimpson, Chao1)%>%
  mutate(sci_name = gsub(' ', "_", sci_name))%>%
  column_to_rownames(., var="sci_name")

p4d <- phylo4d(tree, dat[tree$tip.label,])
barplot(p4d, tree.ladderize =T, tree.ratio = 0.33, grid.vertical = F, show.box = T)


#Test phylosignal
## Moran's I
moran.test <- abouheif.moran(p4d,method="Abouheif", nrepet = 10000, alter = "two-sided")

moran.test
plot(moran.test)

# Adjustment method for multiple comparisons:   none 
# Permutation number:   10000 
# Test         Obs     Std.Obs     Alter    Pvalue
# 1   Observed -0.11203034 -0.71943221 two-sided 0.4635536
# 2    Shannon -0.03208905 -0.01764546 two-sided 0.9874013
# 3 InvSimpson -0.06517691 -0.32430931 two-sided 0.7501250
# 4      Chao1 -0.10370525 -0.64716116 two-sided 0.5142486

##Abouheifs Cmean
abhouief.test <- abouheif.moran(p4d,method="oriAbouheif", nrepet = 10000, alter = "two-sided")
abhouief.test
# Adjustment method for multiple comparisons:   none 
# Permutation number:   10000 
# Test          Obs     Std.Obs     Alter    Pvalue
# 1   Observed -0.087651260 -0.76185291 two-sided 0.4406559
# 2    Shannon  0.001756161  0.02487417 two-sided 0.9815018
# 3 InvSimpson -0.038669641 -0.34354218 two-sided 0.7312269
# 4      Chao1 -0.079122286 -0.68752922 two-sided 0.4864514

##Pagel
###observed
trait <- dat[, 1]
names(trait) <- rownames(dat)

pagel.o <- phylosig(tree,
                trait[tree$tip.label],
                method = "lambda",
                test = TRUE,
                nsim = 10000)

pagel.o

# Phylogenetic signal lambda : 1.00994 
# logL(lambda) : -93.9775 
# LR(lambda=0) : 2.72924 
# P-value (based on LR test) : 0.098526 

###Shannon
trait <- dat[, 2]
names(trait) <- rownames(dat)

pagel.s <- phylosig(tree,
                    trait[tree$tip.label],
                    method = "lambda",
                    test = TRUE,
                    nsim = 10000)

pagel.s

#Phylogenetic signal lambda : 7.40706e-05 
#logL(lambda) : -33.9123 
#LR(lambda=0) : -0.0004728 
#P-value (based on LR test) : 1 

###Inverse Simpson
trait <- dat[, 3]
names(trait) <- rownames(dat)

pagel.is <- phylosig(tree,
                    trait[tree$tip.label],
                    method = "lambda",
                    test = TRUE,
                    nsim = 10000)

pagel.is

#Phylogenetic signal lambda : 0.995502 
#logL(lambda) : -69.3454 
#LR(lambda=0) : 4.23102 
#P-value (based on LR test) : 0.0396917 

###Chao1
trait <- dat[, 4]
names(trait) <- rownames(dat)

pagel.c1 <- phylosig(tree,
                     trait[tree$tip.label],
                     method = "lambda",
                     test = TRUE,
                     nsim = 10000)

pagel.c1

#Phylogenetic signal lambda : 7.40706e-05 
#logL(lambda) : -106.745 
#LR(lambda=0) : -0.000668286 
#P-value (based on LR test) : 1 




###non-animal symbionts
er.na.r <-estimate_richness(not.anm.sym.physeq.MS.r, measures=c("Observed", "InvSimpson", "Shannon", "Chao1", "ACE"))%>%
  rownames_to_column(var = "guid")%>%
  left_join(all_meta)  

phylo.lm.df <- er.na.r

dat <- phylo.lm.df%>%
  dplyr::select(guid, Observed, Shannon,InvSimpson, Chao1)%>%
  mutate(guid = gsub('\\.', "_", guid))%>%
  column_to_rownames(., var="guid")

p4d <- phylo4d(tree, dat[tree$tip.label,])
barplot(p4d, tree.ladderize =T, tree.ratio = 0.33, grid.vertical = F, show.box = T)


#Test phylosignal
## Moran's I
moran.test <- abouheif.moran(p4d,method="Abouheif", nrepet = 10000)

moran.test
plot(moran.test)

#       Test        Obs  Std.Obs   Alter     Pvalue
# 1   Observed 0.31746862 4.320256 greater 0.01249875
# 2    Shannon 0.09838103 1.126797 greater 0.13038696
# 3 InvSimpson 0.30846063 3.358918 greater 0.01939806
# 4      Chao1 0.30067951 3.078854 greater 0.01479852

##Abouheifs Cmean
abhouief.test <- abouheif.moran(p4d,method="oriAbouheif", nrepet = 10000)
abhouief.test

#       Test       Obs  Std.Obs   Alter     Pvalue
# 1   Observed 0.3292629 3.730388 greater 0.00999900
# 2    Shannon 0.1156806 1.020393 greater 0.15688431
# 3 InvSimpson 0.3195391 3.083084 greater 0.01869813
# 4      Chao1 0.3119237 2.939826 greater 0.01579842

##Pagel
###observed
trait <- dat[, 1]
names(trait) <- rownames(dat)

pagel.o <- phylosig(tree,
                    trait[tree$tip.label],
                    method = "lambda",
                    test = TRUE,
                    nsim = 10000)

pagel.o

# Phylogenetic signal lambda : 0.660072 
# logL(lambda) : -148.775 
# LR(lambda=0) : 7.95021 
# P-value (based on LR test) : 0.00480817 

###Shannon
trait <- dat[, 2]
names(trait) <- rownames(dat)

pagel.s <- phylosig(tree,
                    trait[tree$tip.label],
                    method = "lambda",
                    test = TRUE,
                    nsim = 10000)

pagel.s

# Phylogenetic signal lambda : 0.541971 
# logL(lambda) : -42.785 
# LR(lambda=0) : 1.48728 
# P-value (based on LR test) : 0.222639 

###Inverse Simpson
trait <- dat[, 3]
names(trait) <- rownames(dat)

pagel.is <- phylosig(tree,
                     trait[tree$tip.label],
                     method = "lambda",
                     test = TRUE,
                     nsim = 10000)

pagel.is

#Phylogenetic signal lambda : 0.987917 
#logL(lambda) : -93.6934 
#LR(lambda=0) : 18.2556 
#P-value (based on LR test) : 1.9316e-05 

###Chao1
trait <- dat[, 4]
names(trait) <- rownames(dat)

pagel.c1 <- phylosig(tree,
                     trait[tree$tip.label],
                     method = "lambda",
                     test = TRUE,
                     nsim = 10000)

pagel.c1

#Phylogenetic signal lambda : 0.91593 
#logL(lambda) : -155.558 
#LR(lambda=0) : 11.3185 
#P-value (based on LR test) : 0.000767384 



#Phylogenetic model

# Bayesian phylogenetic models

#First, do we need phylogenetic models? Revell 2010 says we need to test for signal in residuals of Y~X but instead most folks text for signal in X and Y independently
dat.rs <- phylo.lm.df%>%
  dplyr::select(sci_name, Observed, Shannon,InvSimpson, Chao1,  MDS1, MDS2, Mass, Hand.Wing.Index, Max.Latitude, Range.Size, Trophic.Level, Trophic.Niche, Migration, Order1, sex, weight)%>%
  group_by(sci_name)%>%
  slice_sample(n=1)%>% #randomly sample single individuals of each subspecies
  ungroup()%>%
  mutate(Mass = ifelse(sci_name =="Antigone canadensis canadensis", 3400, ifelse(sci_name =="Antigone canadensis tabida", 5300, Mass)),
         Max.Latitude = ifelse(sci_name =="Antigone canadensis canadensis", 75, ifelse(sci_name =="Antigone canadensis tabida", 45, Max.Latitude)),
         Migration = as.factor(ifelse(sci_name =="Antigone canadensis canadensis", 3, ifelse(sci_name =="Antigone canadensis tabida", 2, Migration))),
                               weight = if_else(is.na(weight), Mass, weight),
         Observed = round(Observed),
         Hand.Wing.Index.z = scale(Hand.Wing.Index),
         Mass.z = scale(Mass),
         weight.z = scale(weight),
         Max.Lat.z = scale(Max.Latitude),
         Range.size.z = scale(Range.Size),
         species = gsub( " ", "_", sci_name, fixed=T))


#test residuals
lm.sh <- lm(Shannon ~ weight.z + Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z, data=dat.rs)


res.sh <- as.vector(lm.sh$residuals)
names(res.sh) <- dat.rs$sci_name
res.sh

#Test phylosignal
p4d <- phylo4d(tree, res.sh[tree$tip.label])

moran.test <- abouheif.moran(p4d,method="Abouheif")

moran.test
plot(moran.test)

abouheif.test <- abouheif.moran(p4d,method="oriAbouheif")

abouheif.test

plot(abouheif.test)


pagel.test <- phylosig(tree,
                 res.sh[tree$tip.label],
                 method = "lambda",
                 test = TRUE,
                 nsim = 10000)
pagel.test

plot(pagel.test)

#Shannon Diversity model has no phylogenetic signal

#have to do this with separate individuals
#df <- read_csv("data/phylo_df_unmerged.csv") #this has seperate individual community metrics

dat<- phylo.lm.df%>%
  dplyr::select(sci_name, Observed, Shannon,InvSimpson, Chao1,  MDS1, MDS2, Mass, Hand.Wing.Index, Max.Latitude, Range.Size, Trophic.Level, Trophic.Niche, Migration, Order1, sex, weight)%>%
  mutate(Mass = ifelse(sci_name =="Antigone canadensis canadensis", 3400, ifelse(sci_name =="Antigone canadensis tabida", 5300, Mass)),
         Max.Latitude = ifelse(sci_name =="Antigone canadensis canadensis", 75, ifelse(sci_name =="Antigone canadensis tabida", 45, Max.Latitude)),
         Migration = as.factor(ifelse(sci_name =="Antigone canadensis canadensis", 3, ifelse(sci_name =="Antigone canadensis tabida", 2, Migration))),
         weight = if_else(is.na(weight), Mass, weight),
         Observed = round(Observed),
         Hand.Wing.Index.z = scale(Hand.Wing.Index),
         Mass.z = scale(Mass),
         weight.z = scale(weight),
         Max.Lat.z = scale(Max.Latitude),
         Range.size.z = scale(Range.Size),
         species = gsub( " ", "_", sci_name, fixed=T))

dat$sci_name<-gsub("Antigone_canadensis_canadensis", "Antigone_canadensis canadensis", dat$sci_name)
dat$sci_name<-gsub("Antigone_canadensis_tabida", "Antigone_canadensis tabida", dat$sci_name)


#Tests of phylogenetic signal showed non on residuals or on Y for this model. Running without phylogeny
m1.g <- brm(
  formula = brms::bf(Shannon ~ 1 + weight.z + sex + Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z + (1|species)),
  data = dat,
  family = student(),
  cores = 4,
  chains = 4,
  thin = 10, #chop out transitions?
  warmup = 20000, #half of iterations
  iter = 40000,
  #prior = c(
  #     prior(normal(0, 10), "Intercept"),
  #     prior(cauchy(0, 2.5), "sigma")),
  
  prior = c(
    prior(student_t(3, 1.3, 2.5), "Intercept"),
    prior(student_t(3, 0, 20), "sigma"),
    prior(student_t(3, 0, 2.5), class ="b"),
    prior(cauchy(0,10), class = "sd")),
  save_pars = save_pars(all = TRUE), #need this for loo comparison
  control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
  set.seed(2020)
)

m1.g

conditional_effects(m1.g)

mcmc_plot(m1.g, variable = c("b_weight.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), prob_outer = 0.75)+
  scale_y_discrete(breaks =c("b_weight.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), labels = c("mass", "hand-wing index", "herbivore","omnivore", "partial migrant", "long distance migrant", "max latitude", "range size"))


m1 <- brm(
  formula = brms::bf(Shannon ~ 1 + weight.z + sex+ Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z + (1|species)),
  data = dat,
  family = Gamma(link="log"),
  cores = 8,
  chains = 4,
  thin = 10, #chop out transitions?
  warmup = 20000, #half of iterations
  iter = 40000,
  #prior = c(
  #     prior(normal(0, 10), "Intercept"),
  #     prior(cauchy(0, 2.5), "sigma")),
  
  prior = c(
    prior(student_t(3, 1.3, 2.5), "Intercept"),
    prior(student_t(3, 0, 2.5), class ="b"),
    prior(gamma(0.01,0.01),class="shape")),
  save_pars = save_pars(all = TRUE), #need this for loo comparison
  control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
  set.seed(2020)
)


#Test for overdispersion                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "b_Order1Cuculiformes", "b_Order1Falconiformes", "b_Order1Gaviiformes", "b_Order1Passeriformes", "b_Order1Pelecaniformes", "b_Order1Piciformes", "b_Order1Strigiformes"
rd <- glm(Observed ~ ., data = dat, family = poisson)
dispersiontest(rd,trafo=1)
#not overdispersed



#First, do we need phylogenetic models? Revell 2010 says we need to test for signal in residuals of Y~X but instead most folks text for signal in X and Y independently

#test residuals
lm.sh <- glm(Observed ~ weight.z + Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z, data=dat.rs, family = gaussian)

res.sh <- as.vector(lm.sh$residuals)
names(res.sh) <- dat.rs$sci_name
res.sh


#Test phylosignal
p4d <- phylo4d(tree, res.sh[tree$tip.label])

moran.test <- abouheif.moran(p4d,method="Abouheif")

moran.test
plot(moran.test)

abouheif.test <- abouheif.moran(p4d,method="oriAbouheif")

abouheif.test

plot(abouheif.test)


pagel.test <- phylosig(tree,
                       res.sh[tree$tip.label],
                       method = "lambda",
                       test = TRUE,
                       nsim = 10000)
pagel.test

plot(pagel.test)



m2 <- brm(
  formula = brms::bf(Observed ~ 1 + weight.z + sex +Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z + (1|species)),
  data = dat,
  #data2 = list(v = v),
  family = poisson(),
  cores = 8,
  chains = 4,
  thin = 10, #chop out transitions?
  warmup = 20000, #half of iterations
  iter = 40000,
  #prior = c(
  #     prior(normal(0, 10), "Intercept"),
  #     prior(cauchy(0, 2.5), "sigma")),
  
  prior = c(
    prior(student_t(3, 1.3, 2.5), "Intercept"),
    #prior(student_t(3, 0, 20), "sigma"),
    prior(student_t(3, 0, 2.5), class ="b"),
    prior(student_t(3, 0, 20), class = "sd")),
  save_pars = save_pars(all = TRUE), #need this for loo comparison
  control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
  set.seed(2020)
)

conditional_effects(m2)

mcmc_plot(m2, variable = c("b_weight.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), prob_outer = 0.95)+
  scale_y_discrete(breaks =c("b_weight.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), labels = c("mass", "hand-wing index", "herbivore","omnivore", "partial migrant", "long distance migrant", "max latitude", "range size"))+
  theme(aspect.ratio = 1)

m2.nb <- brm(
  formula = brms::bf(Observed ~ 1  + weight.z + sex + Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z +  (1|species)),
  data = dat,
  #data2 = list(v = v),
  family = negbinomial(),
  cores = 8,
  chains = 4,
  thin = 10, #chop out transitions?
  warmup = 20000, #half of iterations
  iter = 40000,
  #prior = c(
  #     prior(normal(0, 10), "Intercept"),
  #     prior(cauchy(0, 2.5), "sigma")),
  
  prior = c(
    prior(student_t(3, 1.3, 2.5), "Intercept"),
    #prior(student_t(3, 0, 20), "sigma"),
    prior(student_t(3, 0, 2.5), class ="b"),
    prior(student_t(3, 0, 20), class = "sd")),
  save_pars = save_pars(all = TRUE), #need this for loo comparison
  control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
  set.seed(2020)
)

conditional_effects(m2.nb)

mcmc_plot(m2.nb, variable = c("b_weight.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), prob_outer = 0.95)+
  scale_y_discrete(breaks =c("b_weight.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), labels = c("mass", "hand-wing index", "herbivore","omnivore", "partial migrant", "long distance migrant", "max latitude", "range size"))


m2.g <- brm(
  formula = brms::bf(Observed ~ 1  + Mass.z + Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z +  (1|species)),
  data = dat,
 # data2 = list(v = v),
  family = student(),
  cores = 8,
  chains = 4,
  thin = 10, #chop out transitions?
  warmup = 20000, #half of iterations
  iter = 40000,
  #prior = c(
  #     prior(normal(0, 10), "Intercept"),
  #     prior(cauchy(0, 2.5), "sigma")),
  
  prior = c(
    prior(student_t(3, 1.3, 2.5), "Intercept"),
    prior(student_t(3, 0, 20), "sigma"),
    prior(student_t(3, 0, 2.5), class ="b"),
    prior(student_t(3, 0, 20), class = "sd")),
  save_pars = save_pars(all = TRUE), #need this for loo comparison
  control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
  set.seed(2020)
)

m2.g
conditional_effects(m2.g)


mcmc_plot(m2.g, variable = c("b_Mass.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), prob_outer = 0.95)+
  scale_y_discrete(breaks =c("b_Mass.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), labels = c("mass", "hand-wing index", "herbivore","omnivore", "partial migrant", "long distance migrant", "max latitude", "range size"))



#First, do we need phylogenetic models? Revell 2010 says we need to test for signal in residuals of Y~X but instead most folks text for signal in X and Y independently


#test residuals
lm.sh <- lm(Chao1 ~ weight.z + Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z, data=dat.rs)


res.sh <- as.vector(lm.sh$residuals)
names(res.sh) <- dat.rs$sci_name
res.sh

#Test phylosignal
p4d <- phylo4d(tree, res.sh[tree$tip.label])

moran.test <- abouheif.moran(p4d,method="Abouheif")

moran.test
plot(moran.test)

abouheif.test <- abouheif.moran(p4d,method="oriAbouheif")

abouheif.test

plot(abouheif.test)


pagel.test <- phylosig(tree,
                       res.sh[tree$tip.label],
                       method = "lambda",
                       test = TRUE,
                       nsim = 10000)
pagel.test

plot(pagel.test)


#again no signal in residuals so ok to remove
m3 <- brm(
  formula = brms::bf(Chao1 ~ 1  + weight.z + sex+ Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z + (1|species)),
  data = dat,
  #data2 = list(v = v),
  family = Gamma(link="log"),
  cores = 8,
  chains = 4,
  thin = 10, #chop out transitions?
  warmup = 20000, #half of iterations
  iter = 40000,
  #prior = c(
  #     prior(normal(0, 10), "Intercept"),
  #     prior(cauchy(0, 2.5), "sigma")),
  init = 0,
  prior = c(
    prior(normal(0,2),class="Intercept"),
    prior(normal(0,2),class="b"),
    prior(gamma(0.01,0.01),class="shape")),
  save_pars = save_pars(all = TRUE), #need this for loo comparison
  control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
  set.seed(2020)
)

conditional_effects(m3)
mcmc_plot(m3, variable = c("b_weight.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), prob_outer = 0.95)+
  scale_y_discrete(breaks =c("b_weight.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), labels = c("mass", "hand-wing index", "herbivore","omnivore", "partial migrant", "long distance migrant", "max latitude", "range size"))

#test residuals
lm.sh <- lm(MDS1 ~ Mass.z + Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z, data=dat.rs)


res.sh <- as.vector(lm.sh$residuals)
names(res.sh) <- dat.rs$sci_name
res.sh

#Test phylosignal
p4d <- phylo4d(tree, res.sh[tree$tip.label])

moran.test <- abouheif.moran(p4d,method="Abouheif")

moran.test
plot(moran.test)

abouheif.test <- abouheif.moran(p4d,method="oriAbouheif")

abouheif.test

plot(abouheif.test)


pagel.test <- phylosig(tree,
                       res.sh[tree$tip.label],
                       method = "lambda",
                       test = TRUE,
                       nsim = 10000)
pagel.test

plot(pagel.test)


m4 <- brm(
  formula = brms::bf(MDS1 ~ 1  + weight.z + Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z + (1|species)),
  data = dat,
  #data2 = list(v = v),
  family = student(),
  cores = 8,
  chains = 4,
  thin = 10, #chop out transitions?
  warmup = 20000, #half of iterations
  iter = 40000,
  #prior = c(
  #     prior(normal(0, 10), "Intercept"),
  #     prior(cauchy(0, 2.5), "sigma")),
  
  prior = c(
    prior(student_t(3, 1.3, 2.5), "Intercept"),
    prior(student_t(3, 0, 20), "sigma"),
    prior(student_t(3, 0, 2.5), class ="b"),
    prior(student_t(3, 0, 20), class = "sd")),
  save_pars = save_pars(all = TRUE), #need this for loo comparison
  control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
  set.seed(2020)
)

conditional_effects(m4)
mcmc_plot(m4, variable = c("b_weight.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), prob_outer = 0.75)+
  scale_y_discrete(breaks =c("b_weight.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), labels = c("mass", "hand-wing index", "herbivore","omnivore", "partial migrant", "long distance migrant", "max latitude", "range size"))


#test residuals
lm.sh <- lm(MDS2 ~ Mass.z + Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z, data=dat.rs)


res.sh <- as.vector(lm.sh$residuals)
names(res.sh) <- dat.rs$sci_name
res.sh

#Test phylosignal
p4d <- phylo4d(tree, res.sh[tree$tip.label])

moran.test <- abouheif.moran(p4d,method="Abouheif")

moran.test
plot(moran.test)

abouheif.test <- abouheif.moran(p4d,method="oriAbouheif")

abouheif.test

plot(abouheif.test)


pagel.test <- phylosig(tree,
                       res.sh[tree$tip.label],
                       method = "lambda",
                       test = TRUE,
                       nsim = 10000)
pagel.test

plot(pagel.test)


m5 <- brm(
  formula = brms::bf(MDS2 ~ 1  + weight.z + Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z + (1|species)),
  data = dat,
  #data2 = list(v = v), #no phylosignal in residuals
  family = student(),
  cores = 8,
  chains = 4,
  thin = 10, #chop out transitions?
  warmup = 20000, #half of iterations
  iter = 40000,
  #prior = c(
  #     prior(normal(0, 10), "Intercept"),
  #     prior(cauchy(0, 2.5), "sigma")),
  
  prior = c(
    prior(student_t(3, 1.3, 2.5), "Intercept"),
    prior(student_t(3, 0, 20), "sigma"),
    prior(student_t(3, 0, 2.5), class ="b"),
    prior(cauchy(0,10), class = "sd")),
  save_pars = save_pars(all = TRUE), #need this for loo comparison
  control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
  set.seed(2020)
)

conditional_effects(m5)
mcmc_plot(m5, variable = c("b_weight.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), prob_outer = 0.75)+
  scale_y_discrete(breaks =c("b_weight.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), labels = c("mass", "hand-wing index", "herbivore","omnivore", "partial migrant", "long distance migrant", "max latitude", "range size"))


m6 <- brm(
  formula = brms::bf(InvSimpson ~ 1  + Mass.z + Hand.Wing.Index.z  + Trophic.Level + Migration + Range.size.z + Max.Lat.z + (1|species) + (1 | gr(sci_name, cov = v))),
  data = dat,
  data2 = list(v = v),
  family = skew_normal(),
  cores = 8,
  chains = 4,
  thin = 10, #chop out transitions?
  warmup = 20000, #half of iterations
  iter = 40000,
  #prior = c(
  #     prior(normal(0, 10), "Intercept"),
  #     prior(cauchy(0, 2.5), "sigma")),
  
  prior = c(
    prior(student_t(3, 1.3, 2.5), "Intercept"),
    prior(student_t(3, 0, 20), "sigma"),
    prior(student_t(3, 0, 2.5), class ="b"),
    prior(student_t(3, 0, 20), class = "sd")),
  save_pars = save_pars(all = TRUE), #need this for loo comparison
  control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
  set.seed(2020)
)

summary(m6)
conditional_effects(m6)
mcmc_plot(m6, variable = c("b_Mass.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"),  prob_outer = 0.95)+
  scale_y_discrete(breaks =c("b_Mass.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), labels = c("mass", "hand-wing index", "herbivore","omnivore", "partial migrant", "long distance migrant", "max latitude", "range size"))


m6.ln <- brm(
  formula = brms::bf(InvSimpson ~ 1  + Mass.z + Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z + (1|species) + (1 | gr(sci_name, cov = v))),
  data = dat,
  data2 = list(v = v),
  family = lognormal(),
  cores = 8,
  chains = 4,
  thin = 10, #chop out transitions?
  warmup = 20000, #half of iterations
  iter = 40000,
  #prior = c(
  #     prior(normal(0, 10), "Intercept"),
  #     prior(cauchy(0, 2.5), "sigma")),
  
  prior = c(
    prior(student_t(3, 1.3, 2.5), "Intercept"),
    #prior(student_t(3, 0, 20), "sigma"),
    prior(student_t(3, 0, 2.5), class ="b"),
    prior(student_t(3, 0, 20), class = "sd")),
  save_pars = save_pars(all = TRUE), #need this for loo comparison
  control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
  set.seed(2020)
)

summary(m6.ln)
conditional_effects(m6.ln)
mcmc_plot(m6.ln, variable = c("b_Mass.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"),  prob_outer = 0.95)+
  scale_y_discrete(breaks =c("b_Mass.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), labels = c("mass", "hand-wing index", "herbivore","omnivore", "partial migrant", "long distance migrant", "max latitude", "range size"))



m6.g <- brm(
  formula = brms::bf(InvSimpson ~ 1  + weight.z + sex + Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z + (1|species) + (1 | gr(sci_name, cov = v))),
  data = dat,
  data2 = list(v = v),
  family = student(),
  cores = 8,
  chains = 4,
  thin = 10, #chop out transitions?
  warmup = 20000, #half of iterations
  iter = 40000,
  #prior = c(
  #     prior(normal(0, 10), "Intercept"),
  #     prior(cauchy(0, 2.5), "sigma")),
  
  prior = c(
    prior(student_t(3, 1.3, 2.5), "Intercept"),
    prior(student_t(3, 0, 20), "sigma"),
    prior(student_t(3, 0, 2.5), class ="b"),
    prior(student_t(3, 0, 20), class = "sd")),
  save_pars = save_pars(all = TRUE), #need this for loo comparison
  control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
  set.seed(2020)
)


summary(m6.g)
conditional_effects(m6.g)

mcmc_plot(m6.g, variable = c("b_weight.z", "b_sexmale", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"),  prob_outer = 0.75)+
  scale_y_discrete(breaks =c("b_weight.z", "b_Hand.Wing.Index.z", "b_Trophic.LevelHerbivore","b_Trophic.LevelOmnivore", "b_Migration2", "b_Migration3", "b_Max.Lat.z", "b_Range.size.z"), labels = c("mass", "hand-wing index", "herbivore","omnivore", "partial migrant", "long distance migrant", "max latitude", "range size"))



#PGLS


# Calculate the phylogenetic correlation matrix
v <- vcv(tree,  model = "Brownian", corr = TRUE)


dat<- phylo.lm.df%>%
  dplyr::select(sci_name, Observed, Shannon,InvSimpson, Chao1, Kipps.Distance, Mass, Wing.Length, Centroid.Latitude)%>%
  mutate(sci_name2 = gsub( " ", "_", sci_name, fixed=T))%>%
  column_to_rownames(., var="sci_name2")


# Fit the GLS model using the lm.gls() function
#z <- lm.gls(InvSimpson~Kipps.Distance, data=dat, v = v)


pglsModel <- gls(log10(InvSimpson) ~ log10(Kipps.Distance), correlation = corBrownian(phy = tree),
                 data = dat, method = "ML")
summary(pglsModel)

# Extract columns
is <- log10(dat[, "InvSimpson"])
kd <- log10(dat[, "Kipps.Distance"])
hist(kd)
hist(is)

plot(is ~ kd)
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])

#check residuals
plot(pglsModel, abline=c(0,0)) #better



pglsModel <- gls(log10(InvSimpson) ~ log10(Mass), correlation = corBrownian(phy = tree),
                 data = dat, method = "ML")
summary(pglsModel)

# Extract columns
is <- dat[, "InvSimpson"]
ms <- dat[, "Mass"]
hist(kd)
hist(ms)

plot(log10(is) ~ log10(ms))
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])

#check residuals
plot(pglsModel, abline=c(0,0)) #pretty good



pglsModel <- gls(NMDS1 ~ Kipps.Distance, correlation = corBrownian(phy = tree),
                 data = dat, method = "ML")
summary(pglsModel)

# Extract columns
sh <- dat[, "NMDS1"]
wl <- dat[, "Wing.Length"]
hist(sh)
hist(wl)

plot(log10(sh) ~ log10(wl))
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])

#check residuals
plot(pglsModel, abline=c(0,0)) #pretty good

tree2<-ladderize(tree)
ggplot(phylo.lm.df, aes(x=factor(sci_name, level = tree2$tip.label), y=Observed, color=Trophic.Level))+
  geom_point()+
  labs(y=NULL)+
  theme(axis.text.x = element_text(size=10),
        aspect.ratio = 0.5)
