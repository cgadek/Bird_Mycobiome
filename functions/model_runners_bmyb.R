# Run models functions####

# Shannon####

run_my_bmyc_shannon_models <- function(data, name, seed) {
  require(tidyverse)
  require(brms)
  
  if (is.null(seed))
    seed <- NA
  
f1 <- brm(
    formula = brms::bf(Shannon ~ 1),
    data = dat,
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
      prior(student_t(3, 1.3, 2.5), "Intercept")),
    save_pars = save_pars(all = TRUE), #need this for loo comparison
    control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
    set.seed(2020)
  )


f2 <- brm(
  formula = brms::bf(Shannon ~ 1 + (1|species)),
  data = dat,
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
    prior(student_t(3, 1.3, 2.5), "Intercept")),
  save_pars = save_pars(all = TRUE), #need this for loo comparison
  control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
  set.seed(2020)
)

f3 <- brm(
  formula = brms::bf(Shannon ~ 1 + weight.z + sex+ Hand.Wing.Index.z + Kipps.Distance.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z + (1|species)),
  data = dat,
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
    prior(student_t(3, 0, 2.5), class ="b")),
  save_pars = save_pars(all = TRUE), #need this for loo comparison
  control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
  set.seed(2020)
)
  

  #compare model with LOOIC
  m.comp<-LOO(f1,f2,f3,moment_match=T)
  
  #Get Bayes R^2
  br1<-bayes_R2(f1)
  br2<-bayes_R2(f2)
  br3<-bayes_R2(f3)
  
  br2L <- list(br1, br2, br3)
  
  t <- tibble(dataset = c(rep(name,times=3),"LOO_comp", "BayesR2"),
              model_set = c(1:3, "model_compare", "BayesR2"),
              m = list(f1, f2, f3, m.comp, br2L))
  
  assign(name, t, envir=globalenv())
} 

# Species richness####
run_my_bmyc_sprch_models(data=dat, name="species_rich_mods",  seed=2020)

run_my_bmyc_sprch_models <- function(data, name, seed) {
  require(tidyverse)
  require(brms)
  
  if (is.null(seed))
    seed <- NA
  
  f1 <- brm(
    formula = brms::bf(Observed ~ 1),
    data = dat,
    #data2 = list(v = v),
    family = poisson(link = log),
    cores = 8,
    chains = 4,
    thin = 10, #chop out transitions?
    warmup = 20000, #half of iterations
    iter = 40000,
    #prior = c(
    #     prior(normal(0, 10), "Intercept"),
    #     prior(cauchy(0, 2.5), "sigma")),
    
    prior = c(
      prior(student_t(3, 1.3, 2.5), "Intercept")),
    save_pars = save_pars(all = TRUE), #need this for loo comparison
    control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
    set.seed(2020)
  )
  
  
  f2 <- brm(
    formula = brms::bf(Observed ~ 1 + (1|species)),
    data = dat,
    #data2 = list(v = v),
    family = poisson(link = log),
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
      prior(normal(0,10), class = "sd")),
    save_pars = save_pars(all = TRUE), #need this for loo comparison
    control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
    set.seed(2020)
  )
  
  f3 <- brm(
    formula = brms::bf(Observed ~ 1 + weight.z + sex + Hand.Wing.Index.z + Kipps.Distance.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z + (1|species)),
    data = dat,
    #data2 = list(v = v),
    family = poisson(link = log),
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
      prior(normal(0,10), class = "sd")),
    save_pars = save_pars(all = TRUE), #need this for loo comparison
    control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
    set.seed(2020)
  )
  
  
  #compare model with LOOIC
  m.comp<-LOO(f1,f2,f3,moment_match=T)
  
  #Get Bayes R^2
  br1<-bayes_R2(f1)
  br2<-bayes_R2(f2)
  br3<-bayes_R2(f3)
  
  br2L <- list(br1, br2, br3)
  
  t <- tibble(dataset = c(rep(name,times=3),"LOO_comp", "BayesR2"),
              model_set = c(1:3, "model_compare", "BayesR2"),
              m = list(f1, f2, f3, m.comp, br2L))
  
  assign(name, t, envir=globalenv())
}

#Chao####
run_my_bmyc_chao_models(data=dat, name="shannon_mods",  seed=2020)
run_my_bmyc_chao_models <- function(data, name, seed) {
  require(tidyverse)
  require(brms)
  
  if (is.null(seed))
    seed <- NA
  
  f1 <- brm(
    formula = brms::bf(Chao1 ~ 1),
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
      prior(gamma(0.01,0.01),class="shape")),
    save_pars = save_pars(all = TRUE), #need this for loo comparison
    control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
    set.seed(2020)
  )
  
  f2 <- brm(
    formula = brms::bf(Chao1 ~ 1  +  (1|species)),
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
      prior(gamma(0.01,0.01),class="shape")),
    save_pars = save_pars(all = TRUE), #need this for loo comparison
    control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
    set.seed(2020)
  )
  
  f3 <- brm(
    formula = brms::bf(Chao1 ~ 1  + weight.z + sex+ Hand.Wing.Index.z + Kipps.Distance.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z + (1|species)),
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
  
  
  #compare model with LOOIC
  m.comp<-LOO(f1,f2,f3,moment_match=T)
  
  #Get Bayes R^2
  br1<-bayes_R2(f1)
  br2<-bayes_R2(f2)
  br3<-bayes_R2(f3)
  
  br2L <- list(br1, br2, br3)
  
  t <- tibble(dataset = c(rep(name,times=3),"LOO_comp", "BayesR2"),
              model_set = c(1:3, "model_compare", "BayesR2"),
              m = list(f1, f2, f3, m.comp, br2L))
  
  assign(name, t, envir=globalenv())
}

# Inverse Simpson####

run_my_bmyc_simp_models(data=dat, name="simpson_mods",  seed=2020)
run_my_bmyc_simp_models <- function(data, name, seed) {
  require(tidyverse)
  require(brms)
  
  if (is.null(seed))
    seed <- NA
  
  f1 <- brm(
    formula = brms::bf(InvSimpson ~ 1),
    data = dat,
    #data2 = list(v = v),
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
      prior(student_t(3, 1.3, 2.5), "Intercept")),
    save_pars = save_pars(all = TRUE), #need this for loo comparison
    control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
    set.seed(2020)
  )
  
  f2 <- brm(
    formula = brms::bf(InvSimpson ~ 1 + (1|species)),
    data = dat,
    #data2 = list(v = v),
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
      prior(student_t(3, 0, 20), class = "sd")),
    save_pars = save_pars(all = TRUE), #need this for loo comparison
    control = list(adapt_delta = 0.995, max_treedepth = 19.5), #adapt_delta the target average proposal acceptance probability during Stan's adaptation period
    set.seed(2020)
  )
  
  f3 <- brm(
    formula = brms::bf(InvSimpson ~ 1  + weight.z + Kipps.Distance.z + Hand.Wing.Index.z + Trophic.Level + Migration + Range.size.z + Max.Lat.z + (1|species)),
    data = dat,
    #data2 = list(v = v),
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
  
  #compare model with LOOIC
  m.comp<-LOO(f1,f2,f3,moment_match=T)
  
  #Get Bayes R^2
  br1<-bayes_R2(f1)
  br2<-bayes_R2(f2)
  br3<-bayes_R2(f3)
  
  br2L <- list(br1, br2, br3)
  
  t <- tibble(dataset = c(rep(name,times=3),"LOO_comp", "BayesR2"),
              model_set = c(1:3, "model_compare", "BayesR2"),
              m = list(f1, f2, f3, m.comp, br2L))
  
  assign(name, t, envir=globalenv())
}


