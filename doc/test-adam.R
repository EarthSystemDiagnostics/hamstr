# test adam
library(tidyverse)
library(baconr)
library(rstan)

# devtools::document()
# pkgbuild::compile_dll(force = TRUE)
# devtools::load_all()
# devtools::install(quick=FALSE)


#name <- "BLOODMA"
#dat <- read.csv(paste0("/Users/andrewdolman/Dropbox/Work/AWI/Data/terrestrial-age-models/terr_14C_min10_dates-2020.03.04_15-19-42/", name, "/", name,".csv"))


#dat <- read.csv(paste0("../envi-age-modelling/working-data/terr_14C_min10_dates-2020.03.04_15-19-42/", name, "/", name,".csv")) %>%
#  filter(depth > 0)

#Load and filter data
all.terr.dat <- readr::read_csv2("doc/Dating_Data.csv") %>%
  filter(complete.cases(age, depth, e.older)) %>%
  mutate(sigma.age = e.older)

all.terr.14C.dat <- all.terr.dat %>%
  filter(age.type == "Radiocarbon years BP") %>%
  group_by(DataName) %>%
  mutate(n.dates = n()) %>%
  ungroup()

name <- "BUROVER"
#name <- "BLACKMA"
dat2 <- all.terr.14C.dat %>%
  filter(DataName == name)

dat1 <- ecustools::CalibrateAge(dat2, age.14C = "age", age.14C.se = "sigma.age") %>%
  filter(complete.cases(age.14C.cal)) %>%
  select(depth, age.14C.cal, age.14C.cal.se, age, sigma.age) %>%
  mutate(age.14C.cal = round(age.14C.cal),
         age.14C.cal.se = round(age.14C.cal.se))

dat1 %>%
  ggplot(aes(x = depth, y = age.14C.cal)) +
  geom_point() +
  geom_line()


options(mc.cores = parallel::detectCores())

adam.fit1 <- adam(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  K = c(3,3,3, 3, 3),
  nu = 6,
  #record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 0, chains = 3)


plot_adam(adam.fit1, plot_diagnostics = T)
plot_adam(adam.fit1, type = "spaghetti", n.iter = 1000, plot_diagnostics = T)



adam.fit2 <- adam(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  #K = c(5,5,5,5),
  K = baconr:::optimal_K(100, 10),
  nu = 6,
  #record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 0, chains = 3)



plot_adam(adam.fit2, plot_diagnostics = T)
plot_adam(adam.fit2, type = "spaghetti", n.iter = 1000, plot_diagnostics = T)


a2 <- rstan::summary(adam.fit2$fit)
a2 <- as_tibble(a2$summary, rownames = "par")


summary(adam.fit2$fit, par = c("alpha[1]", "shape"))$summary

traceplot(adam.fit2$fit, pars = c("shape", "alpha[1]"), inc_warmup = T)

traceplot(adam.fit2$fit, pars = paste0("alpha[", 91:111, "]"), inc_warmup = T)


#shinystan::launch_shinystan(adam.fit3$fit)


adam.fit3 <- adam(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  #top_depth = 1,
  K = baconr:::optimal_K(100, 10),
  nu = 6,
  #record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 1, chains = 3)

plot_adam(adam.fit3, type = "ribbon", plot_diagnostics = TRUE)
plot_adam(adam.fit3, type = "spaghetti", n.iter = 1000, plot_diagnostics = TRUE)


a3 <- rstan::summary(adam.fit3$fit)
a3 <- as_tibble(a3$summary, rownames = "par")

summary(adam.fit3$fit, par = c("alpha[1]", "shape"))$summary
summary(adam.fit3$fit, par = c("w"))$summary

traceplot(adam.fit3$fit, par = c("alpha[1]", "shape"))

traceplot(adam.fit3$fit, par = c("infl_mean", "infl_shape"))
traceplot(adam.fit3$fit, par = c("infl"))

summary(adam.fit3$fit, par = c("infl_mean", "infl_shape"))$summary


baconr:::plot_memory_prior_posterior(adam.fit3)



idx <- as_tibble(alpha_indices(c(10,10,10))[1:3]) %>%
  mutate(alpha_idx = as.character(alpha_idx))


alph <- a3 %>%
  filter(grepl("alpha", par, fixed = TRUE)) %>%
  separate(par, into = c("par", "alpha_idx")) %>%
  left_join(idx, .) %>%
  mutate(lvl = factor(lvl))


alph %>%
  ggplot(aes(x = mean)) +
  geom_histogram(aes()) +
  facet_wrap(~lvl, scales = "free_y")

alph %>%
  filter(lvl == 4) %>%
  ggplot(aes(x = mean)) +
  geom_histogram(aes(fill = as.factor(parent))) +
  facet_wrap(~lvl, scales = "free_y")


alph %>%
  ggplot(aes(x = Rhat)) +
  geom_histogram(aes(fill = lvl)) +
  expand_limits(x = c(0.95, 1.1)) +
  facet_wrap(~lvl, scales = "free_y")




plot_hierarchical_acc_rate(adam.fit2)







CompareAgePDF <- function(age.14C, age.14C.se, curve, t.df = 6, return.type = c("plot", "data")){

  dt_ls <- function(x, df=1, mu=0, sigma=1) 1/sigma * dt((x - mu)/sigma, df)

  cal.dat <- data.frame(age.14C = age.14C, age.14C.se = age.14C.se)

  calib <- ecustools::CalibrateAge(cal.dat,
                                   return.type = "lst",
                                   offset = 0, curve = curve)

  # The summarised calendar ages are appended to the input data
  C14 <- calib$df
  C14$.id <- 1:nrow(cal.dat)

  # These are the full PDFs
  cal.ages <- calib$cal.ages

  cali.pdf.df <- plyr::ldply(1:length(cal.ages), function(i){
    x <- cal.ages[[i]]
    if (is.na(x)==FALSE){data.frame(age = x[[1]]$ageGrid, d = x[[1]]$densities, .id = i)}else{
      data.frame(age = 0, d = 0, .id = i)
    }
  })

  t.dat <- C14 %>%
    group_by(.id) %>%
    do({
      rng <- .$age.14C.cal.se * 5
      age = seq(.$age.14C.cal - rng, .$age.14C.cal + rng, length.out = 100)
      data.frame(
        age = age,
        d = dt_ls(age, df = t.df, mu = .$age.14C.cal, sigma = .$age.14C.cal.se)
      )
    })

  gg <- cali.pdf.df %>%
    ggplot(aes(x = age/1000, y = d, group = .id)) +
    geom_line(aes(colour = "Empirical")) +
    geom_line(data = t.dat, aes(y = d, colour = "t-dist")) +
    labs(colour = "")

  if (return.type == "data"){
    return(list(cal.age.pdf = cali.pdf.df, t.dist.age = t.dat))
  } else if (return.type == "plot") {
    return(gg)
  }



}

p <- CompareAgePDF(dat1$age.14C.cal[1:6], dat1$sigma.age[1:6], curve = "intcal13", t.df = 6)
p + facet_wrap(~.id, scales = "free")


