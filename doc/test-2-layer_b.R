# remotes::install_github("andrewdolman/baconr", ref = "hierarchical-acc.mean", args = "--preclean", build_vignettes = FALSE)
# or
# devtools::install_github("andrewdolman/baconr", ref = "hierarchical-acc.mean", args = "--preclean", build_vignettes = FALSE)

library(rstan)
library(baconr)
library(tidyverse)

# pkgbuild::compile_dll()
# devtools::load_all()
# devtools::install(quick=FALSE)


## New make stan dat
insertLayer <- function(P, after=0, ...) {
  #  P     : Plot object
  # after  : Position where to insert new layers, relative to existing layers
  #  ...   : additional layers, separated by commas (,) instead of plus sign (+)
  
  if (after < 0)
    after <- after + length(P$layers)
  
  if (!length(P$layers))
    P$layers <- list(...)
  else
    P$layers <- append(P$layers, list(...), after)
  
  return(P)
}


name <- "CARPLAKE"
#dat <- read.csv(paste0("/Users/andrewdolman/Dropbox/Work/AWI/Data/terrestrial-age-models/terr_14C_min10_dates-2020.03.04_15-19-42/", name, "/", name,".csv")) %>%
dat <- read.csv(paste0("../envi-age-modelling/working-data/terr_14C_min10_dates-2020.03.04_15-19-42/", name, "/", name,".csv")) %>%
  filter(depth > 0)

dat <- ecustools::CalibrateAge(dat, age.14C = "age", age.14C.se = "error") %>%
  filter(complete.cases(age.14C.cal))

dat %>%
  ggplot(aes(x = depth, y = age.14C.cal)) +
  geom_point() +
  geom_line()

dat %>% 
  select(depth, age.14C.cal, age.14C.cal.se) %>% 
  mutate(age.14C.cal = round(age.14C.cal),
         age.14C.cal.se = round(age.14C.cal.se)) %>% 
  dput()

dat1 <- structure(list(depth = c(0, 185, 243.5, 262, 304, 321, 435.5, 485, 
                                 515, 556, 605, 652, 705, 775),
                       age.14C.cal = c(0, 6620, 9810, 10789,  9249, 11201, 
                                       19399, 22055, 25878, 25302, 25239, 
                                       22230, 30470, 36844),
                       age.14C.cal.se = c(50, 64, 162, 183, 146, 601, 451, 125, 
                                          368, 437, 440, 212, 232, 611)),
                  class = "data.frame", row.names = c(NA, -14L))


dat1 %>%
  ggplot(aes(x = depth, y = age.14C.cal)) +
  geom_point() +
  geom_line()

dat1 <- dat1 %>% 
  arrange(age.14C.cal.se)

make_stan_dat_sarnie(depth = dat1$depth, obs_age = dat1$age.14C.cal, obs_err = dat1$age.14C.cal.se)


acc.mean <- coef(MASS::rlm(age.14C.cal~depth, data = dat1))[2]

acc.mean

options(mc.cores = parallel::detectCores())

bacon.fit1 <- stan_bacon(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  K = 100,
  nu = 6,
  acc_mean = 40,
  mem_mean = 0.7, mem_strength = 4,
  chains = 3)

plot_stan_bacon(bacon.fit1, n.iter = 100, plot_priors = F)


sarn.fit1 <- sarnie(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  #bottom_depth = 800,
  K = 100,
  K1 = 1,
  nu = 6,
  record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 0, chains = 3)



plot_stan_bacon(sarn.fit1, n.iter = 100, plot_priors = T)


sarn.fit2 <- sarnie(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  K = 150,
  K1 = 21,
  nu = 6,
  record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 0, chains = 3)

p2 <- plot_stan_bacon(sarn.fit2, n.iter = 100, plot_priors = F)
p2


summary(sarn.fit2$fit, par = c("record_acc_mean", "record_acc_shape"))$summary
traceplot(sarn.fit2$fit, par = c("record_acc_mean", "record_acc_shape"))

#shinystan::launch_shinystan(sarn.fit3$fit)


sarn.fit3 <- sarnie(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  K = 100,
  K1 = 10,
  nu = 6,
  record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 1, chains = 3)


p.single.level <- plot_stan_bacon(sarn.fit1, 1000, plot_priors = T)
p.multi.level <- plot_stan_bacon(sarn.fit2, 1000, plot_priors = T)
p.multi.level.infl <- plot_stan_bacon(sarn.fit3, 1000, plot_priors = T)





# # infl errors
# 
# s3 <- summary(sarn.fit3$fit) 
# s3 <- as_tibble(s3$summary, rownames = "par")
# 
# infl <- s3 %>%
#   filter(grepl("infl[", par, fixed = T),
#          grepl("err_infl", par, fixed = T) == FALSE)
# 
# dat.infl <- tibble(depth = sarn.fit3$data$depth, obs_age = sarn.fit3$data$obs_age,
#                    obs_err = sarn.fit3$data$obs_err,
#                    infl = infl$mean) %>%
#   mutate(infl.err = obs_err + obs_err * infl)
# 

