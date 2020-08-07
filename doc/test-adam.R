# test adam
library(tidyverse)
library(baconr)
library(rstan)

# devtools::document()
# pkgbuild::compile_dll(force = TRUE)
# devtools::load_all()
# devtools::install(quick=FALSE)


name <- "BLACKMA"
#dat <- read.csv(paste0("/Users/andrewdolman/Dropbox/Work/AWI/Data/terrestrial-age-models/terr_14C_min10_dates-2020.03.04_15-19-42/", name, "/", name,".csv"))


dat <- read.csv(paste0("../envi-age-modelling/working-data/terr_14C_min10_dates-2020.03.04_15-19-42/", name, "/", name,".csv")) %>%
  filter(depth > 0)

dat <- ecustools::CalibrateAge(dat, age.14C = "age", age.14C.se = "error") %>%
  filter(complete.cases(age.14C.cal))

dat %>%
  ggplot(aes(x = depth, y = age.14C.cal)) +
  geom_point() +
  geom_line()

dat1 <- dat %>%
  select(depth, age.14C.cal, age.14C.cal.se) %>%
  mutate(age.14C.cal = round(age.14C.cal),
         age.14C.cal.se = round(age.14C.cal.se)) #%>%

dat1 %>%
  ggplot(aes(x = depth, y = age.14C.cal)) +
  geom_point() +
  geom_line()

dat1 <- dat1 %>%
  arrange(age.14C.cal.se)

stan_dat <- make_stan_dat_adam(depth = dat1$depth, obs_age = dat1$age.14C.cal, obs_err = dat1$age.14C.cal.se)


acc.mean <- coef(MASS::rlm(age.14C.cal~depth, data = dat1))[2]

acc.mean

options(mc.cores = parallel::detectCores())

# bacon.fit1 <- stan_bacon(
#   depth = dat1$depth,
#   obs_age = dat1$age.14C.cal,
#   obs_err = dat1$age.14C.cal.se,
#   K = 100,
#   nu = 6,
#   acc_mean = acc.mean,
#   mem_mean = 0.7, mem_strength = 4,
#   chains = 3)
# 
# plot_adam(bacon.fit1, n.iter = 100, plot_diagnostics = F)


adam.fit2 <- adam(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  #K = c(5,5,5,5),
  K = optimal_K(700, 10),
  nu = 6,
  #record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 0, chains = 3)



plot_adam(adam.fit2, plot_diagnostics = T)
plot_adam(adam.fit2, type = "spaghetti", n.iter = 100, plot_diagnostics = T)




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
  K = c(10, 10, 10),
  nu = 6,
  #record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 1, chains = 3)

plot_adam(adam.fit3, type = "ribbon", plot_diagnostics = TRUE)
plot_adam(adam.fit3, type = "spaghetti", n.iter = 100, plot_diagnostics = TRUE)


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









