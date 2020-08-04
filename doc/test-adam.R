# test adam
library(tidyverse)
library(baconr)
library(rstan)

#adam_mod <- stan_model(file = "inst/stan/adam.stan")

# stan_dat <- make_stan_dat_adam(depth = MSB2K$depth,
#                                obs_age = MSB2K$age,
#                                obs_err = MSB2K$error,
#                                #K1 = 10, K = 100,
#                                nu = 6)
#
# inits <- get_inits_adam(stan_dat)
#
# n.chains = 3
#
# inits <- rep(list(inits), n.chains)
#
#
# fit <- sampling(adam_mod,
#                 data = stan_dat, init = inits, chains = n.chains,
#                        verbose = FALSE)
#
#
# fit
#
# traceplot(fit, par = "alpha[1]")

##########

name <- "ALUT"
dat <- read.csv(paste0("/Users/andrewdolman/Dropbox/Work/AWI/Data/terrestrial-age-models/terr_14C_min10_dates-2020.03.04_15-19-42/", name, "/", name,".csv"))


#dat <- read.csv(paste0("../envi-age-modelling/working-data/terr_14C_min10_dates-2020.03.04_15-19-42/", name, "/", name,".csv")) %>%
#  filter(depth > 0)

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

make_stan_dat_adam(depth = dat1$depth, obs_age = dat1$age.14C.cal, obs_err = dat1$age.14C.cal.se)


acc.mean <- coef(MASS::rlm(age.14C.cal~depth, data = dat1))[2]

acc.mean

options(mc.cores = parallel::detectCores())

bacon.fit1 <- stan_bacon(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  K = 100,
  nu = 6,
  acc_mean = acc.mean,
  mem_mean = 0.7, mem_strength = 4,
  chains = 3)

plot_stan_bacon(bacon.fit1, n.iter = 100, plot_priors = F)


adam.fit1 <- adam(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  #bottom_depth = 800,
  K = c(10, 10),
  nu = 6,
  record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 0, chains = 3)

plot_stan_bacon(adam.fit1, n.iter = 1000, plot_priors = F)

a1 <- rstan::summary(adam.fit1$fit)
a1 <- as_tibble(a1$summary, rownames = "par")



adam.fit2b <- adam(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
 # K = c(3,3,3,3,3,3),
  K = c(9, 9, 9),
  nu = 6,
  record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 0, chains = 3)

plot_stan_bacon(adam.fit2b, n.iter = 1000, plot_priors = F)

a2 <- rstan::summary(adam.fit2$fit)
a2 <- as_tibble(a2$summary, rownames = "par")

plot_stan_bacon(adam.fit2, n.iter = 1000, plot_priors = F)
plot_stan_bacon(adam.fit2b, n.iter = 1000, plot_priors = F)

a2b <- rstan::summary(adam.fit2b$fit)
a2b <- as_tibble(a2b$summary, rownames = "par")



summary(adam.fit2$fit, par = c("alpha[1]", "shape"))$summary

traceplot(adam.fit2b$fit, par = c("alpha[1]", "shape"))

traceplot(adam.fit2b$fit, par = c("alpha[700]"))

#shinystan::launch_shinystan(adam.fit3$fit)


adam.fit3 <- adam(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  K = c(9,9,9),
  nu = 6,
  record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 1, chains = 3)

plot_stan_bacon(adam.fit3, n.iter = 1000, plot_priors = F)

a3 <- rstan::summary(adam.fit3$fit)
a3 <- as_tibble(a3$summary, rownames = "par")

summary(adam.fit3$fit, par = c("alpha[1]", "shape"))$summary
summary(adam.fit3$fit, par = c("w"))$summary

traceplot(adam.fit3$fit, par = c("alpha[1]", "shape"))

baconr:::plot_memory_prior_posterior(adam.fit3)



idx <- as_tibble(alpha_indices(c(9,9,9))[1:3]) %>%
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

alph %>%
  mutate(alpha_idx = as.numeric(alpha_idx)) %>%
  ggplot(aes(x = alpha_idx, y = mean, colour = lvl)) +
  geom_point()






