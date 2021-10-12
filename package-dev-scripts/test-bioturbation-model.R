

## Fit Bayesian model

library(cmdstanr)
library(posterior)
library(hamstr)
library(tidyverse)



MSB2K_cal <- calibrate_14C_age(MSB2K, age.14C = "age", age.14C.se = "error")

stan_dat <- make_stan_dat_hamstr(
  depth = MSB2K_cal$depth,
  obs_age = MSB2K_cal$age.14C.cal,
  obs_err = MSB2K_cal$age.14C.cal.se,
  min_age = 0,
  mem_mean = 0.5,
  mem_strength = 10,
  scale_R  =1,
  inflate_errors = 0,
  pad_top_bottom = 0

)

stan_dat <- append(stan_dat, list(
  n_ind = rep(1, nrow(MSB2K_cal)),
  L_prior_mean = 10,
  L_prior_shape = 2,
  model_bioturbation = 1
))


stan_dat <- stan_dat[-22]

str(stan_dat)


# compile model ------
file <- file.path("inst/stan/hamstr-bioturbation.stan")
mod <- cmdstan_model(file)

tmp1 <- mod$sample(data = stan_dat, chains = 3, parallel_chains = 3,#, iter_warmup = 100, iter_sampling = 100,
                   init = 1)

stan_posterior <- tmp1$draws(variables = c("L", "alpha[1]", "age_het", "bt_age")) %>% as_draws_df()

sl <- stan_posterior %>%
  select(starts_with("L")) %>%
  pivot_longer(everything()) %>%
  arrange(name)

hist(sl$value, freq = FALSE, xlim = c(0, 50))
lines(seq(0, 50), dgamma(0:50, stan_dat$L_prior_shape, stan_dat$L_prior_shape/stan_dat$L_prior_mean), col = "red")

sl <- stan_posterior %>%
  select(starts_with("bt_age")) %>%
  pivot_longer(everything()) %>%
  arrange(name)


age.dat <- filter(mutate(MSB2K_cal, dpt = 1:nrow(MSB2K_cal)),
                  dpt < 5)

sl %>%
  mutate(dpt = readr::parse_number(name)) %>%
  filter(dpt < 5) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 50) +
  facet_wrap(~dpt) +
  geom_vline(data = age.dat, aes(xintercept = age.14C.cal), colour = "red") +
  geom_vline(data = age.dat, aes(xintercept = age.14C.cal + age.14C.cal.se), colour = "blue")+
  geom_vline(data = age.dat, aes(xintercept = age.14C.cal - age.14C.cal.se), colour = "blue")


# file2 <- file.path("inst/stan/hamstr.stan")
# mod2 <- cmdstan_model(file2)
# mod2$sample(data = stan_dat, chains = 3)

###

ham.bt.1 <- hamstr(
  MSB2K_cal$depth, obs_age = MSB2K_cal$age.14C.cal, obs_err = MSB2K_cal$age.14C.cal.se,
  n_ind = rep(100, nrow(MSB2K_cal)),

  model_bioturbation = 1,
  chains = 1
)

plot(ham.bt.1)

rstan::summary(ham.bt.1$fit, pars = "L")



### Test on multicore data -------

library(replicatecoredata)

muc <- (replicatecoredata::MUC.14C) %>%
  arrange(depth_cm)


ham.muc <- hamstr(muc$depth_cm, obs_age = muc$age.14C.cal, obs_err = muc$age.14C.cal.se,
                  cores = 3)
plot(ham.muc)

ham.muc.err <- hamstr(muc$depth_cm, obs_age = muc$age.14C.cal, obs_err = muc$age.14C.cal.se,
                      inflate_errors = TRUE,
                      cores = 3)
plot(ham.muc.err)


ham.muc.bt <- hamstr(muc$depth_cm, obs_age = muc$age.14C.cal, obs_err = muc$age.14C.cal.se,
                     n_ind = muc$n.forams,
                      model_bioturbation = 1,
                    L_prior_shape = 2,
                     #K = c(3,3,3,3),
                     cores = 3)

plot(ham.muc.bt)

rstan::summary(ham.muc.bt$fit, pars = "L")


stan_dat <- make_stan_dat_hamstr(
  depth = muc$depth_cm,
  obs_age = muc$age.14C.cal,
  obs_err = muc$age.14C.cal.se,
  min_age = 0,
  mem_mean = 0.5,
  mem_strength = 10,
  scale_R  =1,
  inflate_errors = 0,
  pad_top_bottom = 0

)

stan_dat <- append(stan_dat, list(
  n_ind = muc$n.forams,
  L_prior_mean = 10,
  L_prior_shape = 2
))


stan_dat <- stan_dat[-22]


stan.bt.muc <- mod$sample(data = stan_dat, chains = 3, parallel_chains = 3,#, iter_warmup = 100, iter_sampling = 100,
                   init = 1)

stan_posterior <- stan.bt.muc$draws(variables = c("L", "alpha[1]", "age_het", "bt_age")) %>% as_draws_df()


sl <- stan_posterior %>%
  select(starts_with("L")) %>%
  pivot_longer(everything()) %>%
  arrange(name)

hist(sl$value, freq = FALSE, xlim = c(0, 150), breaks = 50)
lines(seq(0, 150), dgamma(0:150, stan_dat$L_prior_shape, stan_dat$L_prior_shape/stan_dat$L_prior_mean), col = "red")

sl <- stan_posterior %>%
  select(starts_with("bt_age")) %>%
  pivot_longer(everything()) %>%
  arrange(name)


age.dat <- filter(mutate(muc, dpt = 1:nrow(muc)),
                  dpt < 5)

sl %>%
  mutate(dpt = readr::parse_number(name)) %>%
  filter(dpt < 5) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 50) +
  facet_wrap(~dpt) +
  geom_vline(data = age.dat, aes(xintercept = age.14C.cal), colour = "red") +
  geom_vline(data = age.dat, aes(xintercept = age.14C.cal + age.14C.cal.se), colour = "blue")+
  geom_vline(data = age.dat, aes(xintercept = age.14C.cal - age.14C.cal.se), colour = "blue")


ham.bt.muc <- list(dat = stan_dat, fit = rstan::read_stan_csv(stan.bt.muc$output_files()))
#class(ham.bt.muc) <- append(class(ham.bt.muc), c("hamstr_fit"))
class(ham.bt.muc) <- c("hamstr_fit", "list")

plot(ham.bt.muc)


### Test on boxcore data ---------
bc <- (replicatecoredata::OR1_1218.14C) %>%
  arrange(depth_cm)


bc %>%
  ggplot(aes(x = depth_cm, y = age.14C.cal, colour = factor(subcore),
             shape = factor(n.forams))) +
  geom_point() +
  facet_wrap(~tube)

bc.sub <- bc #%>%
 #filter(tube == 8) %>%
 #filter(n.forams == 200)


ham.bc <- hamstr(bc.sub$depth_cm, obs_age = bc.sub$age.14C.cal, obs_err = bc.sub$age.14C.cal.se,
                 cores = 3)

plot(ham.bc)

ham.bc.err <- hamstr(bc.sub$depth_cm, obs_age = bc.sub$age.14C.cal, obs_err = bc.sub$age.14C.cal.se,
                     cores = 3,
                     inflate_errors = TRUE)

plot(ham.bc.err)


ham.bt.bc <- hamstr(bc.sub$depth_cm, obs_age = bc.sub$age.14C.cal, obs_err = bc.sub$age.14C.cal.se,
                    n_ind = bc.sub$n.forams,
                    #n_ind = rep(1, nrow(bc.sub)),
                    model_bioturbation = TRUE,
                    #L_prior_shape = 2000,
                 cores = 3)

plot(ham.bt.bc)

plot(ham.bt.bc, type = "age_") +
  geom_point(data = bc, aes(x = depth_cm, y = age.14C.cal))

plot(ham.bc, type = "age_") +
  geom_point(data = bc, aes(x = depth_cm, y = age.14C.cal))



rstan::summary(ham.bt.bc$fit, pars = "L")



tmp <- as.data.frame(ham.bt.bc$fit, pars = c("bt_age")) %>%
  as_tibble() %>%
  dplyr::mutate(iter = 1:nrow(.)) %>%
  pivot_longer(cols = -iter) %>%
  mutate(dpt = readr::parse_number(name))


tmp <- bc.sub %>%
  arrange(depth_cm) %>%
  mutate(dpt = 1:n()) %>%
  left_join(., tmp)


tmp %>%
  ggplot(aes(x = depth_cm, group = as.factor(dpt), y = value, fill = factor(n.forams))) +
  geom_violin() +
  facet_wrap(~depth_cm, scales = "free")



hist(tmp$`bt_age[1]`)

stan_dat <- make_stan_dat_hamstr(
  depth = bc$depth_cm,
  obs_age = bc$age.14C.cal,
  obs_err = bc$age.14C.cal.se,
  min_age = 0,
  mem_mean = 0.5,
  mem_strength = 10,
  scale_R  =1,
  inflate_errors = 0,
  pad_top_bottom = 0,
  model_bioturbation = FALSE

)

stan_dat <- append(stan_dat, list(
  n_ind = bc$n.forams,
  L_prior_mean = 10,
  L_prior_shape = 2
))


stan_dat <- stan_dat[-22]


stan.bt.bc <- mod$sample(data = stan_dat, chains = 3, parallel_chains = 3,#, iter_warmup = 100, iter_sampling = 100,
                          init = 1)

stan_posterior <- stan.bt.bc$draws(variables = c("L", "alpha[1]", "age_het", "bt_age")) %>% as_draws_df()


sl <- stan_posterior %>%
  select(starts_with("L")) %>%
  pivot_longer(everything()) %>%
  arrange(name)

hist(sl$value, freq = FALSE, xlim = c(0, 50), breaks = 50)
lines(seq(0, 150), dgamma(0:150, stan_dat$L_prior_shape, stan_dat$L_prior_shape/stan_dat$L_prior_mean), col = "red")

sl <- stan_posterior %>%
  select(starts_with("bt_age")) %>%
  pivot_longer(everything()) %>%
  arrange(name)


age.dat <- filter(mutate(bc, dpt = 1:nrow(bc)),
                  dpt < 5)

sl %>%
  mutate(dpt = readr::parse_number(name)) %>%
  filter(dpt < 5) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 50) +
  facet_wrap(~dpt) +
  geom_vline(data = age.dat, aes(xintercept = age.14C.cal), colour = "red") +
  geom_vline(data = age.dat, aes(xintercept = age.14C.cal + age.14C.cal.se), colour = "blue")+
  geom_vline(data = age.dat, aes(xintercept = age.14C.cal - age.14C.cal.se), colour = "blue")


















