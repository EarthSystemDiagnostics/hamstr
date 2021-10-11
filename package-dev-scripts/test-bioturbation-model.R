

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
  L_prior_mean = 20,
  L_prior_shape = 2,
  extra_var = 1
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

s <- rgamma(1e04, 2, rate = 2/10)
hist(s, 51)
abline(v = mean(s), col = "red")
#abline(v = median(s), col = "red")
abline(v = 10, col = "red")

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

s <- rgamma(1e03, shape = 1, rate = 1/200)
hist(s-100, 50)
abline(v = mean(s))
abline(v = 100, col = "Red")

# file2 <- file.path("inst/stan/hamstr.stan")
# mod2 <- cmdstan_model(file2)
# mod2$sample(data = stan_dat, chains = 3)

