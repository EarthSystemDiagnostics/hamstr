## Fit Bayesian model -----

library(cmdstanr)
library(posterior)
library(hamstr)
library(tidyverse)



### Test on boxcore data ---------
bc <- (replicatecoredata::OR1_1218.14C) %>%
  arrange(depth_cm)


bc %>%
  ggplot(aes(x = depth_cm, y = age.14C.cal, colour = factor(subcore),
             shape = factor(n.forams))) +
  geom_point() +
  facet_wrap(~tube)

bc.sub <- bc %>%
  filter(tube == 1) %>%
  #filter(n.forams == 200) %>%
  arrange(depth_cm) %>%
  mutate(dpt = 1:n())#%>%

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
                    model_bioturbation = TRUE,
                    inflate_errors = FALSE,
                    L_prior_mean = 10,
                    L_prior_shape = 2,
                 cores = 3)


ham.bt.bc.fix <- hamstr(bc.sub$depth_cm, obs_age = bc.sub$age.14C.cal, obs_err = bc.sub$age.14C.cal.se,
                    n_ind = bc.sub$n.forams,
                    #n_ind = rep(200, nrow(bc.sub)),
                    model_bioturbation = TRUE,
                    inflate_errors = FALSE,
                    L_prior_mean = 10,
                    L_prior_shape = 0,
                    cores = 3)


ham.bt.bc.sp <- hamstr(bc.sub$depth_cm, obs_age = bc.sub$age.14C.cal, obs_err = bc.sub$age.14C.cal.se,
                    n_ind = bc.sub$n.forams,
                    #n_ind = rep(200, nrow(bc.sub)),
                    model_bioturbation = TRUE,
                    inflate_errors = FALSE,
                    L_prior_mean = 10,
                    L_prior_shape = 2000,
                    cores = 3)

#hist(rgamma(1e03, 2000, scale = 10/2000))


p1 <- plot(ham.bc, type = "age")
p2 <- plot(ham.bt.bc, type = "age")
p3 <- plot(ham.bt.bc.sp, type = "age")
p4 <- plot(ham.bt.bc.fix, type = "age")

egg::ggarrange(p1, p2, p3)

pairs(ham.bt.bc$fit, pars = c("L",
                              #"alpha[1]", "x[1]",
                              "age0", "R", "bt_age[1]"))

pairs(ham.bt.bc$fit, pars = c("R", "age0"))
pairs(ham.bt.bc$fit, pars = c("L", "bt_age[1]", "bt_age[2]"))

#shinystan::launch_shinystan(ham.bt.bc$fit)

post <- as.data.frame(ham.bt.bc.sp$fit, pars = c("L", "bt_age")) %>%
  as_tibble() %>%
  dplyr::mutate(iter = 1:nrow(.)) %>%
  pivot_longer(cols = -iter) %>%
  mutate(dpt = readr::parse_number(name))

prior <- tibble(
  x = seq(5, 15, length.out = 1000),
  y = dgamma(x, ham.bt.bc.sp$data$L_prior_shape,
             rate = ham.bt.bc.sp$data$L_prior_shape/ham.bt.bc.sp$data$L_prior_mean)
)

post %>%
  filter(name == "L[1]") %>%
  ggplot(aes(x = value, after_stat(density))) +
  geom_histogram(alpha = 0.5, aes(fill = "Posterior", colour = "Posterior")) +
  geom_line(data = prior, aes(x =x, y = y, colour = "Prior", fill = "Prior")) +
  theme_bw()

rstan::summary(ham.bt.bc.sp$fit, pars = "L")$summary


tmp <- bc.sub %>%
  left_join(., filter(post, name != "L[1]"))


p2 +
  geom_violin(data = tmp,
              aes(x = depth_cm, y = value, group = as.factor(dpt)),
              fill = "pink", colour = "pink", scale = "area",
              position = position_identity(), alpha = 0.15)



# ## compile model ------
# file <- file.path("inst/stan/hamstr.stan")
# mod <- cmdstan_model(file)
#
#
# stan_dat <- make_stan_dat_hamstr(
#   depth = bc.sub$depth_cm,
#   obs_age = bc.sub$age.14C.cal,
#   obs_err = bc.sub$age.14C.cal.se,
#   n_ind = bc.sub$n.forams,
#   min_age = 0,
#   mem_mean = 0.5,
#   mem_strength = 10,
#   scale_R = 1,
#   inflate_errors = 1,
#   pad_top_bottom = 0,
#   model_bioturbation = TRUE,
#   L_prior_mean = 10,
#   L_prior_shape = 20
# )
#
#
# stan_dat <- stan_dat[-26]
#
#
# stan.bt.bc <- mod$sample(data = stan_dat, chains = 3, parallel_chains = 3)
#
# stan_posterior <- stan.bt.bc$draws(variables = c("L", "alpha[1]", "age_het", "bt_age")) %>% as_draws_df()
#
#
# sl <- stan_posterior %>%
#   select(starts_with("L")) %>%
#   pivot_longer(everything()) %>%
#   arrange(name)
#
# hist(sl$value, freq = FALSE, xlim = c(0, 50), breaks = 50)
# lines(seq(0, 150), dgamma(0:150, stan_dat$L_prior_shape, stan_dat$L_prior_shape/stan_dat$L_prior_mean), col = "red")
#
# sl <- stan_posterior %>%
#   select(starts_with("bt_age")) %>%
#   pivot_longer(everything()) %>%
#   arrange(name)
#
#
# age.dat <- filter(mutate(bc, dpt = 1:nrow(bc)),
#                   dpt < 5)
#
# sl %>%
#   mutate(dpt = readr::parse_number(name)) %>%
#   filter(dpt < 5) %>%
#   ggplot(aes(x = value)) +
#   geom_histogram(bins = 50) +
#   facet_wrap(~dpt) +
#   geom_vline(data = age.dat, aes(xintercept = age.14C.cal), colour = "red") +
#   geom_vline(data = age.dat, aes(xintercept = age.14C.cal + age.14C.cal.se), colour = "blue")+
#   geom_vline(data = age.dat, aes(xintercept = age.14C.cal - age.14C.cal.se), colour = "blue")
#
#
#















