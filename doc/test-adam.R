# test adam
library(tidyverse)
library(baconr)
library(rstan)

# build -----

# devtools::document()
# pkgbuild::compile_dll(force = TRUE)
# devtools::load_all()

# devtools::install(quick = FALSE)
# devtools::install(quick = FALSE, dependencies = FALSE)

#sm <- stan_model("inst/stan/adam.stan")
#stan_data <- baconr::make_stan_dat_adam(depth = dat1$depth, obs_age = dat1$age.14C.cal, obs_err = dat1$age.14C.cal.se)

#rstan::sampling(sm, data = stan_data, iter = 10)


#Load and filter data
all.terr.dat <- readr::read_csv2("doc/Dating_Data.csv") %>%
  filter(complete.cases(age, depth, e.older)) %>%
  mutate(sigma.age = e.older)

all.terr.14C.dat <- all.terr.dat %>%
  filter(age.type == "Radiocarbon years BP") %>%
  group_by(DataName) %>%
  mutate(n.dates = n()) %>%
  ungroup()


#name <- "BEEFPAST"
name <- "BRUCHBG1"

dat2 <- all.terr.14C.dat %>%
  filter(DataName == name)

dat1 <- calibrate_14C_age(dat2, age.14C = "age", age.14C.se = "sigma.age") %>%
  filter(complete.cases(age.14C.cal)) %>%
  select(depth, age.14C.cal, age.14C.cal.se, age, sigma.age) %>%
  mutate(age.14C.cal = round(age.14C.cal),
         age.14C.cal.se = round(age.14C.cal.se))

dat1 %>%
  ggplot(aes(x = depth, y = age.14C.cal)) +
  geom_point() +
  geom_line()


#options(mc.cores = parallel::detectCores())

options(mc.cores = 4)

adam.fit1 <- adam(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  K = c(100),
  nu = 6,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 0, chains = 3)


ad_K100_rib <- plot_adam(adam.fit1, plot_diagnostics = T)
ad_K100_spag <- plot_adam(adam.fit1, type = "spaghetti", n.iter = 100, plot_diagnostics = T)

ggsave("doc/ad_K100_rib.png", ad_K100_rib, width = 8, height = 6)
ggsave("doc/ad_K100_spag.png", ad_K100_spag, width = 8, height = 6)

traceplot(adam.fit1$fit, pars = c("shape", "alpha[1]"), inc_warmup = T)



adam.fit2 <- adam(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  K = baconr:::optimal_K(100, 10),
  nu = 6,
  shape = 1.5,
  mem_mean = 0.5, mem_strength = 2,
  scale_R = FALSE,
  inflate_errors = FALSE,
  chains = 3)

# warmup sample
# chain:1 53.281 82.095
# chain:2 53.627 73.500
# chain:3 55.919 73.026

# warmup sample
# chain:1 37.423 40.540
# chain:2 34.730 37.017
# chain:3 38.712 29.908

get_elapsed_time(adam.fit2$fit)

ad_K10_100_rib <- plot_adam(adam.fit2, plot_diagnostics = T) 

ad_K10_100_spag <- plot_adam(adam.fit2, type = "spaghetti", n.iter = 1000, plot_diagnostics = T)

ad_K10_100_spag_nod <- plot_adam(adam.fit2, type = "spaghetti", n.iter = 100, plot_diagnostics = F)

plot_hierarchical_acc_rate(adam.fit2)

plot_hierarchical_acc_rate(adam.fit2) +
  facet_grid(lvl~.)

plot_acc_mean_prior_posterior(adam.fit2)

ggsave("doc/ad_K10_100_rib.png", ad_K10_100_rib, width = 8, height = 6)
ggsave("doc/ad_K10_100_spag.png", ad_K10_100_spag, width = 8, height = 8)

ggsave("doc/ad_K10_100_spag_nod.png", ad_K10_100_spag_nod, width = 8, height = 6)


a2 <- rstan::summary(adam.fit2$fit)
a2 <- as_tibble(a2$summary, rownames = "par")


summary(adam.fit2$fit, par = c("alpha[1]", "beta"))$summary

traceplot(adam.fit2$fit, pars = c("shape", "alpha[1]"), inc_warmup = T)

traceplot(adam.fit2$fit, pars = paste0("alpha[", 91:111, "]"), inc_warmup = T)

traceplot(adam.fit2$fit, pars = c("shape", "alpha[1]"), inc_warmup = T)

traceplot(adam.fit2$fit, pars = c("shape"), inc_warmup = T)

traceplot(adam.fit2$fit, pars = c("w", "R"), inc_warmup = T)

stan_dat <- make_stan_dat_adam(depth = dat1$depth, obs_age = dat1$age.14C.cal, obs_err = dat1$age.14C.cal.se,
                               K=c(10, 10),
                               inflate_errors = 1)


p_idx <- stan_dat$parent
p_idx[p_idx == 0] <- NA

alph <- a2 %>%
  filter(grepl("alpha", par)) %>%
  rename(alpha = mean) %>%
  select(alpha) %>%
  mutate(parent_alpha = alpha[p_idx])

bet <- a2 %>%
  filter(grepl("beta", par))%>%
  rename(beta = mean) %>%
  select(beta)

shp <- a2 %>%
  filter(grepl("shape_vec", par))%>%
  rename(shape = mean) %>%
  select(shape)

prs <- bind_cols(alph, bet) %>%
  mutate(shape = c(NA, shp$shape)) %>%
  mutate(beta_calc = shape / parent_alpha)


# extract a few iterations

as_tibble(as.data.frame(adam.fit2$fit, pars = c("alpha")))


prs <- tibble(
shape = c(NA, rstan::extract(adam.fit2$fit, "shape_vec")[[1]][1,]),
alpha = rstan::extract(adam.fit2$fit, "alpha")[[1]][1,],
beta = rstan::extract(adam.fit2$fit, "beta")[[1]][1,]
) %>%
  mutate(parent_alpha = alpha[p_idx])%>%
  mutate(beta_calc = shape / parent_alpha)





## with 3 layers

adam.fit2b <- adam(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  K = baconr:::optimal_K(1000, 10),
  nu = 6,
  #acc_mean_prior = acc.mean,
  #record_prior_acc_shape_mean = 1.5,
  #record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 0, chains = 3)


ad_K10_100_1000_rib <- plot_adam(adam.fit2b, plot_diagnostics = T)
ad_K10_100_1000_spag <- plot_adam(adam.fit2b, type = "spaghetti", n.iter = 100, plot_diagnostics = T)

ad_K10_100_1000_spag_nod <- plot_adam(adam.fit2b, type = "spaghetti", n.iter = 100, plot_diagnostics = F)

ad_K10_100_1000_HA <- plot_hierarchical_acc_rate(adam.fit2b)


ggsave("doc/ad_K10_100_1000_rib.png", ad_K10_100_1000_rib, width = 8, height = 6)
ggsave("doc/ad_K10_100_1000_spag.png", ad_K10_100_1000_spag, width = 8, height = 6)

ggsave("doc/ad_K10_100_1000_spag_nod.png", ad_K10_100_1000_spag_nod, width = 8, height = 6)


ggsave("doc/ad_K10_100_1000_HA.png", ad_K10_100_1000_HA, width = 6, height = 4)






#shinystan::launch_shinystan(adam.fit3$fit)
## adam.fit3 -----

adam.fit3 <- adam(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  K = baconr:::optimal_K(100, 10),
  shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  infl_shape_shape = 1,
  infl_shape_mean = 1,
  #infl_mean_shape = 1,
  #infl_mean_mean = 3,
  inflate_errors = TRUE, chains = 3)

plot_adam(adam.fit3, type = "ribbon", plot_diagnostics = TRUE)
plot_adam(adam.fit3, type = "spaghetti", n.iter = 1000, plot_diagnostics = TRUE)

print(adam.fit3$fit, pars = c("infl", "infl_mean",  "infl_shape"#, "infl_sigma"
                              ))
traceplot(adam.fit3$fit, pars = c("infl_mean", "infl_shape"#, "infl_sigma"
                                  ), inc_warmup = T)


plot_infl_prior_posterior(adam.fit3)

pairs(adam.fit3$fit, pars = c("infl_mean", "infl_shape"))

a3 <- rstan::summary(adam.fit3$fit)
a3 <- as_tibble(a3$summary, rownames = "par")

summary(adam.fit3$fit, par = c("alpha[1]", "shape"))$summary
summary(adam.fit3$fit, par = c("w"))$summary

traceplot(adam.fit3$fit, par = c("alpha[1]", "shape"))

traceplot(adam.fit3$fit, par = c("infl_mean", "infl_shape"))
traceplot(adam.fit3$fit, par = c("infl"))

summary(adam.fit3$fit, par = c("infl_mean", "infl_shape"))$summary


as.data.frame(adam.fit3$fit, pars = c("infl[1]", "obs_err_infl[1]")) %>% 
  as_tibble() %>% 
  mutate(obs_err = dat1$age.14C.cal.se[1],
         tst = obs_err + `infl[1]`) %>% 
  ggplot(aes(x = `obs_err_infl[1]`, y = `infl[1]`)) +
  geom_point()



###

name <- "BIBER"

dat3 <- all.terr.14C.dat %>%
  filter(DataName == name)

dat3 <- ecustools::CalibrateAge(dat3, age.14C = "age", age.14C.se = "sigma.age") %>%
  filter(complete.cases(age.14C.cal)) %>%
  select(depth, age.14C.cal, age.14C.cal.se, age, sigma.age) %>%
  mutate(age.14C.cal = round(age.14C.cal),
         age.14C.cal.se = round(age.14C.cal.se))

dat3 %>%
  ggplot(aes(x = depth, y = age.14C.cal)) +
  geom_point() +
  geom_line()


adam.fit4 <- adam(
  depth = dat3$depth,
  obs_age = dat3$age.14C.cal,
  obs_err = dat3$age.14C.cal.se,
  #top_depth = 1,
  K = 100,
  nu = 6,
  #acc_mean_prior = acc.mean,

  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 0, chains = 3)

plot_adam(adam.fit4, type = "ribbon", plot_diagnostics = TRUE)
plot_adam(adam.fit4, type = "spaghetti", n.iter = 1000, plot_diagnostics = TRUE)



adam.fit4a <- adam(
  depth = dat3$depth,
  obs_age = dat3$age.14C.cal,
  obs_err = dat3$age.14C.cal.se,
  #top_depth = 1,
  K = 100,
  nu = 6,
  #acc_mean_prior = acc.mean,

  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 1, chains = 3)

plot_adam(adam.fit4a, type = "ribbon", plot_diagnostics = TRUE)
plot_adam(adam.fit4a, type = "spaghetti", n.iter = 1000, plot_diagnostics = TRUE)



adam.fit4b <- adam(
  depth = dat3$depth,
  obs_age = dat3$age.14C.cal,
  obs_err = dat3$age.14C.cal.se,
  #top_depth = 1,
  K = baconr:::optimal_K(100, 10),
  nu = 6,
  #acc_mean_prior = acc.mean,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 1, chains = 3)

plot_adam(adam.fit4b, type = "ribbon", plot_diagnostics = TRUE)
plot_adam(adam.fit4b, type = "spaghetti", n.iter = 1000, plot_diagnostics = TRUE)

##


#

gamma_ex <- crossing(x = seq(-0.1, 4, length.out = 1000), shape = c(0.75, 1.5, 15), mu = 1) %>%
  mutate(scale = mu / shape) %>%
  group_by(shape) %>%
  do({
    tibble(
      x = .$x,
      d = dgamma(x, shape = .$shape[1], scale = .$scale[1])
    )
  })

fig_gamma_infl <- gamma_ex %>%
  ggplot(aes(x = x, y = d, colour = factor(shape), group = shape)) +
  geom_line() +
  theme_bw() +
  labs(colour = "shape", x = "Error inflation")

ggsave("doc/fig_gamma_infl.png", fig_gamma_infl, width = 6, height = 4)


#




idx <- as_tibble(baconr:::alpha_indices(c(3,3,3,3,3))[1:3]) %>%
  mutate(alpha_idx = as.character(alpha_idx))


alph <- a2 %>%
  #as_tibble(., rownames = "par") %>%
  filter(grepl("alpha", par, fixed = TRUE)) %>%
  separate(par, into = c("par", "alpha_idx")) %>%
  left_join(idx, .) %>%
  mutate(lvl = factor(lvl))


alph %>%
  ggplot(aes(x = mean)) +
  geom_histogram(aes()) +
  facet_wrap(~lvl, scales = "free_y", labeller = label_both)

alph %>%
  filter(lvl == 5) %>%
  ggplot() +
  geom_histogram(aes(x = mean, after_stat(density), fill = as.factor(parent))) +
  facet_wrap(~parent) +
  theme(legend.position = "none")


dat <- tibble(
  x = seq(0, 600, length.out = 1000),
  d = dgamma(x, shape = 1.5, rate = 1.5 / adam.fit2$data$acc_mean_prior)
)

alph %>%
  filter(lvl == tail(alph$lvl, 1)) %>%
  ggplot() +
  geom_histogram(aes(x = mean, after_stat(density))) +
  geom_line(data = dat, aes(x = x, y = d, colour = "Prior")) +
  facet_wrap(~lvl, scales = "free_y")


alph %>%
  ggplot(aes(x = Rhat)) +
  geom_histogram(aes(fill = lvl)) +
  expand_limits(x = c(0.95, 1.1)) +
  facet_wrap(~lvl, scales = "free_y")



