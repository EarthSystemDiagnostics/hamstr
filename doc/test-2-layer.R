library(rstan)
library(baconr)
library(tidyverse)


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

# make_stan_dat_K1 <- function(depth, obs_age, obs_err,
#                              K1 = 10, K = 100, nu = 6,
#                              record_prior_acc_mean_mean = 20,
#                              record_prior_acc_mean_shape = 1.5,
#                              record_prior_acc_shape_mean = 1.5,
#                              record_prior_acc_shape_shape = 1.5,
#                              mem_mean = 0.7, mem_strength = 4,
#                              inflate_errors = 0) {
#   
#   # Pretty sure nu = 6 is equivalent to the default parameterisation of t.a = 3, t.b = 4 in Bacon
#   # when t.b = t.a +1, there are 2a degrees of freedom, and the distribution is symetric.
#   
#   l <- c(as.list(environment()))
#   
#   # Transformed arguments
#   l$N <- length(depth)
#   
#   stopifnot(K%%K1 == 0)
#   l$whichK1 = rep(1:K1, each = K / K1)
#   
#   l$mem_alpha = mem_strength * mem_mean
#   l$mem_beta = mem_strength * (1-mem_mean)
#   
#   l$mem_mean = mem_mean
#   l$mem_strength = mem_strength
#   
#   # Set start depth to 5% less than first depth observation, but do not allow negative depths
#   depth_range <- diff(range(depth))
#   buff_5 <- 0.05 * depth_range
#   strt_dpth <- depth[1] - buff_5
#   strt_dpth[strt_dpth < 0] <- 0
#   end_dpth <- tail(depth, 1) + buff_5
#   depth_range = end_dpth - strt_dpth
#   
#   
#   l$c = 1:K
#   l$delta_c = depth_range / K
#   l$c_depth_bottom = l$delta_c * l$c + strt_dpth
#   l$c_depth_top = c(strt_dpth, l$c_depth_bottom[1:(K-1)])
#   
#   # Index for which sections the target depth is in
#   l$which_c = sapply(l$depth, function(d) which.max((l$c_depth_bottom < d) * (l$c_depth_bottom - d) ))
#   
#   
#   return(l)
# }
# 
# GetInits <- function(dat){
#   l <- list(
#     R = runif(1, 0.1, 0.9),
#     alpha = (abs(rnorm(dat$K, dat$record_prior_acc_mean_mean, dat$record_prior_acc_mean_mean/3))),
#     record_acc_mean = (abs(rnorm(1, dat$record_prior_acc_mean_mean, dat$record_prior_acc_mean_mean/3))),
#     section_acc_mean = as.array(abs(rnorm(dat$K1, dat$record_prior_acc_mean_mean, dat$record_prior_acc_mean_mean/3))),
#     
#     record_acc_shape = rnorm(1, 1.5, 1.5/3),
#     
#     age0 = runif(1, 0, dat$obs_age)
#   )
#   
#   # need to make this conditional and make sure initial values are arrays! 
#   if (dat$inflate_errors == 1){
#     l$infl_mean = as.array(abs(rnorm(1, 0, 0.1)))
#     l$infl_sd = as.array(abs(rnorm(1, 0, 0.1)))
#     l$infl = abs(rnorm(dat$N, 0, 0.1))   
#   } else {
#     l$infl_mean = numeric(0)
#     l$infl_sd = numeric(0)
#     l$infl = numeric(0)   
#   }
#   
#   return(l)
#   
# }


name <- "BAMBILI2"
#dat <- read.csv(paste0("/Users/andrewdolman/Dropbox/Work/AWI/Data/terrestrial-age-models/terr_14C_min10_dates-2020.03.04_15-19-42/", name, "/", name,".csv")) %>%
dat <- read.csv(paste0("../envi-age-modelling/working-data/terr_14C_min10_dates-2020.03.04_15-19-42/", name, "/", name,".csv")) %>%
  filter(depth > 0)

dat <- ecustools::CalibrateAge(dat, age.14C = "age", age.14C.se = "error") %>%
  filter(complete.cases(age.14C.cal))

dat %>%
  ggplot(aes(x = depth, y = age.14C.cal)) +
  geom_point() +
  geom_line()

acc.mean <- coef(MASS::rlm(age~depth, data = dat))[2]


acc.mean


#options(mc.cores = parallel::detectCores())
#options(mc.cores = 1)

# dat.fit.2l <- make_stan_dat_K1(
#   depth = dat$depth,
#   obs_age = dat$age.14C.cal,
#   obs_err = dat$age.14C.cal.se,
#   K = 100,
#   K1 = 10,
#   nu = 6,
#   record_prior_acc_mean_mean = acc.mean,
#   record_prior_acc_mean_shape = 1.5,
#   record_prior_acc_shape_mean = 1.5,
#   record_prior_acc_shape_shape = 1.5,
#   mem_mean = 0.7, mem_strength = 4,
#   inflate_errors = 1)

sarn.fit1 <- sarnie(
  depth = dat$depth,
  obs_age = dat$age.14C.cal,
  obs_err = dat$age.14C.cal.se,
  K = 100,
  K1 = 10,
  nu = 6,
  record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 1, chains = 3)



#two_l_mod <- rstan::stan_model("inst/stan/bacon_2_layer_hier_acc_mean.stan")

twoFit.fit <- rstan::sampling(two_l_mod, data = dat.fit.2l, chains = 3, 
                              init = inits, iter = 4000,
                              control = list(adapt_delta = 0.8))

twoFit <- list(fit = twoFit.fit, data = dat.fit.2l)

s1 <- summary(twoFit$fit)
s1 <- as_tibble(s1$summary, rownames = "par")

samples <- as.data.frame(twoFit.fit) %>% as_tibble()

table(is.na(s1$Rhat))
range(s1$Rhat, na.rm = T)
range(s1$n_eff)

options(width = 85)
print(twoFit$fit, par = c("record_acc_mean"))

print(twoFit$fit, par = c("section_acc_mean"))

print(twoFit$fit, par = c("infl_shape", "infl_mean"))
print(twoFit$fit, par = c("infl"))

print(twoFit$fit, par = c("w", "R"))

set.seed(20170406)
p.multilevel <- plot_stan_bacon(twoFit, 1000, plot_priors = F)
#p.multilevel

# infl errors
infl <- s1 %>%
  filter(grepl("infl[", par, fixed = T),
         grepl("err_infl", par, fixed = T) == FALSE)

dat.infl <- tibble(depth = twoFit$data$depth, obs_age = twoFit$data$obs_age,
                   obs_err = twoFit$data$obs_err,
                   infl = infl$mean) %>%
  mutate(infl.err = obs_err + obs_err * infl)

insertLayer(p.multilevel, after = 1,
            geom_linerange(data = dat.infl,
                           aes(x = depth,
                               y = obs_age,
                               ymax = obs_age + 2*infl.err,
                               ymin = obs_age - 2*infl.err),
                           colour = "green", group = NA)) +
  labs(subtitle = paste0("Inflation factor = ",
                         paste0(round(infl$mean, 1), collapse = ", ")))



#ggsave(paste0("doc/", name, "multilevel.png"), p.multilevel)

rstan::traceplot(twoFit$fit, par = "infl", inc_warmup = T)

rstan::traceplot(twoFit$fit, par = "section_acc_mean", inc_warmup = F)

rstan::traceplot(twoFit$fit, par = "record_acc_shape", inc_warmup = T)

d2 <- s1 %>%
  filter(grepl("alpha", par)) %>%
  filter(grepl("beta", par) == FALSE) %>%
  mutate(depth = 0.5+(1:100)/10)

s1 %>%
  filter(grepl("section_acc_mean", par)) %>%
  filter(grepl("beta", par) == FALSE) %>%
  mutate(depth = 1:10) %>%
  ggplot(aes(x = depth, y = mean)) +
  geom_point() +
  geom_line() +
  geom_line(data = d2, colour = "Red")


# mean infl is not neccesarily the same as mean(shape) / mean(beta)
samples %>%
  select(starts_with("infl_")) %>%
  mutate(mean = infl_shape / infl_beta) %>%
  summarise_if(is.numeric, mean) %>%
  mutate(mean2 = infl_shape / infl_beta)



#######

dat.fit <- stan_bacon_hier(
  depth = dat$depth,
  obs_age = dat$age.14C.cal,
  obs_err = dat$age.14C.cal.se,
  K = 100,
  nu = 6,
  acc_mean_mean = 20, acc_mean_shape = 1.5, acc_alpha = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  iter = 2000, chains = 3)

set.seed(20170406)
p.single.level <- plot_stan_bacon(dat.fit, 1000)
p.single.level
ggsave(paste0("doc/", name, "singlelevel.png"), p.single.level)


##

library(rbacon)

thick <- diff(range(dat$depth)) / 100

b1 <- Bacon(core = "BLUELAKE", coredir = "../envi-age-modelling/working-data/terr_14C_min10_dates-2020.03.04_15-19-42",
      thick = thick)





a <- 0.05
b <- 2.25
x <- seq(-1, 20, length.out = 1000)
d <- dgamma(x, a, scale = b)
d2 <- dcauchy(x, scale = 2.5)
plot(x, d, type = "l")
lines(x ,d2)
abline(v = 0)
abline(v = b/a)




