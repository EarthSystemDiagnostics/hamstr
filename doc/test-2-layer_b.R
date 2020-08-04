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


name <- "AREO"
#dat <- read.csv(paste0("/Users/andrewdolman/Dropbox/Work/AWI/Data/terrestrial-age-models/terr_14C_min10_dates-2020.03.04_15-19-42/", name, "/", name,".csv")) %>%
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
  #dput()



# dat1 <- structure(list(depth = c(52.5, 62.5, 67.5, 73.5, 74.5, 84.5,
#                                  96.5, 106.5, 107.5, 113.5, 117.5, 157.5), age.14C.cal = c(6599,
#                                                                                            17329, 18814, 21828, 21568, 34747, 36012, 39192, 35321, 41948,
#                                                                                            42868, 45432), age.14C.cal.se = c(82, 86, 50, 81, 380, 190, 736,
#                                                                                                                              429, 505, 1068, 679, 1121)), class = "data.frame", row.names = c(NA,
#                                                                                                                                                                                               -12L))


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
  acc_mean = acc.mean,
  mem_mean = 0.7, mem_strength = 4,
  chains = 3)

plot_stan_bacon(bacon.fit1, n.iter = 100, plot_priors = F)


sarn.fit1 <- sarnie(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  #bottom_depth = 800,
  K = 100,
  K1 = 10,
  nu = 6,
  record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 0, chains = 3)

plot_stan_bacon(sarn.fit1, n.iter = 100, plot_priors = T)

s1 <- rstan::summary(sarn.fit1$fit)
s1 <- as_tibble(s1$summary, rownames = "par")

plot(sarn.fit1$fit, par = "section_acc_mean")
plot(sarn.fit1$fit, par = "section_acc_mean")


sarn.fit2 <- sarnie(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  K = 1000,
  K1 = 10,
  nu = 6,
  record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 0, chains = 3)

plot_stan_bacon(sarn.fit2, n.iter = 100, plot_priors = T)



summary(sarn.fit2$fit, par = c("record_acc_mean", "record_acc_shape"))$summary
traceplot(sarn.fit2$fit, par = c("record_acc_mean", "record_acc_shape"))

hist(extract(sarn.fit2$fit)$record_acc_mean, freq = F)
lines(0:100, dgamma(0:100,
                  shape = sarn.fit2$data$record_prior_acc_mean_shape,
                  rate = sarn.fit2$data$record_prior_acc_mean_shape / sarn.fit2$data$record_prior_acc_mean_mean),
     type = "l", col = "red")


#shinystan::launch_shinystan(sarn.fit3$fit)


sarn.fit3 <- sarnie(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  K = 1000,
  K1 = 10,
  nu = 6,
  record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 1, chains = 3)

plot_stan_bacon(sarn.fit3, n.iter = 100, plot_priors = T)

summary(sarn.fit3$fit, par = c("record_acc_mean", "record_acc_shape"))$summary

traceplot(sarn.fit3$fit, par = c("record_acc_mean", "record_acc_shape"))



p.single.level <- plot_stan_bacon(sarn.fit1, 100, plot_priors = T)
p.multi.level <- plot_stan_bacon(sarn.fit2, 100, plot_priors = T)
p.multi.level.infl <- plot_stan_bacon(sarn.fit3, 100, plot_priors = T)



plot_acc_rate_pars(sarn.fit2)


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






