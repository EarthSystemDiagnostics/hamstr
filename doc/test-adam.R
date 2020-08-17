# test adam
library(tidyverse)
library(baconr)
library(rstan)

# devtools::document()
# pkgbuild::compile_dll(force = TRUE)
# devtools::load_all()
# devtools::install(quick = FALSE)
 

#Load and filter data
all.terr.dat <- readr::read_csv2("doc/Dating_Data.csv") %>%
  filter(complete.cases(age, depth, e.older)) %>%
  mutate(sigma.age = e.older)

all.terr.14C.dat <- all.terr.dat %>%
  filter(age.type == "Radiocarbon years BP") %>%
  group_by(DataName) %>%
  mutate(n.dates = n()) %>%
  ungroup()

name <- "MOKHOVOE"
#name <- "HANGING"

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


#options(mc.cores = parallel::detectCores())

options(mc.cores = 3)

adam.fit1 <- adam(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  K = c(100),
  nu = 6,
  #record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
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
  #record_prior_acc_mean_mean = 20,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  inflate_errors = 0, chains = 3)



ad_K10_100_rib <- plot_adam(adam.fit2, plot_diagnostics = T)
ad_K10_100_spag <- plot_adam(adam.fit2, type = "spaghetti", n.iter = 1000, plot_diagnostics = T)

ad_K10_100_spag_nod <- plot_adam(adam.fit2, type = "spaghetti", n.iter = 100, plot_diagnostics = F)

plot_hierarchical_acc_rate(adam.fit2)


ggsave("doc/ad_K10_100_rib.png", ad_K10_100_rib, width = 8, height = 6)
ggsave("doc/ad_K10_100_spag.png", ad_K10_100_spag, width = 8, height = 8)

ggsave("doc/ad_K10_100_spag_nod.png", ad_K10_100_spag_nod, width = 8, height = 6)


a2 <- rstan::summary(adam.fit2$fit)
a2 <- as_tibble(a2$summary, rownames = "par")


summary(adam.fit2$fit, par = c("alpha[1]", "shape"))$summary

traceplot(adam.fit2$fit, pars = c("shape", "alpha[1]"), inc_warmup = T)

traceplot(adam.fit2$fit, pars = paste0("alpha[", 91:111, "]"), inc_warmup = T)

traceplot(adam.fit2$fit, pars = c("shape", "alpha[1]"), inc_warmup = T)

traceplot(adam.fit2$fit, pars = c("w", "R"), inc_warmup = T)


## with 3 layers

adam.fit2b <- adam(
  depth = dat1$depth,
  obs_age = dat1$age.14C.cal,
  obs_err = dat1$age.14C.cal.se,
  K = baconr:::optimal_K(1000, 10),
  nu = 6,
  #record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
  record_prior_acc_shape_mean = 1.5,
  record_prior_acc_shape_shape = 1.5,
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

print(adam.fit3$fit, pars = c("infl_mean", "infl_shape", "infl_sigma", "infl"))
traceplot(adam.fit3$fit, pars = c("infl_mean", "infl_shape", "infl_sigma"), inc_warmup = T)


a3 <- rstan::summary(adam.fit3$fit)
a3 <- as_tibble(a3$summary, rownames = "par")

summary(adam.fit3$fit, par = c("alpha[1]", "shape"))$summary
summary(adam.fit3$fit, par = c("w"))$summary

traceplot(adam.fit3$fit, par = c("alpha[1]", "shape"))

traceplot(adam.fit3$fit, par = c("infl_mean", "infl_shape"))
traceplot(adam.fit3$fit, par = c("infl"))

summary(adam.fit3$fit, par = c("infl_mean", "infl_shape"))$summary


baconr:::plot_memory_prior_posterior(adam.fit3)



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
  #record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
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
  #record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
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
  #record_prior_acc_mean_mean = acc.mean,
  record_prior_acc_mean_shape = 1.5,
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

wa <- c(1:length(dat1$age.14C.cal))
wa <- c(1, 2, 5, 11)
p_age_pdf <- CompareAgePDF(dat1$age.14C.cal[wa], dat1$sigma.age[wa], curve = "intcal13", t.df = 6) + 
  facet_wrap(~.id, scales = "free") + theme_bw() + theme(panel.grid.minor = element_blank())

ggsave("doc/age_pdf.png", p_age_pdf, width = 6, height = 4, dpi = 300)



