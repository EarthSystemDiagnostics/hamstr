test_that("set seed works", {
  
  library(hamstr)
  
  hamstr_fit_1 <- hamstr(depth = 1:10,
                         obs_age = 1:10,
                         obs_err = rep(1, 10), 
                         # the seed argument for the sampler is set here so that
                         # this example always returns the same numerical result
                         stan_sampler_args = list(seed = 1, iter = 20, cores = 1))
  
  hamstr_fit_2 <- hamstr(depth = 1:10,
                         obs_age = 1:10,
                         obs_err = rep(1, 10), 
                         # the seed argument for the sampler is set here so that
                         # this example always returns the same numerical result
                         stan_sampler_args = list(seed = 1, iter = 20, cores = 1))
  
  
  testthat::expect_equal(summary(hamstr_fit_1)$n_eff, summary(hamstr_fit_2)$n_eff)
  
})


test_that("sample_posterior = FALSE works", {
  
  library(hamstr)
  
  hamstr_fit_1 <- hamstr(depth = 1:10,
                         obs_age = 1:10,
                         obs_err = rep(1, 10),
                         sample_posterior = FALSE)
  
  testthat::expect_true(is.na(hamstr_fit_1$fit))
  
  p1 <- plot(hamstr_fit_1, plot_diagnostics = FALSE)
  
  testthat::expect_equal(class(p1), c("gg", "ggplot"))
  
})



test_that("posterior and plotting functions work", {
  
  library(hamstr)
  
  hamstr_fit_1 <- hamstr(depth = 1:10,
                         obs_age = seq(1000, 10000, length.out = 10),
                         obs_err = rep(100, 10), 
                         # the seed argument for the sampler is set here so that
                         # this example always returns the same numerical result,
                         model_bioturbation = TRUE,
                         L_prior_shape = 2,
                         n_ind = rep(10, 10),
                         stan_sampler_args = list(seed = 1, iter = 20, cores = 1))
  
  
  # plotting age models
  p1 <- plot(hamstr_fit_1)
  p2 <- plot(hamstr_fit_1, "age")
  p3 <- plot(hamstr_fit_1, summarise = FALSE, n.iter = 5,
             plot_diagnostics = FALSE)
  
  testthat::expect_equal(class(p1), c("patchwork", "plot_filler", "gg", "ggplot"))
  testthat::expect_equal(class(p2), c("gg", "ggplot"))
  testthat::expect_equal(class(p3), c("gg", "ggplot"))
  
  # plotting accumulation rates
  p_acc <- plot(hamstr_fit_1, "acc_rates", tau = 2)
  p_h_acc <- plot(hamstr_fit_1, "hier_acc_rates")
  
  testthat::expect_equal(class(p_acc),  c("patchwork", "gg", "ggplot"))
  testthat::expect_equal(class(p_h_acc),  c("gg", "ggplot"))
 
  # plotting priors and posteriors
  p_acc_pr <- plot(hamstr_fit_1, "acc_mean_pr")
  p_L <- plot(hamstr_fit_1, "L")
  
  testthat::expect_equal(class(p_acc_pr),  c("gg", "ggplot"))
  testthat::expect_equal(class(p_L),  c("gg", "ggplot"))
  
  
  
  s1 <- summary(hamstr_fit_1)
  p1 <- predict(hamstr_fit_1)
  
  testthat::expect_equal(class(s1), c("tbl_df", "tbl", "data.frame"))
  testthat::expect_equal(class(p1), c("tbl_df", "tbl", "data.frame"))
  
  testthat::expect_equal(nrow(s1), 9)
  testthat::expect_equal(nrow(p1), 360)
  
  
  
  s1 <- summary(hamstr_fit_1, probs = c(0.23))
  
  spars <- summary(hamstr_fit_1, "pars")
  
  p1 <- predict(hamstr_fit_1, depth = c(3.4, 5.7))
  sp1 <- summary(p1)
  
  testthat::expect_equal(unique(p1$depth), c(3.4, 5.7))
  testthat::expect_equal(unique(sp1$depth), c(3.4, 5.7))
  
  testthat::expect_equal(names(s1),
                         c("depth", "idx", "par", "mean",
                           "se_mean", "sd", "23%", "n_eff", "Rhat"))
  
  
  testthat::expect_equal(spars$Parameter, c("alpha[1]", "R", "w",  "L[1]"))
  
  
  p_comp <- plot(hamstr_fit_1, type = "PDF", cal_curve = "marine20")
  
  
})

test_that("inflate_errors", {
  
  hamstr_fit_1 <- hamstr(depth = 1:10,
                         obs_age = seq(1000, 10000, length.out = 10),
                         obs_err = rep(100, 10), 
                         # the seed argument for the sampler is set here so that
                         # this example always returns the same numerical result,
                        hamstr_control = hamstr_control(inflate_errors = TRUE),
                        stan_sampler_args = list(seed = 1, iter = 20, cores = 1))
  
  p1 <- plot(hamstr_fit_1, plot_diagnostics = FALSE)
  
  p2 <- hamstr:::plot_infl_prior_posterior(hamstr_fit_1)
  
  testthat::expect_equal(class(p1), c("gg", "ggplot"))
  testthat::expect_equal(class(p2), c("patchwork", "plot_filler", "gg", "ggplot"))
  
})


test_that("displacement modelling works", {
  
  library(hamstr)
  
  hamstr_fit_1 <- hamstr(depth = 1:10,
                         obs_age = seq(1000, 10000, length.out = 10),
                         obs_err = rep(100, 10), 
                         # the seed argument for the sampler is set here so that
                         # this example always returns the same numerical result,
                         model_displacement = TRUE,
                         stan_sampler_args = list(seed = 1, iter = 20, cores = 1))
  
  p1 <- plot(hamstr_fit_1)
  p2 <- plot(hamstr_fit_1, "age")
  
  p_acc <- plot(hamstr_fit_1, "acc_rates")
  p_h_acc <- plot(hamstr_fit_1, "hier_acc_rates")
  
  p_D <- plot(hamstr_fit_1, "D")
  
  
  testthat::expect_equal(class(p1), c("patchwork", "plot_filler", "gg", "ggplot"))
  testthat::expect_equal(class(p2), c("gg", "ggplot"))
  
  testthat::expect_equal(class(p_acc),  c("patchwork", "gg", "ggplot"))
  
  testthat::expect_equal(class(p_h_acc),  c("gg", "ggplot"))
  testthat::expect_equal(class(p_D),  c("gg", "ggplot"))
  
  
  
  s1 <- summary(hamstr_fit_1)
  p1 <- predict(hamstr_fit_1)
  
  testthat::expect_equal(class(s1), c("tbl_df", "tbl", "data.frame"))
  testthat::expect_equal(class(p1), c("tbl_df", "tbl", "data.frame"))
  
  testthat::expect_equal(nrow(s1), 9)
  testthat::expect_equal(nrow(p1), 360)
  
})


testthat::test_that("compare_14C_PDF", {
  
  p_14C_PDF <- compare_14C_PDF(c(1000, 10000), c(20, 50))
  testthat::expect_equal(class(p_14C_PDF), c("gg", "ggplot"))
  
})