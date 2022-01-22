test_that("set seed works", {
  
  library(hamstr)
  
  hamstr_fit_1 <- hamstr(depth = 1:10,
                         obs_age = 1:10,
                         obs_err = rep(1, 10), 
                         # the seed argument for the sampler is set here so that
                         # this example always returns the same numerical result
                         stan_sampler_args = list(seed = 1, iter = 20))
  
  hamstr_fit_2 <- hamstr(depth = 1:10,
                         obs_age = 1:10,
                         obs_err = rep(1, 10), 
                         # the seed argument for the sampler is set here so that
                         # this example always returns the same numerical result
                         stan_sampler_args = list(seed = 1, iter = 20))
  
  
  testthat::expect_equal(summary(hamstr_fit_1)$n_eff, summary(hamstr_fit_2)$n_eff)
  
})



test_that("posterior and plotting functions work", {
  
  library(hamstr)
  
  hamstr_fit_1 <- hamstr(depth = 1:10,
                         obs_age = 1:10,
                         obs_err = rep(1, 10), 
                         # the seed argument for the sampler is set here so that
                         # this example always returns the same numerical result,
                         model_bioturbation = TRUE,
                         L_prior_shape = 0,
                         n_ind = rep(10, 10),
                         stan_sampler_args = list(seed = 1, iter = 20))
  
  p1 <- plot(hamstr_fit_1)
  p2 <- plot(hamstr_fit_1, "age")
  
  p_acc <- plot(hamstr_fit_1, "acc_rates")
  
  p_L <- plot(hamstr_fit_1, "L")
  
  
  testthat::expect_equal(class(p1), c("gg", "ggplot", "ggarrange"))
  testthat::expect_equal(class(p2), c("gg", "ggplot"))
  
  testthat::expect_equal(class(p_acc),  c("gg", "ggplot", "ggarrange"))
  testthat::expect_equal(class(p_L),  c("gg", "ggplot"))
  
  
  
  s1 <- summary(hamstr_fit_1)
  p1 <- predict(hamstr_fit_1)
  
  testthat::expect_equal(class(s1), c("tbl_df", "tbl", "data.frame"))
  testthat::expect_equal(class(p1), c("tbl_df", "tbl", "data.frame"))
  
  testthat::expect_equal(nrow(s1), 9)
  testthat::expect_equal(nrow(p1), 360)
  
  
})