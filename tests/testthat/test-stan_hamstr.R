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
