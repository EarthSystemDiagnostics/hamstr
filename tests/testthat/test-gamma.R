test_that("gamma_sigma_shape works", {

  library(hamstr)

  testthat::expect_equal(gamma_sigma_shape(mean = 10, sigma = 2)$shape, 25)
  testthat::expect_equal(gamma_sigma_shape(mean = 10, shape = 25)$sigma, 2)
  testthat::expect_equal(gamma_sigma_shape(mean = 1, sigma = 25)$mode, 0)
  testthat::expect_equal(gamma_sigma_shape(mode = 1, shape = 2)$mean, 2)

  testthat::expect_error(gamma_sigma_shape(mode = 0, sigma = 25))
  testthat::expect_length(gamma_sigma_shape(shape = 10, mean = 25), 5)

})

