#' stan_bacon_dev
#' @md
#' @useDynLib baconr, .registration = TRUE
#' @exportPattern("^[[:alpha:]]+")
#' @inheritParams rstan::sampling
#' @inheritParams baconr::stan_bacon
#'   
#' @return Returns a list composed of the output from the Stan sampler .$fit,
#'   and the list of data passed to the sampler
#' @export
#' 
#' @examples
stan_bacon <- function(depth, obs_age, obs_err, K = 10, nu = 6,
                  acc_mean, acc_alpha = 1.5,
                  mem_mean = 0.7, mem_strength = 4,
                  iter = 2000, chains = 4){
  
  stan_dat <- make_stan_dat(depth, obs_age, obs_err, K, nu,
                         acc_mean, acc_alpha,
                         mem_mean, mem_strength)
  
  fit <- rstan::sampling(stanmodels$bacon_dev, 
              data = stan_dat, iter = iter, chains = chains,
              verbose = FALSE)
  
  return(list(fit=fit, data=stan_dat))

  }
