#' make_stan_dat_sarnie
#'
#' @inheritParams sarnie
#'
#' @return a list of data, and parameters to be passed as data to the Stan sampler
#' @export
#'
#' @examples
#' make_stan_dat_sarnie(depth = MSB2K$depth,
#'               obs_age = MSB2K$age,
#'               obs_err = MSB2K$error,
#'               K1 = 10, K = 100,
#'               nu = 6)
make_stan_dat_sarnie <- function(depth, obs_age, obs_err,
                             K1 = 10, K = 100, nu = 6,
                             record_prior_acc_mean_mean = 20,
                             record_prior_acc_mean_shape = 1.5,
                             record_prior_acc_shape_mean = 1.5,
                             record_prior_acc_shape_shape = 1.5,
                             mem_mean = 0.7, mem_strength = 4,
                             inflate_errors = 0) {
  
  
  l <- c(as.list(environment()))
  
  # make sure in depth order
  ord <- order(l$depth)
  l$depth <- l$depth[ord]
  l$obs_age <- l$obs_age[ord]
  l$obs_err <- l$obs_err[ord]
  
  # Transformed arguments
  l$N <- length(l$depth)
  
  stopifnot(K%%K1 == 0)
  l$whichK1 = rep(1:K1, each = K / K1)
  
  l$mem_alpha = mem_strength * mem_mean
  l$mem_beta = mem_strength * (1-mem_mean)
  
  l$mem_mean = mem_mean
  l$mem_strength = mem_strength
  
  # Set start depth to 5% less than first depth observation, and DO allow negative depths
  depth_range <- diff(range(l$depth))
  buff_5 <- 0.05 * depth_range
  strt_dpth <- l$depth[1] - buff_5
  #strt_dpth[strt_dpth < 0] <- 0
  end_dpth <- tail(l$depth, 1) + buff_5
  depth_range = end_dpth - strt_dpth
  
  
  l$c = 1:K
  l$delta_c = depth_range / K
  l$c_depth_bottom = l$delta_c * l$c + strt_dpth
  l$c_depth_top = c(strt_dpth, l$c_depth_bottom[1:(K-1)])
  
  # Index for which sections the target depth is in
  l$which_c = sapply(l$depth, function(d) which.max((l$c_depth_bottom < d) * (l$c_depth_bottom - d) ))
  
  return(l)
}


get_inits_sarnie <- function(dat){
  l <- list(
    R = runif(1, 0.1, 0.9),
    alpha = (abs(rnorm(dat$K, dat$record_prior_acc_mean_mean, dat$record_prior_acc_mean_mean/3))),
    record_acc_mean = (abs(rnorm(1, dat$record_prior_acc_mean_mean, dat$record_prior_acc_mean_mean/3))),
    section_acc_mean = as.array(abs(rnorm(dat$K1, dat$record_prior_acc_mean_mean, dat$record_prior_acc_mean_mean/3))),
    
    record_acc_shape = rnorm(1, 1.5, 1.5/3),
    
    age0 = rnorm(1, min(dat$obs_age), 2)
  )
  
  # need to make this conditional and make sure initial values are arrays! 
  if (dat$inflate_errors == 1){
    l$infl_mean = as.array(abs(rnorm(1, 0, 0.1)))
    l$infl_sd = as.array(abs(rnorm(1, 0, 0.1)))
    l$infl = abs(rnorm(dat$N, 0, 0.1))   
  } else {
    l$infl_mean = numeric(0)
    l$infl_sd = numeric(0)
    l$infl = numeric(0)   
  }
  
  return(l)
  
}
