#' make_stan_dat
#'
#' @inheritParams stan_bacon
#'
#' @return a list of data, and parameters to be passed as data to the Stan sampler
#' @export
#'
#' @examples
#' # Get number of sections K, so that they will be ~ 5cm
#' K_for_5cm <- round(diff(range(MSB2K$depth)) / 5)
#' 
#' make_stan_dat(depth = MSB2K$depth,
#'               obs_age = MSB2K$age,
#'               obs_err = MSB2K$error,
#'               K = K_for_5cm, nu = 6,
#'               acc_mean = 20, acc_var = "default",
#'               mem_mean = 0.7, mem_strength = 4)


make_stan_dat <- function(depth, obs_age, obs_err, K = 10, nu = 6,
                          acc_mean, acc_alpha = 1.5,
                          mem_mean = 0.7, mem_strength = 4, ...) {
  
  # Pretty sure nu = 6 is equivalent to the default parameterisation of t.a = 3, t.b = 4 in Bacon
  # when t.b = t.a +1, there are 2a degrees of freedom, and the distribution is symetric.
  
  l <- c(as.list(environment()), list(...))
  
  # Transformed arguments
  l$N <- length(l$depth)
  
  l$acc_beta =  acc_alpha / acc_mean
  
  l$mem_alpha = mem_strength * mem_mean
  l$mem_beta = mem_strength * (1-mem_mean)
  
  l$mem_mean = mem_mean
  l$mem_strength = mem_strength
  l$c = 1:K
  
  
  # Set start depth to 5% less than first depth observation, but do not allow negative depths
  depth_range <- diff(range(depth))
  buff_5 <- 0.05 * depth_range
  strt_dpth <- depth[1] - buff_5
  strt_dpth[strt_dpth < 0] <- 0
  end_dpth <- tail(depth, 1) + buff_5
  depth_range = end_dpth - strt_dpth
  l$delta_c = depth_range / K
  l$c_depth_bottom = l$delta_c * l$c + strt_dpth
  l$c_depth_top = c(strt_dpth, l$c_depth_bottom[1:(K-1)])
  
  # Index for which sections the target depth is in
  l$which_c = sapply(l$depth, function(d) which.max((l$c_depth_bottom < d) * (l$c_depth_bottom - d) ))
  return(l)
} 
