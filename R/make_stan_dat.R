#' make_stan_dat
#'
#' @inheritParams stan_bacon
#'
#' @return
#' @export
#'
#' @examples
make_stan_dat <- function(depth, obs_age, obs_err, K = 10, nu = 6,
                          acc_mean, acc_var = "default",
                          mem_mean = 0.7, mem_strength = 4, ...) {
  
  # Pretty sure nu = 6 is equivalent to the default parameterisation of t.a = 3, t.b = 4 in Bacon
  # when t.b = t.a +1, there are 2a degrees of freedom, and the distribution is symetric.
  
  l <- c(as.list(environment()), list(...))
  
  # Transformed arguments
  l$N <- length(l$depth)
  
  if (acc_var == "default") {
    l$acc_var = acc_mean^2 / 1.5}else{
      l$acc_var = acc_var
    }
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
