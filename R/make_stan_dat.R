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
                       mem_mean = 0.7, mem_strength = 4) {
  
  # Pretty sure nu = 6 is equivalent to the default parameterisation of t.a = 3, t.b = 4 in Bacon
  # when t.b = t.a +1, there are 2a degrees of freedom, and the distribution is symetric.
  
  l <- list(depth=depth, obs_age=obs_age, obs_err=obs_err)
  l$N = length(depth)
  l$K = K
  l$nu = nu
  l$acc_mean = acc_mean
  if (acc_var == "default") {
    l$acc_var = acc_mean^2 / 1.5}else{
      l$acc_var = acc_var
    }
  # l$mem_alpha = mem_strength * mem_mean
  # l$mem_beta = mem_strength * (1-mem_mean)
  l$mem_mean = mem_mean
  l$mem_strength = mem_strength
  l$c = 1:K
  l$delta_c = diff(range(depth)) / K
  l$c_depth_bottom = l$delta_c * l$c + depth[1]
  l$c_depth_top = c(0, l$c_depth_bottom[1:(K-1)])
  l$which_c = sapply(l$depth, function(d) which.max((l$c_depth_bottom < d) * (l$c_depth_bottom - d) ))
  return(l)
} 
