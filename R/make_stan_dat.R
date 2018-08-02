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


#####

make_stan_dat_hi <- function(depth, obs_age, obs_err,
                             hiatus_depth = NULL, hiatus_length = NULL,
                             hiatus_shape = 1, hiatus_interval = 0.1,
                             K = 10, nu = 6,
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
 
  # Hiatus data
  l$hiatus_depth = as.array(hiatus_depth)
  l$hiatus_length = as.array(hiatus_length)
  
  l$N_hiatus = length(hiatus_depth)
  l$hiatus_shape = hiatus_shape
  l$hiatus_beta = as.array(hiatus_shape / hiatus_length)
  
  stopifnot(length(hiatus_depth) == length(hiatus_length))
  if (any(hiatus_depth < min(depth))) stop("One or more hiatus is shallower than the minimum core depth")
  if (any(hiatus_depth > max(depth))) stop("One or more hiatus is deeper than the maximum core depth")
  
  # Set start depth to 5% less than first depth observation, but do not allow negative depths
  depth_range <- diff(range(depth))
  buff_5 <- 0.05 * depth_range
  strt_dpth <- depth[1] - buff_5
  strt_dpth[strt_dpth < 0] <- 0
  end_dpth <- tail(depth, 1) + buff_5
  depth_range = end_dpth - strt_dpth
  #l$delta_c = depth_range / K
  
  pos <- sort(c(strt_dpth, hiatus_depth, hiatus_depth+hiatus_interval, end_dpth))
  dfs <- diff(pos)
  prps <- dfs / sum(dfs)
  brks <- ceiling(prps * K)
  #print(pos)
  #print(brks)
  dpths <- unique(unlist(sapply(1:(length(brks)), function(i) seq(pos[i], pos[i+1], length.out = brks[i]))))
  #print(dpths)
  l$c_depth_top = dpths[1:(length(dpths)-1)]
  l$c_depth_bottom = dpths[2:length(dpths)]
  
  l$K = length(l$c_depth_bottom)
  
  l$delta_c = diff(dpths)
  
  # Index for which sections the target depth is in
  l$which_c = sapply(l$depth, function(d) which.max((l$c_depth_bottom < d) * (l$c_depth_bottom - d) ))
  
  # Index for which sections the hiatuses are in
  l$which_c_hiatus = as.array(sapply(l$hiatus_depth, function(d) which.max((l$c_depth_top < d) * (l$c_depth_top - d))))
  l$hiatus_factor = rep(1, l$K)
  l$hiatus_factor[l$which_c_hiatus] <- 0
  
  return(l)
} 

# hi_ages <- MSB2K$age
# hi_ages[MSB2K$depth > 40] <- hi_ages[MSB2K$depth > 40] + 1000 
# 
# hi_dat <- make_stan_dat_hi(depth = MSB2K$depth,
#                            obs_age = hi_ages,
#                            obs_err = MSB2K$error,
#                            hiatus_depth = c(39.5),
#                            hiatus_length = c(1000),
#                            K = 25, nu = 6,
#                            acc_mean = 20, 
#                            mem_mean = 0.7, mem_strength = 4)
# 
# 
# fit <- rstan::stan(file = "exec/bacon_hiatuses.stan", 
#                    data = hi_dat, iter = 1000, chains = 3,
#                    verbose = FALSE)
# 
# 
# fit_obj <- list(fit=fit, data=hi_dat)
# 
# plot_stan_bacon(fit_obj, plot_priors = T)# + geom_vline(xintercept = hi_dat$c_depth_top)

