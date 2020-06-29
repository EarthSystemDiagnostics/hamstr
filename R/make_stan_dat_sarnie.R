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
#'               
# acc_mean_hyp = record_prior_acc_mean_mean 
# acc_shape_hyp = record_prior_acc_mean_shape
# record_acc_mean_shape_mean 
# record_prior_acc_shape_shape

make_stan_dat_sarnie <- function(depth, obs_age, obs_err,
                                 top_depth = NULL, bottom_depth = NULL,
                                 pad_top_bottom = FALSE,
                                 K1 = 10, K = 100, nu = 6,
                                 record_prior_acc_mean_mean = 20,
                                 record_prior_acc_mean_shape = 1.5,
                                 record_prior_acc_shape_mean = 1.5,
                                 record_prior_acc_shape_shape = 1.5,
                                 section_acc_shape = 1.5,
                                 mem_mean = 0.7, mem_strength = 4,
                                 inflate_errors = 0) {
  
  l <- c(as.list(environment()))
  
  ord <- order(depth)
  
  l$depth <- depth[ord]
  l$obs_age <- obs_age[ord]
  l$obs_err <- obs_err[ord]
  
  if (pad_top_bottom == TRUE){
     # Set start depth to 5% less than first depth observation, and DO allow negative depths
    depth_range <- diff(range(l$depth))
    buff <- 0.05 * depth_range
  } else {
    buff <- 0
  }
 
  if (is.null(top_depth)) l$top_depth <- l$depth[1] - buff
  
  if (is.null(bottom_depth)) l$bottom_depth <- tail(l$depth, 1) + buff
  
  depth_range <- l$bottom_depth - l$top_depth
  
  if(l$top_depth > min(l$depth)) stop("top_depth must be above or equal to the shallowest data point")
  if(l$bottom_depth < max(l$depth)) stop("bottom_depth must be deeper or equal to the deepest data point")
  
  
  # Transformed arguments
  l$N <- length(l$depth)
  
  stopifnot(l$N == length(obs_err), l$N == length(obs_age))
  
  l$K <- K
  l$c <- 1:l$K
  
  l$K1 <- K1
  
  if (K1 == 1){
    l$whichK1 <- rep(1, K)} else {
      l$whichK1 = cut(seq_along(1:K), K1, labels = FALSE)
      }
  
  l$mem_alpha = mem_strength * mem_mean
  l$mem_beta = mem_strength * (1-mem_mean)
  
  l$mem_mean = mem_mean
  l$mem_strength = mem_strength
  
  l$delta_c = depth_range / l$K
  l$c_depth_bottom = l$delta_c * l$c + l$top_depth
  l$c_depth_top = c(l$top_depth, l$c_depth_bottom[1:(l$K-1)])
  
  l$modelled_depths <- c(l$c_depth_top[1], l$c_depth_bottom)
  
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
