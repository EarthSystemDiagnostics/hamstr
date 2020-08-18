
#' Estimate Optimal Hierarchical Structure
#'
#' @param K_tot 
#'
#' @return
#' @keywords internal
optimal_K <- function(K_tot, target_K_per_lvl = 10){
 
   n_lvls <- ceiling(log(K_tot, base = target_K_per_lvl))
   K_per_lvl <- round(K_tot^(1/n_lvls))
   
   K <- rep(K_per_lvl, n_lvls)
  
   return(K)
}


# Make index creating functions for K levels
#' Create alpha level indices
#'
#' @inheritParams adam
#'
#' @return a list
#' @keywords internal
alpha_indices <- function(K){

  # prepend 1 for the single overall mean alpha 
  K <- c(1, K)

  # number of sections at each level
  nK <- cumprod(K)

  
  K <- c(0, K)
  nK <- c(0, nK)

  alpha_idx <- 1:sum(nK)

  nLevels <- length(K)-1

  # which level is each parameter
  lvl <- unlist(lapply(seq_along(nK[-1]), function(i) rep(i, times = nK[i+1])))

  # index the parent of each parameter
  parent <- c(rep(0, K[2]), unlist(lapply(alpha_idx[1:sum(nK[1:nLevels])], function(i) rep(i, K[lvl[i]+2]))))

  list(alpha_idx=alpha_idx, lvl=lvl, parent=parent, nK = nK[-1])
}


#' Make the data object required by the Stan program
#'
#' @inheritParams adam
#'
#' @return a list of data and parameters to be passed as data to the Stan sampler
#' @export
#'
#' @examples
#' make_stan_dat_adam(depth = MSB2K$depth,
#'               obs_age = MSB2K$age,
#'               obs_err = MSB2K$error,
#'               nu = 6)
make_stan_dat_adam <- function(depth, obs_age, obs_err,
                               top_depth = NULL, bottom_depth = NULL,
                                 pad_top_bottom = FALSE,
                                 K = c(10, 10), nu = 6,
                                 acc_mean_prior = NULL,
                               shape = 1.5,
                                 #record_prior_acc_shape_mean = 1.5,
                                 #record_prior_acc_shape_shape = 1.5,
                                 # section_acc_shape = 1.5,
                                 mem_mean = 0.7, mem_strength = 4,
                                 inflate_errors = 0) {

  l <- c(as.list(environment()))
  
  if (is.null(acc_mean_prior)){
    
    d <- data.frame(depth = depth, obs_age = obs_age)
    acc_mean <- coef(MASS::rlm(obs_age~depth, data = d))[2]
    
    # if negative replace with 20
    if (acc_mean <= 0) {
      warning("Estimated mean accumulation rate is negative - using value = 20")
      acc_mean <- 20
      }
    l$acc_mean_prior <- acc_mean
    
  }

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

  alpha_idx <- alpha_indices(K)

  l$K_tot <- sum(alpha_idx$nK)
  l$K_fine <- tail(alpha_idx$nK, 1)
  l$c <- 1:l$K_fine

  l$mem_alpha = mem_strength * mem_mean
  l$mem_beta = mem_strength * (1-mem_mean)

  l$mem_mean = mem_mean
  l$mem_strength = mem_strength

  l$delta_c = depth_range / l$K_fine
  l$c_depth_bottom = l$delta_c * l$c + l$top_depth
  l$c_depth_top = c(l$top_depth, l$c_depth_bottom[1:(l$K_fine-1)])

  l$modelled_depths <- c(l$c_depth_top[1], l$c_depth_bottom)

  # Index for which sections the target depth is in
  l$which_c = sapply(l$depth, function(d) which.max((l$c_depth_bottom < d) * (l$c_depth_bottom - d) ))

  l <- append(l, alpha_idx)
  
  #l$K_lvls <- length(l$K)
  #l$K_idx <- l$lvl - 1

  return(l)
}


#' Calculated depth of section boundary at all hierarchical levels
#'
#' @param stan_dat 
#'
#' @return
#' @keywords internal
hierarchical_depths <- function(stan_dat){
  d_range <- diff(range(stan_dat$modelled_depths))
  min_d <- min(stan_dat$modelled_depths)
  
  delta_d <- d_range / stan_dat$nK
  
  lapply(stan_dat$nK[-1], function(x) {
    delta_d <- d_range / x
    c(min_d, delta_d * 1:x + min_d)
  })
}



#' Create Random Initial Values for the adam Stan Model
#'
#' @param stan_dat stan data for the adam model
#'
#' @return a list of lists
#' @keywords internal
get_inits_adam <- function(stan_dat){
  
  l <- list(
    R = runif(1, 0.1, 0.9),
    
    # create starting alpha values +- 3 SD from the overal prior mean (but always +ve)
    alpha = with(stan_dat, abs(rnorm(K_tot, acc_mean_prior, acc_mean_prior/3))),
    #record_acc_mean = (abs(rnorm(1, stan_dat$acc_mean_prior, stan_dat$acc_mean_prior/3))),

    #shape = abs(rnorm(1, 1.5, 1.5/3)),

    age0 = rnorm(1, min(stan_dat$obs_age), 2)
  )

  # need to make this conditional and make sure initial values are arrays!
  if (stan_dat$inflate_errors == 1){
    l$infl_mean = as.array(abs(rnorm(1, 0, 0.1)))
    l$infl_sd = as.array(abs(rnorm(1, 0, 0.1)))
    l$infl = abs(rnorm(stan_dat$N, 0, 0.1))
  } else {
    l$infl_mean = numeric(0)
    l$infl_sd = numeric(0)
    l$infl = numeric(0)
  }

  return(l)
}

