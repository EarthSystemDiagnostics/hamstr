#' Create K structure from K_tot and target_K_per_lvl
#'
#' @param K_tot
#' @param target_K_per_lvl
#'
#' @return
#' @keywords internal
GetK <- function(K_tot, target_K_per_lvl = 10){

  K_per_lvl <- seq(floor(target_K_per_lvl/2), 2*target_K_per_lvl, by = 1)

  K_per_lvl <- K_per_lvl[K_per_lvl>1]

  n_lvls <- unlist(
    lapply(K_per_lvl,
           function(x) seq(floor(log(K_tot, base = x)/2),
                           ceiling(2*log(K_tot, base = x)))
    )
  )

  df <- expand.grid(K_per_lvl = unique(K_per_lvl), n_lvls = unique(n_lvls))

  df$K_fine <- with(df, K_per_lvl^n_lvls)

  idx <- which.min(abs(df$K_fine - K_tot))

  n_lvls <- df[idx, "n_lvls"]
  K_per_lvl <- df[idx, "K_per_lvl"]

  K <- rep(K_per_lvl, n_lvls)

  message(cat(cumprod(K)))

  return(K)
}


#' Adjust numbers of splits per level
#'
#' @param K_fine
#'
#' @return
#' @keywords internal
AdjustK <- function(K_fine, base){

  a <- log(K_fine, base)
  nLevels <- floor(a)

  K <- rep(base, nLevels)

  tot <- cumprod(K)[nLevels]

  i <- 0
  while(tot < K_fine){
    K[nLevels - i] <- K[nLevels - i] +1
    i <- i + 1
    if (i >= nLevels) i <- 0
    tot <- cumprod(K)[nLevels]
  }

  if (cumprod(K)[nLevels] > K_fine) {
    if (i == 0) i <- nLevels
    K[nLevels - i+1] <- K[nLevels - i+1] -1
  }

  K

}

#' Default K structure
#'
#' @param K_fine total number of sections at the finest resolution
#'
#' @return a vector
#' @export
#' @examples
#' default_K(100)
#' default_K(500)

default_K <- function(K_fine){

  bar <- function(x, y){
    abs(y - x^x)
  }

  base <- round(optimize(bar, c(1, 10), y = K_fine)$minimum)

  AdjustK(K_fine, base)

}


#' Calculate number of parameters being estimated for a given hierarchical structure
#'
#' @param base Number of new sections per per section
#' @param n Number of hierarchical levels
#'
#' @return
#' @export
#'
#' @examples
#' hierarchy_efficiency(10, 3)
hierarchy_efficiency <- function(base, n){

  nTot <- (base^1 - base^(n+1)) / (1-base)

  nFine <- base^n

  eff <- nTot / nFine

  return(c(nFine = nFine, nTot = nTot, eff = eff))

}


# Make index creating functions for K levels
#' Create alpha level indices
#'
#' @inheritParams hamstr
#'
#' @return a list
#' @keywords internal
alpha_indices <- function(K){

  K <- eval(K)

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

#' Convert between parametrisations of the gamma distribution
#'
#' @param mean 
#' @param sigma 
#' @param shape 
#'
#' @return a list
#' @keywords internal
#'
#' @examples
#' gamma_sigma_shape(mean = 10, sigma = 2)
gamma_sigma_shape <- function(mean, sigma=NULL, shape=NULL){
  if (is.null(sigma)){
    rate <- shape / mean
    sigma <- sqrt(shape/(rate^2))
  } else if (is.null(shape)){
    rate <- mean / sigma^2
    shape <- mean^2 / sigma^2
  }
  return(list(mean = mean, rate = rate, shape = shape, sigma = sigma))
}


#' Make the data object required by the Stan program
#'
#' @inheritParams hamstr
#'
#' @return a list of data and parameters to be passed as data to the Stan sampler
#' @export
#'
#' @examples
#' make_stan_dat_hamstr(depth = MSB2K$depth,
#'               obs_age = MSB2K$age,
#'               obs_err = MSB2K$error)
make_stan_dat_hamstr <- function(depth=NULL, obs_age=NULL, obs_err=NULL,
                                 n_ind = NULL,
                                 min_age = NULL,
                                 top_depth=NULL, bottom_depth=NULL,
                                 pad_top_bottom=NULL,
                                 K=NULL, nu=NULL,
                                 acc_mean_prior=NULL,
                                 acc_shape=NULL,
                                 scale_shape = NULL,
                                 mem_mean=NULL, mem_strength=NULL,
                                 scale_R=NULL,
                                 inflate_errors=NULL,
                                 infl_sigma_sd=NULL,
                                 infl_shape_shape=NULL, infl_shape_mean=NULL,
                                 model_bioturbation = NULL,
                                 L_prior_mean = NULL,
                                 L_prior_shape = NULL,
                                 L_prior_sigma = NULL,
                                 ...) {

  l <- c(as.list(environment()))

  # get defaults
  default.args <- formals(hamstr)
  default.arg.nms <- names(default.args)


  # Overwrite the defaults with non-null passed arguments

  l <- l[lapply(l, is.null) == FALSE]

  default.args[names(l)] <- l

  l <- default.args



  if (is.null(l$acc_mean_prior)){

    d <- data.frame(depth = depth, obs_age = obs_age)
    acc_mean <- stats::coef(MASS::rlm(obs_age~depth, data = d))[2]

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

  if (model_bioturbation == TRUE){
    
    # check parameters
      if (length(c(L_prior_sigma, L_prior_shape)) > 1) 
        stop("Specify only one of either L_prior_sigma or L_prior_shape. 
             The other will be calculated.")
      
      if (length(c(L_prior_sigma, L_prior_shape)) == 0) 
        stop("One of either L_prior_sigma or L_prior_shape must be specified. 
             Set either to 0 to impose a fixed mixing depth.")
      
    l$n_ind <- n_ind[ord]
    
    if (is.null(L_prior_shape)) {
      if (L_prior_sigma == 0) L_prior_shape <- 0 else 
        L_prior_shape <- gamma_sigma_shape(L_prior_mean, L_prior_sigma)$shape
    }
    
  } else if(model_bioturbation == FALSE){
    l$n_ind <- rep(1, length(l$obs_age))
  }


  if (is.null(infl_sigma_sd)){
    l$infl_sigma_sd <- 10 * mean(obs_err)
  }

  if (l$min_age > min(l$obs_age)) {
    warning("min_age is older than minimum obs_age")
  }

  if (l$pad_top_bottom == TRUE){
    # Set start depth to 5% less than first depth observation, and DO allow negative depths
    depth_range <- diff(range(l$depth))
    buff <- 0.05 * depth_range
  } else {
    buff <- 0
  }

  if (is.null(top_depth)) l$top_depth <- l$depth[1] - buff

  if (is.null(bottom_depth)) l$bottom_depth <- utils::tail(l$depth, 1) + buff

  depth_range <- l$bottom_depth - l$top_depth

  if(l$top_depth > min(l$depth)) stop("top_depth must be above or equal to the shallowest data point")
  if(l$bottom_depth < max(l$depth)) stop("bottom_depth must be deeper or equal to the deepest data point")


  if (is.null(K)){
    K_fine <- l$bottom_depth - l$top_depth
    if (K_fine > 900) K_fine <- 900
    l$K <- default_K(K_fine)
  }


  # Transformed arguments
  l$N <- length(l$depth)

  stopifnot(l$N == length(obs_err), l$N == length(obs_age), l$N == length(l$n_ind))

  alpha_idx <- alpha_indices(l$K)

  l$K_tot <- sum(alpha_idx$nK)
  l$K_fine <- utils::tail(alpha_idx$nK, 1)
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

  l$n_lvls <- length(l$K)
  l$scale_shape = as.numeric(l$scale_shape)
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


#' Create Random Initial Values for the hamstr Stan Model
#'
#' @param stan_dat stan data for the hamstr model
#'
#' @return a list of lists
#' @keywords internal
get_inits_hamstr <- function(stan_dat){

  l <- list(
    R = stats::runif(1, 0.1, 0.9),

    # create starting alpha values +- 3 SD from the overal prior mean (but always +ve)
    alpha = with(stan_dat, abs(stats::rnorm(K_tot, acc_mean_prior, acc_mean_prior/3))),
    #record_acc_mean = (abs(rnorm(1, stan_dat$acc_mean_prior, stan_dat$acc_mean_prior/3))),

    age0 = stats::rnorm(1, min(stan_dat$obs_age), 2)
  )

  # need to make this conditional and make sure initial values are arrays!
  if (stan_dat$inflate_errors == 1){
    l$infl_mean = as.array(abs(stats::rnorm(1, 0, 0.1)))
    l$infl_shape = as.array(1+abs(stats::rnorm(1, 0, 0.1)))
    l$infl = abs(stats::rnorm(stan_dat$N, 0, 0.1))
  } else {
    l$infl_mean = numeric(0)
    l$infl_sd = numeric(0)
    l$infl = numeric(0)
  }

  return(l)
}

