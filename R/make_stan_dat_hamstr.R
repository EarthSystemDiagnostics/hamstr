#' Create K structure from K_tot and target_K_per_lvl
#'
#' @param K_tot total number of required sections at highest resolution
#' @param target_K_per_lvl approximate number of sections per level
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
#' @param K_fine total number of required sections at highest resolution
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
#' @keywords internal
#' @examples
#' default_K(100)
#' default_K(500)

default_K <- function(K_fine){

  bar <- function(x, y){
    abs(y - x^x)
  }

  #base <- round(optimize(bar, c(1, 10), y = K_fine)$minimum)
  
  base <- 2

  AdjustK(K_fine, base)

}


#' Calculate number of parameters being estimated for a given hierarchical structure
#'
#' @param base Number of new sections per per section
#' @param n Number of hierarchical levels
#'
#' @return
#' @keywords internal
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
#' @param mean Mean of gamma distribution
#' @param mode Mode of the gamma distribution
#' @param sigma Standard deviation of gamma distribution
#' @param shape Shape of gamma distribution
#'
#' @return a list
#' @keywords internal
#'
#' @examples
#' gamma_sigma_shape(mean = 10, sigma = 2)
gamma_sigma_shape <- function(mean = NULL, mode = NULL, sigma=NULL, shape=NULL){
  
  if (is.null(mean) & is.null(mode)) stop("One of either the mean or mode must be specified")
  if (is.null(shape) & is.null(sigma)) stop("One of either the shape or sigma must be specified")
  
  if (is.null(mean)==FALSE & is.null(mode)==FALSE) stop("Only one of the mean and mode can be specified")
  if (is.null(shape)==FALSE & is.null(sigma)==FALSE) stop("Only one of the shape and sigma can be specified")
  
  if (is.null(mean)){
    if (is.null(shape) == FALSE){
      if (shape <= 1) stop("Gamma cannot be defined by mode and shape if shape <= 1")
      mean <- (shape * mode) / (shape -1)
    } else if (is.null(shape)){
      
      # from mode and sigma
      rate = (mode + sqrt(mode^2 + 4*sigma^2)) / (2 * sigma^2)
      shape = 1 + mode * rate
      
      if (shape <= 1) stop("No solution for Gamma with this mode and sigma")
      mean <- (shape * mode) / (shape -1)
    }
  }
  
  if (is.null(sigma)){
    rate <- shape / mean
    sigma <- sqrt(shape/(rate^2))
  } else if (is.null(shape)){
    rate <- mean / sigma^2
    shape <- mean^2 / sigma^2
  }
  
  if (is.null(mode)){
   mode <- (shape - 1) / rate 
  }
  
  return(list(mean = mean, mode = mode, rate = rate, shape = shape, sigma = sigma))
}


#' Make the data object required by the Stan sampler
#'
#' @param ... Arguments passed from \code{\link{hamstr}}
#'
#' @return a list of data and parameters to be passed as data to the Stan sampler
#' @keywords internal
make_stan_dat_hamstr <- function(...) {

  # take the calling environment and create the data required for stan
  l <- c(as.list(parent.frame()))

  # get defaults
  default.args <- formals(hamstr)
  default.arg.nms <- names(default.args)


  # Overwrite the defaults with non-null passed arguments

  l <- l[lapply(l, is.null) == FALSE]

  default.args[names(l)] <- l

  l <- default.args
  
  
  # expand hamstr_control
  hc.default.args <- formals(hamstr_control)
  hc.default.arg.nms <- names(hc.default.args)

  hc <- l$hamstr_control
  hc <- hc[lapply(hc, is.null) == FALSE]

  hc.default.args[names(hc)] <- hc
  
  hc <- hc.default.args
  
  l <- append(l, hc)
  
  l <- l[names(l)!= "hamstr.control"]
  

 if (is.null(l$acc_mean_prior)){

    d <- data.frame(depth = l$depth, obs_age = l$obs_age)
    acc_mean <- stats::coef(MASS::rlm(obs_age~depth, data = d))[2]

    acc_mean <- signif(acc_mean, 2)
    
    # if negative replace with 20
    if (acc_mean <= 0) {
      warning("Estimated mean accumulation rate is negative - using value = 20")
      acc_mean <- 20
      }
    l$acc_mean_prior <- acc_mean
  }


  ord <- order(l$depth)

  l$depth <- l$depth[ord]
  l$obs_age <- l$obs_age[ord]
  l$obs_err <- l$obs_err[ord]


    if (l$model_bioturbation == TRUE){

      # check parameters
      if (is.null(l$L_prior_sigma) == FALSE)
        message("L_prior_shape is being overriden by L_prior_sigma.")

      if (length(c(l$L_prior_sigma, l$L_prior_shape)) == 0)
        stop("One of either L_prior_sigma or L_prior_shape must be specified.
             Set either to 0 to impose a fixed mixing depth.")

      if ((length(l$n_ind) == 1 | length(l$n_ind) == length(l$obs_age)) == FALSE)
        stop("n_ind must be either a single value or a vector the same length as obs_age")

      if (length(l$n_ind == 1)) l$n_ind <- rep(l$n_ind, length(l$obs_age))

      l$n_ind <- l$n_ind[ord]

      if (is.null(l$L_prior_sigma) == FALSE) {
        if (l$L_prior_sigma == 0) l$L_prior_shape <- 0 else
          l$L_prior_shape <- gamma_sigma_shape(mean = l$L_prior_mean,
                                               sigma = l$L_prior_sigma)$shape
      }

    } else if(l$model_bioturbation == FALSE){
      l$n_ind <- numeric(0)
    }


    if (is.null(l$infl_sigma_sd)){
      l$infl_sigma_sd <- 10 * mean(l$obs_err)
    }

    if (l$min_age > min(l$obs_age)) {
      warning("min_age is older than minimum obs_age")
    }

    # if (l$pad_top_bottom == TRUE){
    #   # Set start depth to 5% less than first depth observation, and DO allow negative depths
    #   depth_range <- diff(range(l$depth))
    #   buff <- 0.05 * depth_range
    # } else {
    #   buff <- 0
    # }

    if (is.null(l$top_depth)) l$top_depth <- l$depth[1] #- buff

    if (is.null(l$bottom_depth)) l$bottom_depth <- utils::tail(l$depth, 1) #+ buff

    depth_range <- l$bottom_depth - l$top_depth

    if(l$top_depth > min(l$depth)) stop("top_depth must be above or equal to the shallowest data point")
    if(l$bottom_depth < max(l$depth)) stop("bottom_depth must be deeper or equal to the deepest data point")


    if (is.null(l$K)){
      K_fine_1 <- l$bottom_depth - l$top_depth

      # set resolution so that there are only 10
      # sections between the median spaced 2 data points
      min.d.depth <- median(diff(sort(unique(l$depth))))
      K_fine_2 <- round(16 * K_fine_1 / min.d.depth )

      K_fine <- min(c(K_fine_1, K_fine_2))

      # prevent default values higher than 900
      if (K_fine > 900) K_fine <- 900

      l$K <- default_K(K_fine)
      }


    # Transformed arguments
    l$N <- length(l$depth)

    stopifnot(l$N == length(l$obs_err), l$N == length(l$obs_age))#, l$N == length(l$n_ind))

    alpha_idx <- alpha_indices(l$K)

    l$K_tot <- sum(alpha_idx$nK)
    l$K_fine <- utils::tail(alpha_idx$nK, 1)
    l$c <- 1:l$K_fine

    l$mem_alpha = l$mem_strength * l$mem_mean
    l$mem_beta = l$mem_strength * (1-l$mem_mean)

    l$mem_mean = l$mem_mean
    l$mem_strength = l$mem_strength

    l$delta_c = depth_range / l$K_fine
    l$c_depth_bottom = l$delta_c * l$c + l$top_depth
    l$c_depth_top = c(l$top_depth, l$c_depth_bottom[1:(l$K_fine-1)])

    l$modelled_depths <- c(l$c_depth_top[1], l$c_depth_bottom)

    # Index for which sections the target depth is in
    l$which_c = sapply(l$depth, function(d) which.max((l$c_depth_bottom < d) * (l$c_depth_bottom - d) ))

    l <- append(l, alpha_idx)

    l$n_lvls <- length(l$K)
    l$scale_shape = as.numeric(l$scale_shape)
    l$model_bioturbation = as.numeric(l$model_bioturbation)
    l$model_displacement = as.numeric(l$model_displacement)
    l$smooth_s = as.numeric(l$smooth_s)
    l$model_hiatus = as.numeric(l$model_hiatus)
    if (is.null(l$H_top)) l$H_top = l$top_depth
    if (is.null(l$H_bottom)) l$H_bottom = l$bottom_depth
    
    #l$K_idx <- l$lvl - 1

    l$smooth_i <- get_smooth_i(l, l$L_prior_mean)
    l$I <- nrow(l$smooth_i)

  return(l)
}


get_smooth_i <- function(d, w){
  
  w <- (w / d$delta_c )
  
  ri <- (-floor(w/2):floor(w/2))
  
  mi <- sapply(d$which_c, function(x) x + ri)
  
  mi[mi <= 0] <- abs(mi[mi <= 0]) +1  
  
  #mi[mi > d$K_fine] <- mi[mi > d$K_fine]-(2*(mi[mi > d$K_fine] - d$K_fine)-1)
  mi[mi > d$K_fine] <- 2*d$K_fine - mi[mi > d$K_fine] +1
  
  if (any(mi > d$K_fine)) stop("Acc rate smoothing index > K_fine")
  if (any(mi < 1)) stop("Acc rate smoothing index < 1")
  
  if (is.matrix(mi) == FALSE) mi <- rbind(mi)
  
  mi
  
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

  # robust lm to get age0
  d <- data.frame(depth = stan_dat$depth, obs_age = stan_dat$obs_age)
  rlm1 <- MASS::rlm(obs_age~depth, data = d)
  sigma <- summary(rlm1)$sigma
  
  
  l <- list(
    R = stats::runif(1, 0.1, 0.9),
    
    # create starting alpha values +- 3 SD from the overal prior mean (but always +ve)
    alpha = with(stan_dat, abs(stats::rnorm(K_tot, acc_mean_prior, acc_mean_prior/3))),
    #record_acc_mean = (abs(rnorm(1, stan_dat$acc_mean_prior, stan_dat$acc_mean_prior/3))),
    
    age0 = as.numeric(
      stats::predict(rlm1,
                     newdata = data.frame(depth = stan_dat$top_depth))
    ) + stats::rnorm(1, 0, sigma)
    
  )
  
  
  if (l$age0 < stan_dat$min_age) l$age0 <- stan_dat$min_age + abs(stats::rnorm(1, 0, 2))
  
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

