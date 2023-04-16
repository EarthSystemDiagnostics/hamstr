# new K_fine and parent calculation
#library(tidyverse)




GetWts <- function(a, b){
  
  intvls <- lapply(seq_along(b)[-1], function(i) {
    b[c(i-1, i)]
  })
  
  gaps <- lapply(intvls, function(x) {
    a[a >= x[1] & a <= x[2]]
  })
  
  wts <- sapply(seq_along(intvls), function(i) {
    
    wts <- diff(unique(sort(c(gaps[[i]], intvls[[i]])))) / diff(range(sort(c(gaps[[i]], intvls[[i]]))))
    
    if (length(wts) == 1){rep(wts, 2)}else{wts}
    
  })
  wts
}


GetIndices <- function(nK = NULL, brks = NULL) {
  
  if (is.null(brks)){
    
    lvl <- unlist(lapply(seq_along(nK), function(i) rep(i, times = nK[i])))
    brks <- (lapply(nK, function(x) seq(0, 1, length.out = x + 1)))
  
  }
  
  if (is.null(nK)){
    nK <- sapply(brks, length)-1
    lvl <- unlist(lapply(seq_along(nK), function(i) rep(i, times = nK[i])))
  }
  

  # get the left hand parent
  parenta <- lapply(seq_along(brks)[-1], function(x) {
    sapply(1:(length(brks[[x]]) - 1), function(i) {
      cut(brks[[x]][c(i)], brks[[x - 1]],
          include.lowest = TRUE,
          right = FALSE, 
          labels = FALSE
      )
    })
  })
  
  # get the right hand parent
  parentb <- lapply(seq_along(brks)[-1], function(x) {
    sapply(1:(length(brks[[x]]) - 1), function(i) {
      cut(brks[[x]][c(i+1)], brks[[x - 1]],
          include.lowest = TRUE,
          right = TRUE, 
          labels = FALSE
      )
    })
  })
  
  
  parent1 <- lapply(seq_along(parenta), function(i){
    rbind(parenta[[i]], parentb[[i]])
  })
  
  parent <- append(list(rbind(0,0)), parent1)
  
  cumNcol <- cumsum(sapply(seq_along(parent)[-1], function(i) max(parent[[i-1]])))
  
  # shift the index to account for previous levels
  parent <- lapply(seq_along(parent)[-1], function(x) {
    parent[[x]] + cumNcol[[x-1]] 
  })
  
  # get the left/right weights
  wts <- lapply(seq_along(brks)[-1], function(i){
    GetWts(brks[[i-1]], brks[[i]])
  })
  
 # collapse to single matrix
 parent <- do.call(cbind, parent)
 
 multi_parent_adj <- mean(apply(parent, 2, function(x) abs(diff(x))+1))
 
 wts <- do.call(cbind, wts)
 
 # make weights add to 1
 wts <- apply(wts, MARGIN = 2, function(x) x / sum(x)) 
  
  list(nK = nK,
       #K = sapply(brks, length),
       alpha_idx = 1:sum(nK),
       lvl = lvl, brks = brks,
       multi_parent_adj = multi_parent_adj,
       parent1 = as.numeric(parent[1,]),
       parent2 = as.numeric(parent[2,]),
       wts1 = as.numeric(wts[1,]),
       wts2 = as.numeric(wts[2,])
       )
}


GetBrksHalfOffset <- function(K_fine, K_factor = 2){
  
  db_fine <- 1 / K_fine
  db <- db_fine
  
  brks <- list(
    seq(0, 1, by = db)
  )
  
  n_br <- length(brks[[1]])
  n_sec <- n_br - 1
  
  newbrks <- brks[[1]]
  
  while (n_sec > 3){
    
    strt <- min(brks[[length(brks)]])
    end <- max(brks[[length(brks)]])
    
    n_new <- ceiling((n_sec+1)/K_factor)
     
    # in units of old sections
    l_new <- n_new * K_factor
    l_old <- n_sec 
   
    d_new_old <- l_new - l_old
    
    if (d_new_old %% 2 == 0){
      new_strt <- strt - db * (d_new_old-1)/ 2
    } else {
      new_strt <- strt - db * (d_new_old)/ 2
    }
    
    
    newbrks <- seq(new_strt, by = db*K_factor, length.out = n_new+1)
    
    #newbrks <- unique(newbrks)
    brks <- c(brks, list(newbrks))
    
    db <- K_factor * db
    n_br <- length(newbrks)
    n_sec <- n_br -1
   
    #if (K_factor == n_sec){break} 
   
  }
  
  brks <- c(brks, list(c(newbrks[1], tail(newbrks, 1))))
  brks <- rev(brks)
  
  return(brks)
}

# #debugonce(GetIndices)
# 
# debugonce(GetBrksHalfOffset)
# 
#brks <- GetBrksHalfOffset(K_fine = c(99), K_factor = 2)
#brks
#brks <- rev(brks)

#GetIndices(brks = brks)
#GetIndices(nK = c(2,3,5))


# #debugonce(GetIndices)
# nK <- c(1, 3,5)
# GetIndices(nK)



#hierarchical_depths2(fit_HP2$data)

#' #' Get Nearest Primes to Geometric Series
#' #'
#' #' @param m Terminal number
#' #' @param f Approximate ratio of series
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' #' GetNextMultiPrimes(100, 3)
#' GetNextMultiPrimes <- function(m, f = 2){
#'   x <- 1
#'   np <- 1
#'   
#'   while(np < m){
#'     np <- as.numeric(gmp::nextprime(f*tail(x, 1)-1))
#'     x <- append(x, np)
#'   }
#'   
#'   x <- x[x <= m]
#'   
#'   if (max(x) < 0.5*m){
#'    x <- append(x, m)
#'   }else{
#'     x <- append(head(x, -1), m)
#'   }
#'   
#'   return(x)
#' }


#' Make the data object required by the Stan sampler
#'
#' @param ... Arguments passed from \code{\link{hamstr}}
#'
#' @return a list of data and parameters to be passed as data to the Stan sampler
#' @keywords internal
make_stan_dat_hamstr_prime <- function(...) {
  
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
  
  
  if (is.null(l$top_depth)) l$top_depth <- l$depth[1] #- buff
  
  if (is.null(l$bottom_depth)) l$bottom_depth <- utils::tail(l$depth, 1) #+ buff
  
  depth_range <- l$bottom_depth - l$top_depth
  
  if(l$top_depth > min(l$depth)) stop("top_depth must be above or equal to the shallowest data point")
  if(l$bottom_depth < max(l$depth)) stop("bottom_depth must be deeper or equal to the deepest data point")
  
  
  if (is.null(l$K_fine)){
    
    K_fine_1 <- l$bottom_depth - l$top_depth
    
    # set resolution so that there are only 16
    # sections between the median spaced 2 data points
    median.d.depth <- stats::median(diff(sort(unique(l$depth))))
    K_fine_2 <- round(16 * K_fine_1 / median.d.depth )
    
    K_fine <- min(c(K_fine_1, K_fine_2))
    
    # prevent default values higher than 900
    if (K_fine > 900) K_fine <- 900
    
    l$K_fine <- K_fine
  }
  
  
  if (is.null(l$K_factor)){
    
    l$K_factor <- get_K_factor(l$K_fine)
    
  }
  
  
  # Transformed arguments
  l$N <- length(l$depth)
  
  stopifnot(l$N == length(l$obs_err), l$N == length(l$obs_age))
  
  brks <- GetBrksHalfOffset(K_fine = l$K_fine, K_factor = l$K_factor)
  
  alpha_idx <- GetIndices(brks = brks)
  
  # Shannon <- function(x){
  #   
  #   x <- x[x>0]
  #   N = sum(x)
  #   p = x/N
  #   H = -sum(p*log(p))
  #   H
  # }
  # 
  # EffN <- function(x){
  #   
  #   exp(Shannon(x))
  #   
  # }
   
  
  #l$acc_shape <- l$acc_shape / shp_adj
  
  
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
  
  l$n_lvls <- length(l$nK) -1
  l$scale_shape = as.numeric(l$scale_shape)
  l$model_bioturbation = as.numeric(l$model_bioturbation)
  l$model_displacement = as.numeric(l$model_displacement)
  l$smooth_s = as.numeric(l$smooth_s)
  l$model_hiatus = as.numeric(l$model_hiatus)
  if (is.null(l$H_top)) l$H_top = l$top_depth
  if (is.null(l$H_bottom)) l$H_bottom = l$bottom_depth
  

  if (l$smooth_s == 1){
    l$smooth_i <- get_smooth_i(l, l$L_prior_mean)
    l$I <- nrow(l$smooth_i)
  } 
  
  else {
    l$smooth_i <- rbind(rep(1, l$N))
    l$I <- 1
  }
  
  return(l)
}

#' @title Fit a hamstr Age-Depth Model
#' @description \code{hamstr} is used to fit an age-depth model to a set of
#'   age-control points. Ages should already be on the desired scale, e.g.
#'   calendar ages, and will not be calibrated. The function
#'   \code{\link{calibrate_14C_age}} can be used to calibrate radiocarbon dates prior
#'   to fitting a hamstr model.
#' @param depth depths of observed ages (age control points)
#' @param obs_age observed age at each depth (age control points)
#' @param obs_err error associated with each observed age (1 standard error)
#' @param min_age the minimum age that the first modelled depth can be. Useful
#'   if extrapolating above the shallowest age control point to e.g. the
#'   surface. So set min_age to the year the core was collected. E.g. for a core
#'   collected in 1990, with ages in years BP this would be -40 (present = 1950
#'   by convention). The default value is the current year in calendar age BP
#'   (e.g. -71 for 2021).
#' @param top_depth,bottom_depth the top and bottom depths of the desired
#'   age-depth model. Must encompass the range of the data. Defaults to the
#'   shallowest and deepest data points.
#' @param K_fine the number of sections at the highest resolution of the model.
#' @param K_factor the rate at which the thickness of the sections grows between
#' subsequent levels. 
#' @param acc_mean_prior hyperparameter for the prior on the overall mean
#'   accumulation rate for the record. Units are obs_age / depth. E.g. if depth
#'   is in cm and age in years then the accumulation rate is in years/cm. The
#'   overall mean accumulation rate is given a weak half-normal prior with mean
#'   = 0, SD = 10 * acc_mean_prior. If left blank, acc_mean_prior is set to the
#'   mean accumulation rate estimated by fitting a robust linear model using
#'   \link[MASS]{rlm}.
#' @param acc_shape hyperparameter for the shape of the priors on accumulation
#'   rates. Defaults to 1.5 - as for Bacon 2.2.
#' @param mem_mean hyperparameter; a parameter of the Beta prior distribution on
#'   "memory", i.e. the autocorrelation parameter in the underlying AR1 model.
#'   The prior on the correlation between layers is scaled according to the
#'   thickness of the sediment sections in the highest resolution hierarchical
#'   layer, *delta_c*, which is determined by the total length age-models and
#'   the parameter vector *K*. mem_mean sets the mean value for *R* (defaults to
#'   0.5), while *w* = R^(delta_c)
#' @param mem_strength hyperparameter: sets the strength of the memory prior,
#'   defaults to 10 as in Bacon >= 2.5.1
#' @param model_bioturbation defaults to FALSE. If TRUE, additional uncertainty
#'   in the observed ages due to sediment mixing (bioturbation) is modelled via
#'   a latent variable process. The amount of additional uncertainty is a
#'   function of the mixing depth L, the sedimentation rate, and the number of
#'   particles (e.g. individual foraminifera) per measured date. See description
#'   for details.
#' @param n_ind the number of individual particles (e.g. foraminifera) in each
#'   sample that was dated by e.g. radiocarbon dating. This can be a single
#'   value or a vector the same length as obs_age. Only used if
#'   model_bioturbation = TRUE.
#' @param L_prior_mean mean of the gamma prior on mixing depth, defaults to 10.
#' @param L_prior_shape,L_prior_sigma shape and standard deviation of the gamma
#'   prior on the mixing depth. Set only one of these, the other will be
#'   calculated. Defaults to shape = 2. If either the shape or sigma parameter
#'   is set to zero, the mixing depth is fixed at the value of L_prior_mean,
#'   rather than being sampled.
#' @param model_displacement model additional error on observed ages that does
#' not scale with the number of individual particles in a sample, for example
#' due to incomplete mixing.
#' @param D_prior_scale scale of the half-normal prior on additional error on
#'   observed ages. The mean and standard deviation of a half-normal are equal
#'   to the scale. Units are those of the depth variable, e.g. cm.
#' @param model_hiatus optionally model a hiatus.
#' @param H_top,H_bottom limits to the location of a hiatus. By default these
#'   are set to the top and bottom data points but can be set by the user
#' @param sample_posterior if set to FALSE, hamstr skips sampling the model and
#'   returns only the data, model structure and prior parameters so that data
#'   and prior distributions can be plotted and checked prior to running a
#'   model. Defaults to TRUE
#' @param hamstr_control additional arguments to hamstr useful for debugging or
#' development. See \code{\link{hamstr_control}} for details.
#' @param stan_sampler_args additional arguments to \link[rstan]{sampling} passed as a
#' named list. e.g. list(chains = 8, iter = 4000) to run 8 MCMC chains of 4000
#' iterations instead of the default 4 chains of 2000 iterations.
#' See \code{\link{get_stan_sampler_args}} for details.
#' @return Returns an object of class "hamstr_fit", which is a list composed of
#'   the output from the stan sampler .$fit, and the list of data passed to the
#'   sampler, .$data
#' @export
#'
#' @examples
#' \dontrun{
#'
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error)
#'
#' plot(fit)
#'
#' }
hamstr_prime <- function(depth, obs_age, obs_err,
                   min_age = 1950 - as.numeric(format(Sys.Date(), "%Y")),
                   K_fine = NULL, K_factor = NULL,
                   top_depth = NULL, bottom_depth = NULL,
                   acc_mean_prior = NULL,
                   acc_shape = 1.5,
                   mem_mean = 0.5, mem_strength = 10,
                   model_bioturbation = FALSE,
                   n_ind = NULL,
                   L_prior_mean = 10,
                   L_prior_shape = 2,
                   L_prior_sigma = NULL,
                   model_displacement = FALSE,
                   D_prior_scale = 10,
                   model_hiatus = FALSE,
                   H_top = NULL, H_bottom = NULL,
                   sample_posterior = TRUE,
                   hamstr_control = list(),
                   stan_sampler_args = list()
){
  
  
  if (is.null(K_fine)== FALSE){
     if (K_fine < 2) stop("A minimum of 2 sections are required")
  } 
  
  if (is.null(K_factor)== FALSE){
    if (K_factor < 2) stop("K_factor must be 2 or greater")
    if (K_factor %% 1 > 1e-04) stop("K_factor must be an integer")
  }
  
  stan_dat <- make_stan_dat_hamstr_prime()
  
   
  used_sampler_args <- do.call(get_stan_sampler_args, stan_sampler_args)
  
  
  # set the seed here and not inside get_inits_hamstr, so that the chains are
  # different
  
  set.seed(used_sampler_args$seed)
  
  inits <- replicate(used_sampler_args$chains,
                     list(get_inits_hamstr(stan_dat)))
  
  
  
  args <- list(object = stanmodels$hamstr_prime, data = stan_dat,
               init = inits)
  
  args <- append(args, used_sampler_args)
  
  if (sample_posterior){
    fit <- do.call(rstan::sampling, args)
    
  } else if (sample_posterior == FALSE){
    fit <- NA
  }
  
  stan_dat <- append(stan_dat, used_sampler_args)
  
  info <- list(version = utils::packageVersion("hamstr"),
               time = Sys.time())
  
  out <- list(fit=fit, data=stan_dat, info = info)
  
  class(out) <- append("hamstr_fit", class(out))
  
  
  return(out)
  
}
