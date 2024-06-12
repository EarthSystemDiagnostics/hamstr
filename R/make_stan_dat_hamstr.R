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
  
  
  l$N <- length(l$depth)
  
  # If prior for mean acc rate is provided, estimate one from the data
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
  
  # ensure the age control points are in depth order
  ord <- order(l$depth)
  
  l$depth <- l$depth[ord]
  l$obs_age <- l$obs_age[ord]
  l$obs_err <- l$obs_err[ord]
  
  
  # Displacement / depth uncertainty
  
  if (l$model_displacement == TRUE){
    
    # check parameters
    if (is.null(l$D_prior_sigma) == FALSE)
      message("D_prior_shape is being overriden by D_prior_sigma.")
    
    if (length(c(l$D_prior_sigma, l$D_prior_shape)) == 0)
      stop("One of either D_prior_sigma or D_prior_shape must be specified.
             Set either to 0 to impose a fixed depth uncertainty.")
    
    if (is.null(l$D_prior_sigma) == FALSE) {
      l$D_prior_shape <- ifelse(
        l$D_prior_sigma == 0, 0,
        gamma_sigma_shape(mean = l$D_prior_mean,
                          sigma = l$D_prior_sigma)$shape
      )}
    
    #l$D_prior_mean <- as.array(l$D_prior_mean)
    #l$D_prior_shape <- as.array(l$D_prior_shape)
    
  } else if(l$model_displacement == FALSE){
    #l$D_prior_mean <- numeric(0)
    #l$D_prior_shape <- numeric(0)
  }
  
  
  
  
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
  
  
  # Setup hierarchical structure for modelled sections
  stopifnot(l$N == length(l$obs_err), l$N == length(l$obs_age))
  
  brks <- GetBrksHalfOffset(K_fine = l$K_fine, K_factor = l$K_factor)
  alpha_idx <- GetIndices(brks = brks)
  
  l$K_tot <- sum(alpha_idx$nK)
  l$K_fine <- utils::tail(alpha_idx$nK, 1)
  l$c <- 1:l$K_fine
  
  
  # Transformed arguments
  l$mem_alpha = l$mem_strength * l$mem_mean
  l$mem_beta = l$mem_strength * (1-l$mem_mean)
  
  l$mem_mean = l$mem_mean
  l$mem_strength = l$mem_strength
  
  l$delta_c = depth_range / l$K_fine
  
  # depth at top and bottom of each modelled section
  l$c_depth_bottom = l$delta_c * l$c + l$top_depth
  l$c_depth_top = c(l$top_depth, l$c_depth_bottom[1:(l$K_fine-1)])
  
  l$modelled_depths <- c(l$c_depth_top[1], l$c_depth_bottom)
  
  # Index for which sections the age control points are in
  l$which_c = sapply(l$depth, function(d) which.max((l$c_depth_bottom < d) * (l$c_depth_bottom - d) ))
  
  l <- append(l, alpha_idx)
  
  l$n_lvls <- length(l$nK) -1
  l$scale_shape = as.numeric(l$scale_shape)
  l$model_bioturbation = as.numeric(l$model_bioturbation)
  l$model_displacement = as.numeric(l$model_displacement)
  
  # Model hiatus? and set upper/lower limits for position of hiatus 
  l$model_hiatus = as.numeric(l$model_hiatus)
  if (is.null(l$H_top)) l$H_top = l$top_depth
  if (is.null(l$H_bottom)) l$H_bottom = l$bottom_depth
  
  if (is.null(l$H_max)) l$H_max = ceiling(diff(range(l$obs_age)))
  
  # set scale of smoothing of acc_rates for bioturbation calculation
  l$smooth_s = as.numeric(l$smooth_s)
  
  if (l$smooth_s == 1){
    l$smooth_i <- get_smooth_i(l, l$L_prior_mean)
    l$I <- nrow(l$smooth_i)
  } else {
    l$smooth_i <- rbind(rep(1, l$N))
    l$I <- 1
  }
  
  return(l)
}

# Internal functions ---------

# new K_fine and parent calculation -----

#' Get Weights for Parent Sections
#'
#' @param a Parent breaks
#' @param b Child breaks
#' @keywords internal
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


#' Get Indices Structure For Hamstr Model from a Set of Breaks
#'
#' @param nK List of number of breaks in each level
#' @param brks List of breakpoints in each level
#' @keywords internal
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


#' Get Overlapping Breaks Structure
#'
#' @inheritParams hamstr
#' @keywords internal
GetBrksHalfOffset <- function(K_fine, K_factor){
  
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
  
  brks <- c(brks, list(c(newbrks[1], utils::tail(newbrks, 1))))
  brks <- rev(brks)
  
  # fixes behaviour when K_factor is very large and allows for a flat structure
  if (length(brks[[1]]) == length(brks[[2]])){
    if (all(brks[[1]] == brks[[2]])){
      brks <- brks[2:length(brks)]
    }
  }
  return(brks)
}





#' Get Default K_factor
#'
#' @param K_fine The number of sections at the highest resolution
#' @keywords internal
get_K_factor <- function(K_fine){
  
  bar <- function(x, y){
    abs(y - x^x)
  }

  ceiling(stats::optimize(bar, c(1, (10 + log10(K_fine))), y = K_fine)$minimum)

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
#' \dontrun{
#' gamma_sigma_shape(mean = 10, sigma = 2)
#' }
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
    # from mean and shape
    rate <- shape / mean
    sigma <- sqrt(shape/(rate^2))
  } else if (is.null(shape)){
    # from mean and sigma
    rate <- mean / sigma^2
    shape <- mean^2 / sigma^2
  }

  if (is.null(mode)){
   mode <- (shape - 1) / rate
   if (mode < 0) mode <- 0
  }

  return(list(mean = mean, mode = mode, rate = rate, shape = shape, sigma = sigma))
}




#' Get Indices For Smoothing Accumulation Rates When Estimating L
#'
#' @param d 
#' @param w 
#'
#' @keywords internal
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
#' @return a list of boundaries for the modelled sections
#' @keywords internal
hierarchical_depths <- function(stan_dat) {

  # make backward compatible with branching hamstr
  if (is.null(stan_dat$brks)) {
    get_brks <- function(stan_dat) {
      d_range <- diff(range(stan_dat$modelled_depths))
      min_d <- min(stan_dat$modelled_depths)

      delta_d <- d_range / stan_dat$nK

      lapply(stan_dat$nK[-1], function(x) {
        delta_d <- d_range / x

        c(min_d, delta_d * 1:x + min_d)
      })
    }

    stan_dat$brks <- get_brks(stan_dat)
    return(stan_dat$brks)
  }

  lapply(stan_dat$brks, function(x) {
      rng <- stan_dat$bottom_depth - stan_dat$top_depth

      tcks <- x * rng + stan_dat$top_depth

      tcks[tcks <= stan_dat$bottom_depth &
             tcks >= stan_dat$top_depth]

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

