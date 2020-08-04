# Make index creating functions for K levels
alpha_indices <- function(K){

  K <- c(1, K)

  nK <- cumprod(K)

  K <- c(0, K)
  nK <- c(0, nK)

  alpha_idx <- 1:sum(nK)

  nLevels <- length(K)-1

  lvl <- unlist(lapply(seq_along(nK[-1]), function(i) rep(i, times = nK[i+1])))

  parent <- c(rep(0, K[2]), unlist(lapply(alpha_idx[1:sum(nK[1:nLevels])], function(i) rep(i, K[lvl[i]+2]))))

  list(alpha_idx=alpha_idx, lvl=lvl, parent=parent, nK = nK[-1])

}

#level_indices(nLevels = 3, K_per_level = 3)

#alpha_indices(K = c(3, 3, 3, 3, 3, 3))


#' Title
#'
#' @inheritParams adam
#'
#' @return a list of data, and parameters to be passed as data to the Stan sampler
#' @export
#'
#' @examples
#' @return
#' @export
#'
#' @examples
#' make_stan_dat_adam(depth = MSB2K$depth,
#'               obs_age = MSB2K$age,
#'               obs_err = MSB2K$error,
#'               #K1 = 10, K = 100,
#'               nu = 6)
make_stan_dat_adam <- function(depth, obs_age, obs_err,
                               top_depth = NULL, bottom_depth = NULL,
                                 pad_top_bottom = FALSE,
                                 K = c(10, 10), nu = 6,
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

  return(l)
}

get_inits_adam <- function(dat){
  l <- list(
    R = runif(1, 0.1, 0.9),
    alpha = (abs(rnorm(dat$K_tot, dat$record_prior_acc_mean_mean, dat$record_prior_acc_mean_mean/3))),
    record_acc_mean = (abs(rnorm(1, dat$record_prior_acc_mean_mean, dat$record_prior_acc_mean_mean/3))),

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

