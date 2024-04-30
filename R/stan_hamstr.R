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
#' @param K deprecated, use K_fine and K_factor instead.
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
hamstr <- function(depth, obs_age, obs_err,
                   min_age = 1950 - as.numeric(format(Sys.Date(), "%Y")),
                   K_fine = NULL, K_factor = NULL, K,
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
  
  if (!missing(K)) {
    warning("argument K is deprecated; K_fine has been calculated from K but please use K_fine instead.", 
            call. = FALSE)
    K_fine <- utils::tail(cumprod(K), 1)
  }
  
  if (is.null(K_fine)== FALSE){
    if (K_fine < 2) stop("A minimum of 2 sections are required")
  } 
  
  if (is.null(K_factor)== FALSE){
    if (K_factor < 2) stop("K_factor must be 2 or greater")
    if (K_factor %% 1 > 1e-04) stop("K_factor must be an integer")
  }
  
  stan_dat <- make_stan_dat_hamstr()
  
  
  used_sampler_args <- do.call(get_stan_sampler_args, stan_sampler_args)
  
  
  # set the seed here and not inside get_inits_hamstr, so that the chains are
  # different
  
  set.seed(used_sampler_args$seed)
  
  inits <- replicate(used_sampler_args$chains,
                     list(get_inits_hamstr(stan_dat)))
  
  
  
  args <- list(object = stanmodels$hamstr, data = stan_dat,
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

#' Return a List of hamstr Control Arguments
#'
#' @description Returns a list arguments to pass to hamstr_control which may be
#' useful for debugging or development
#'
#' @param scale_shape scale the shape parameter according to the number of
#'   hierarchical levels, to control the total variance of the alpha
#'   innovations. This defaults to TRUE as of hamstr verion 0.5.
#' @param scale_R logical: Scale AR1 coefficient by delta_c (as in Bacon) or
#'   not. Defaults to TRUE.
#' @param nu degrees of freedom for the Student-t distributed error model.
#'   Defaults to 6, which is equivalent to the default parametrisation of
#'   t.a=3, t.b=4 in Bacon 2.2. Set to a high number to approximate a Gaussian
#'   error model, (nu = 100 should do it).
#' @param smooth_s smooth the sedimentation rate used to calculate additional
#' error from bioturbation by taking a mean across nearby sections
#' @param inflate_errors logical: If set to TRUE, observation errors are
#'   inflated so that data are consistent with a gamma AR1 age-depth model.
#'   This is an experimental feature under active development. Defaults to FALSE.
#' @param infl_sigma_sd hyperparameter: sets the standard deviation of the
#'   half-normal prior on the mean of the additional error terms. Defaults to 10
#'   times the mean observation error in obs_err.
#' @param infl_shape_shape,infl_shape_mean hyperparameters: parametrises the
#'   gamma prior on the shape of the distribution of the additional error terms.
#'   Default to 1, 1.
#' @examples
#' hamstr_control()
#' @return a named list
#' @export
#'
hamstr_control <- function(scale_R = TRUE,
                           nu = 6,
                           scale_shape = TRUE,
                           smooth_s = FALSE,
                           inflate_errors = FALSE,
                           infl_sigma_sd = NULL,
                           infl_shape_shape = 1,
                           infl_shape_mean = 1
                           ){

  l <- c(as.list(environment()))

  return(l)

  }



#' Default Parameters for Sampling Hamstr Models with Stan
#' @description Returns a list of parameters for the Stan sampler
#' @inheritParams rstan::sampling
#'
#' @return a named list
#' @export
#'
#' @examples
#' get_stan_sampler_args()
get_stan_sampler_args <- function(chains = 4,
                             cores = chains,
                             iter = 2000,
                             warmup = floor(iter / 2),
                             thin = 1,
                             seed = sample.int(.Machine$integer.max, 1),
                             check_data = TRUE,
                             sample_file = NULL,
                             diagnostic_file = NULL,
                             verbose = FALSE,
                             algorithm = c("NUTS", "HMC", "Fixed_param"),
                             control = NULL,
                             include = TRUE,
                             open_progress = interactive() &&
                               !isatty(stdout()) &&
                               !identical(Sys.getenv("RSTUDIO"), "1"),
                             show_messages = TRUE,
                             ...) {
  l <- c(as.list(environment()), list(...))

  return(l)
}


# Methods and classes -------

is_hamstr_fit <- function(x) inherits(x, "hamstr_fit")
is_hamstr_interpolated_ages <- function(x) inherits(x, "hamstr_interpolated_ages")

