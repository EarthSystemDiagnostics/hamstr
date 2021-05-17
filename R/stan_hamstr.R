#' hamstr
#' @param depth Depths of observed ages (age control points)
#' @param obs_age Observed age at each depth (age control points)
#' @param obs_err Error associated with each observed age (1 standard error)
#' @param top_depth,bottom_depth The top and bottom depths of the desired
#'   age-depth model. Must encompass the range of the data. Defaults to the
#'   shallowest and deepest data points unless \code{pad_top_bottom = TRUE}
#' @param pad_top_bottom logical, pad the length of the age-depth model by 5% on
#'   each end
#' @param K A vector with the number of new sub-divisions at each hierarchical
#'   level, defaults to c(10, 10), resulting in two levels with 10 and 100
#'   sections respectively.
#' @param nu Degrees of freedom for the Student-t distributed error model.
#'   Defaults to 6, which is equivalent to the default parameterisation of
#'   t.a=3, t.b=4 in Bacon 2.2. Set to a high number to approximate a Gaussian
#'   error model, (nu = 100 should do it).
#' @param acc_mean_prior Hyperparameter for the prior on the overall mean
#'   accumulation rate for the record. Units are obs_age / depth. E.g. if depth
#'   is in cm and age in years then the accumulation rate is in years/cm. The
#'   overall mean accumulation rate is given a weak half-normal prior with mean
#'   = 0, SD = 10 * acc_mean_prior. If left blank, acc_mean_prior is set to the
#'   mean accumulation rate estimated by fitting a robust linear model using
#'   \link[MASS]{rlm}.
#' @param acc_shape Hyperparameter for the shape of the priors on accumulation rates.
#'  Defaults to 1.5 - as for Bacon 2.2.
#' @param scale_shape Scale the shape parameter according to the number of hierarchical
#' levels, to control the total variance of the alpha innovations. This defaults
#' to TRUE as of Hamstr verion 0.5.
#' @param mem_mean Hyperparameter: a parameter of the Beta prior distribution on
#'   "memory", i.e. the autocorrelation parameter in the underlying AR1 model.
#'   The prior on the correlation between layers is scaled according to the
#'   thickness of the sediment sections in the highest resolution hierarchical
#'   layer, *delta_c*, which is determined by the total length age-models and
#'   the parameter vector *K*. mem_mean sets the mean value for *R* (defaults to
#'   0.5), while *w* = R^(delta_c)
#' @param mem_strength Hyperparameter: sets the strength of the memory prior,
#'   defaults to 10 as in Bacon >= 2.5.1
#' @param scale_R logical: Scale AR1 coefficient by delta_c (as in Bacon) or
#' not. Defaults to TRUE.
#' @param inflate_errors logical: If set to TRUE, observation errors are
#'   inflated so that data are consistent with a "Bacon-style" monotonic
#'   age-depth model. This is an experimental feature under active development.
#'   Defaults to FALSE.
#' @param infl_sigma_sd Hyperparameter: sets the standard deviation of the
#' half-normal prior on the mean of the additional error terms. Defaults to 10
#' times the mean observation error in obs_err.
#' @param infl_shape_shape,infl_shape_mean Hyperparameters: parametrises the
#' gamma prior on the shape of the distribution of the additional error terms.
#'  Default to 1, 1.
#' @param ... additional arguments to \link[rstan]{sampling}
#' @inheritParams rstan::sampling
#'
#' @return Returns a list composed of the output from the Stan sampler .$fit,
#'   and the list of data passed to the sampler, .$data
#' @export
#'
#' @examples
#' \dontrun{
#'
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error,
#'   K = c(10, 10), nu = 6,
#'   acc_mean_prior = 20,
#'   mem_mean = 0.5, mem_strength = 10,
#'   inflate_errors = 0,
#'   iter = 2000, chains = 3)
#'
#' print(fit$fit, par = c("record_acc_mean"))
#'
#' plot_hamstr(fit, 100, plot_diagnostics = TRUE)
#'
#' }
hamstr <- function(depth, obs_age, obs_err,
                   K = NULL,
                   top_depth = NULL, bottom_depth = NULL,
                   pad_top_bottom = FALSE,
                   acc_mean_prior = NULL,
                   acc_shape = 1.5,
                   scale_shape = TRUE,
                   mem_mean = 0.5, mem_strength = 10,
                   scale_R = TRUE,
                   nu = 6,
                   inflate_errors = FALSE,
                   infl_sigma_sd = NULL,
                   infl_shape_shape = 1, infl_shape_mean = 1,
                   iter = 2000, chains = 3, ...){


  stan_dat <- make_stan_dat_hamstr(depth = depth,
                                   obs_age = obs_age, obs_err = obs_err,
                                   K=K,
                                   top_depth = top_depth,
                                   bottom_depth = bottom_depth,
                                   pad_top_bottom = pad_top_bottom,
                                   acc_mean_prior = acc_mean_prior,
                                   acc_shape = acc_shape,
                                   scale_shape = scale_shape,
                                   mem_mean=mem_mean, mem_strength=mem_strength,
                                   scale_R = as.numeric(scale_R),
                                   nu=nu,
                                   inflate_errors = as.numeric(inflate_errors),
                                   infl_sigma_sd = infl_sigma_sd,
                                   infl_shape_shape = infl_shape_shape,
                                   infl_shape_mean = infl_shape_mean)

  inits <- get_inits_hamstr(stan_dat)

  inits <- rep(list(inits), chains)

  fit <- rstan::sampling(stanmodels$hamstr,
                         data = stan_dat, init = inits, iter = iter, chains = chains,
                         verbose = FALSE, ...)

  out <- list(fit=fit, data=stan_dat)

  class(out) <- append("hamstr_fit", class(out))


  return(out)

}


# Methods and classes -------

is_hamstr_fit <- function(x) inherits(x, "hamstr_fit")
is_hamstr_interpolated_ages <- function(x) inherits(x, "hamstr_interpolated_ages")



