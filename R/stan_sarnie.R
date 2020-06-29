#' sarnie
#' @param depth Vector of depths from an observed age-depth profile
#' @param obs_age Observed age at each depth
#' @param obs_err Error associated with each observed age (1 standard error)
#' @param K Number of sections into which the profile will be divided, must be a multiple of K1
#' @param K1 Number of coarse level sections
#' @param nu Degrees of freedom for the Student-t distributed error model.
#'   Defaults to 6, which is equivalent to the default parameterisation of
#'   t.a=3, t.b=4 in Bacon 2.2. Set to a high number to approximate a Gaussian
#'   error model, (nu = 100 should do it).
#' @param record_prior_acc_mean_mean hyperprior parameter for the prior on the overall mean acc.rate for 
#' the record, set to the mean acc.rate for the record from a linear model or similar
#' @param record_prior_acc_mean_shape shape for the prior on the overall mean, leave at 1.5 for now 
#' acc.rate for the record
#' @param record_prior_acc_shape_mean hyperprior parameters for the prior on 
#' the shape of the innovations distribution, leave at 1.5 for now
#' @param record_prior_acc_shape_shape hyperprior parameters for the prior on 
#' the shape of the innovations distribution, leave at 1.5 for now
#' @param acc_mean The mean sediment accumulation rate for the Gamma prior on
#'   sedimentation rate innovations
#' @param acc_alpha Sets the alpha (shape) parameter for the Gamma prior on
#'   sedimentation rate innovations. The variance of the innovations acc_var =
#'   acc_mean^2^ / acc_alpha. The default is acc_alpha = 1.5 which is equivalent
#'   to the default shape = 1.5 parameter in Bacon 2.2
#' @param mem_mean Hyper-parameter: a parameter of the Beta prior distribution
#'   on "memory", i.e. the autocorrelation parameter in the underlying AR1
#'   model. The prior on the correlation between layers is scaled according to
#'   the thickness of the layers *delta_c*, which is determined by the total
#'   length of the profile and the parameter *K*. mem_mean sets the mean value
#'   for *R* (defaults to 0.7), while *w* = R^(delta_c)
#' @param mem_strength Hyper-parameter: sets the strength of the memory prior,
#'   defaults to 4 as in Bacon 2.2
#' @param inflate_errors logical, 1 or 0. If set to 1, observation errors are 
#' inflated so that data are consistent with a Bacon type monotonic age-depth model
#' @inheritParams rstan::sampling
#'
#' @return Returns a list composed of the output from the Stan sampler .$fit,
#'   and the list of data passed to the sampler
#' @export
#'
#' @examples
#'
#' fit <- sarnie(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error,
#'   K = 100, K1 = 10, nu = 6,
#'   record_prior_acc_mean_mean = 20,
#'   mem_mean = 0.7, mem_strength = 4,
#'   inflate_errors = 0,
#'   iter = 2000, chains = 3)
#'
#' print(fit$fit, par = c("record_acc_mean"))
#'
#' plot_stan_bacon(fit, 100, plot_priors = FALSE)
#'
sarnie <- function(depth, obs_age, obs_err,
                       K1 = 10, K = 100, nu = 6,
                       record_prior_acc_mean_mean = 20,
                       record_prior_acc_mean_shape = 1.5,
                       record_prior_acc_shape_mean = 1.5,
                       record_prior_acc_shape_shape = 1.5,
                       mem_mean = 0.7, mem_strength = 4,
                       inflate_errors = 0,
                       iter = 2000, chains = 3, ...){
  
  stan_dat <- make_stan_dat_sarnie(depth = depth, obs_age = obs_age, obs_err = obs_err,
                                   K1 = K1, K=K, nu=nu,
                                   record_prior_acc_mean_mean = record_prior_acc_mean_mean,
                                   record_prior_acc_mean_shape = record_prior_acc_mean_shape,
                                   record_prior_acc_shape_mean = record_prior_acc_shape_mean,
                                   record_prior_acc_shape_shape = record_prior_acc_shape_shape,
                                   mem_mean=mem_mean, mem_strength=mem_strength,
                                   inflate_errors = inflate_errors)
  
  inits <- get_inits_sarnie(stan_dat)
  
  inits <- rep(list(inits), chains)
  
  fit <- rstan::sampling(stanmodels$sarnie,
                         data = stan_dat, init = inits, iter = iter, chains = chains,
                         verbose = FALSE, ...)

  return(list(fit=fit, data=stan_dat))
  
}
