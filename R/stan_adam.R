#' adam
#' @param depth Vector of depths from an observed age-depth profile
#' @param obs_age Observed age at each depth
#' @param obs_err Error associated with each observed age (1 standard error)
#' @param top_depth,bottom_depth the top and bottom depths of the age-depth
#'   model. Must encompass the range of the data. Default to the shallowest and
#'   deepest data points unless \code{pad_top_bottom = TRUE}
#' @param pad_top_bottom logical, pad the length of the age-depth model by 5% on
#'   each end
#' @param K Number of sub-sections at each level
#' @param nu Degrees of freedom for the Student-t distributed error model.
#'   Defaults to 6, which is equivalent to the default parameterisation of
#'   t.a=3, t.b=4 in Bacon 2.2. Set to a high number to approximate a Gaussian
#'   error model, (nu = 100 should do it).
#' @param record_prior_acc_mean_mean,record_prior_acc_mean_shape hyperparameters
#'   for the hyperprior on the overall mean accumulation rate for the record.
#'   Units are depth / obs_age. E.g. if depth is in cm and age in kyr then the
#'   accumulation rate is in cm/kyr. ...mean_shape sets the shape of the
#'   hyperprior for the mean. Defaults to 1.5.
#' @param record_prior_acc_shape_mean,record_prior_acc_shape_shape
#'   hyperparameters for the hyperprior on the shape of the distribution of K1
#'   coarse level accumulation rates. Defaults to 1.5, 1.5.
#' @param section_acc_shape hyperparameter for the prior on the shape of the K
#'   fine level accumulation rates. Defaults to 1.5 - as for Bacon 2.2.
#' @param mem_mean Hyperparameter: a parameter of the Beta prior distribution on
#'   "memory", i.e. the autocorrelation parameter in the underlying AR1 model.
#'   The prior on the correlation between layers is scaled according to the
#'   thickness of the layers *delta_c*, which is determined by the total length
#'   of the profile and the parameter *K*. mem_mean sets the mean value for *R*
#'   (defaults to 0.7), while *w* = R^(delta_c)
#' @param mem_strength Hyperparameter: sets the strength of the memory prior,
#'   defaults to 4 as in Bacon 2.2
#' @param inflate_errors logical. If set to TRUE, observation errors are
#'   inflated so that data are consistent with a Bacon type monotonic age-depth
#'   model
#' @inheritParams rstan::sampling
#'
#' @return Returns a list composed of the output from the Stan sampler .$fit,
#'   and the list of data passed to the sampler
#' @export
#'
#' @examples
#' dontrun{
#' 
#' fit <- adam(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error,
#'   K = c(10, 10), nu = 6,
#'   record_prior_acc_mean_mean = 20,
#'   mem_mean = 0.7, mem_strength = 4,
#'   inflate_errors = 0,
#'   iter = 2000, chains = 3)
#'
#' print(fit$fit, par = c("record_acc_mean"))
#'
#' plot_stan_bacon(fit, 100, plot_priors = TRUE)
#' 
#' }
adam <- function(depth, obs_age, obs_err,
                   K = c(10, 10), nu = 6,
                   top_depth = NULL, bottom_depth = NULL,
                   pad_top_bottom = FALSE,
                   record_prior_acc_mean_mean = NULL,
                   record_prior_acc_mean_shape = 1.5,
                   record_prior_acc_shape_mean = 1.5,
                   record_prior_acc_shape_shape = 1.5,
                   section_acc_shape = 1.5,
                   mem_mean = 0.7, mem_strength = 4,
                   inflate_errors = 0,
                   iter = 2000, chains = 3, ...){
  
  
  stan_dat <- make_stan_dat_adam(depth = depth, obs_age = obs_age, obs_err = obs_err,
                                   K=K, nu=nu,
                                   top_depth = top_depth, bottom_depth = bottom_depth,
                                   pad_top_bottom = pad_top_bottom,
                                   record_prior_acc_mean_mean = record_prior_acc_mean_mean,
                                   record_prior_acc_mean_shape = record_prior_acc_mean_shape,
                                   record_prior_acc_shape_mean = record_prior_acc_shape_mean,
                                   record_prior_acc_shape_shape = record_prior_acc_shape_shape,
                                   section_acc_shape = section_acc_shape,
                                   mem_mean=mem_mean, mem_strength=mem_strength,
                                   inflate_errors = inflate_errors)
  
  inits <- get_inits_adam(stan_dat)
  
  inits <- rep(list(inits), chains)
  
  fit <- rstan::sampling(stanmodels$adam,
                         data = stan_dat, init = inits, iter = iter, chains = chains,
                         verbose = FALSE, ...)
  
  out <- list(fit=fit, data=stan_dat)
  
  class(out) <- append(class(out), "adam_fit")
  
  
  return(out)
  
}
