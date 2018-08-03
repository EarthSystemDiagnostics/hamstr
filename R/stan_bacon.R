#' stan_bacon
#' @md
#' @useDynLib baconr, .registration = TRUE
#' @exportPattern("^[[:alpha:]]+")
#' @param depth Vector of depths from an observed age-depth profile
#' @param obs_age Observed age at each depth
#' @param obs_err Error associated with each observed age (1 standard error)
#' @param K Number of sections into which the profile will be divided
#' @param nu Degrees of freedom for the Student-t distributed error model. 
#'   Defaults to 6, which is equivalent to the default parameterisation of 
#'   t.a=3, t.b=4 in Bacon 2.2. Set to a high number to approximate a Gaussian 
#'   error model, (nu = 100 should do it).
#' @param acc_mean the mean sediment accumulation rate for the 
#'   Gamma prior on sedimentation rate innovations
#' @param acc_alpha sets the alpha (shape) parameter for the Gamma prior on 
#'   sedimentation rate innovations. The variance of the innovations 
#'   acc_var = acc_mean^2^ / acc_alpha. The default is acc_alpha = 1.5 
#'   which is equivalent to the default shape = 1.5 parameter in Bacon 2.2
#' @param mem_mean Hyper-parameter: a parameter of the Beta prior distribution 
#'   on "memory", i.e. the autocorrelation parameter in the underlying AR1 
#'   model. The prior on the correlation between layers is scaled according to
#'   the thickness of the layers *delta_c*, which is determined by the total
#'   length of the profile and the parameter *K*. mem_mean sets the mean value
#'   for *R* (defaults to 0.7), while *w* = R^(delta_c) 
#' @param mem_strength Hyper-parameter: sets the strength of the memory prior, 
#'   defaults to 4 as in Bacon 2.2
#' @inheritParams rstan::sampling
#'   
#' @return Returns a list composed of the output from the Stan sampler .$fit,
#'   and the list of data passed to the sampler
#' @export
#' 
#' @examples
stan_bacon <- function(depth, obs_age, obs_err,
                       hiatus_depth = NULL, hiatus_length = NULL,
                       hiatus_shape = 1, hiatus_interval = 0.1,
                       K = 10, nu = 6,
                       acc_mean, acc_alpha = 1.5,
                       mem_mean = 0.7, mem_strength = 4,
                       iter = 2000, chains = 4){
  
  stan_dat <- make_stan_dat(depth, obs_age, obs_err,
                            hiatus_depth, hiatus_length,
                            hiatus_shape, hiatus_interval,
                            K, nu,
                            acc_mean, acc_alpha,
                            mem_mean, mem_strength)
  
  if (is.null(hiatus_depth)){
    fit <- rstan::sampling(stanmodels$bacon, 
                           data = stan_dat, iter = iter, chains = chains,
                           verbose = FALSE)  
  }else if(is.null(hiatus_depth) == FALSE){
    fit <- rstan::sampling(stanmodels$bacon_hiatuses, 
                           data = stan_dat, iter = iter, chains = chains,
                           verbose = FALSE) 
  }

  return(list(fit=fit, data=stan_dat))
  
}
