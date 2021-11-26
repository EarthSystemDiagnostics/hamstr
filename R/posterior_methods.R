# Extractor functions

#' Extract posterior age models or parameters
#'
#' @param object hamstr_fit object
#' @param pars The parameters to extract. If pars = "ages" (the default) a special
#' method is invoke to return the age-depth models. If pars = "pars" the memory parameters
#' and the overall mean accumulation rate are returned. Any other specification
#' invokes the default rstan as.data.frame extract method for the use specified parameters.
#' @param ... Other arguments to rstan::extract
#' @return
#'
#' @importFrom rstan extract
#' @export
#'
extract_hamstr_fit <- function(object, pars = c("ages"), ...){

  switch(pars[1],
         ages = get_posterior_ages(object),
         pars = get_posterior_parameters(object),
                as.data.frame(object$fit, pars = pars, ...)
         )

}


#' Get Posterior Parameters
#'
#' @inheritParams plot_hamstr
#'
#' @return a dataframe/tibble with posterior ages for all iterations after warmup
#' @export
#' @importFrom readr parse_number
#' @importFrom rstan extract
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error)
#'   
#' get_posterior_parameters(fit)
#' }
get_posterior_parameters <- function(hamstr_fit){

  posterior_pars <- as.data.frame(hamstr_fit$fit,
                                  pars = c("R", "w", "alpha[1]")) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(iter = 1:nrow(.))

  return(posterior_pars)

}


#' Get Posterior Age Models
#'
#' @inheritParams plot_hamstr
#'
#' @return a dataframe/tibble with posterior ages for all iterations after warmup
#' @export
#' @importFrom readr parse_number
#' @importFrom rstan extract
#' @examples
#' \dontrun{
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
#' get_posterior_ages(fit)
#' }
get_posterior_ages <- function(hamstr_fit){

  depths <- tibble::tibble(depth = hamstr_fit$data$modelled_depths,
                   idx = 1:length(hamstr_fit$data$modelled_depths))

  posterior_ages <- as.data.frame(hamstr_fit$fit, pars = "c_ages") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(iter = 1:nrow(.)) %>%
    tidyr::gather(par, age, -iter) %>%
    dplyr::mutate(idx = readr::parse_number(par),
           par = "c_ages") %>%
    dplyr::left_join(depths, .data$., by = "idx") %>%
    #dplyr::arrange(.data$par, .data$iter, .data$idx, .data$depth)
    dplyr::select(.data$iter,.data$depth, .data$age) %>%
    dplyr::arrange(.data$iter, .data$depth, .data$age)

  return(posterior_ages)
}


#' Interpolate Posterior Age Model At New Depths
#'
#' @inheritParams plot_hamstr
#' @param depth a vector of depths at which to interpolate the age models.
#' If left NULL, the depths of the age control points are used.
#'
#' @return hamstr_interpolated_ages object
#' @keywords internal
#'
#' @examples
#' \dontrun{
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
#' interpolate_age_models(fit, depth = seq(1000, 15000, by = 1000))
#' }
#'


#' Interpolate Posterior Age Model At New Depths
#'
#' @inheritParams plot_hamstr
#' @param depth a vector of depths at which to interpolate the age models.
#' If left NULL, the depths of the age control points are used.
#'
#' @return hamstr_interpolated_ages object
#' @keywords internal
#'
#' @examples
#' \dontrun{
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
#' interpolate.age.models(fit, depth = seq(1000, 15000, by = 1000))
#' }
#'
interpolate_age_models <- function(hamstr_fit, depth) {

  if (is.numeric(depth) == FALSE) {

    depth <- match.arg(depth, choices = c("modelled", "data"))

    if (depth == "modelled") {
      out <- get_posterior_ages(hamstr_fit)
      return(out)

    } else {

      if (depth == "data") {

        depth <- hamstr_fit$data$depth

      }
    }
  }

  # get posterior age models
  pst_age <- get_posterior_ages(hamstr_fit)

  # use base list, split methods, much faster than dplyr::do
  pst_age_lst <- split(pst_age, pst_age$iter)

  new_pst_age <- lapply(pst_age_lst, function(x) {
    stats::approx(x$depth, x$age, depth)$y
  })

  out <- expand.grid(depth = depth,
                     iter = 1:length(pst_age_lst))

  out$age <- unlist(new_pst_age)

  out <- dplyr::as_tibble(out)
  out <- out[, c(2, 1, 3)]

  class(out) <- append("hamstr_interpolated_ages", class(out))

  return(out)

}


#' Summarise Interpolated Posterior Age Models
#'
#' @param new_ages Object of class "hamstr_interpolated_ages"
#'
#' @return data.frame / tibble
#' @keywords internal
summarise_new_ages <- function(new_ages){

  new_ages_sum <- new_ages %>%
    dplyr::group_by(depth) %>%
    dplyr::summarise(mean = mean(age),
              #se_mean = NA,
              sd = stats::sd(age),
              `2.5%` = stats::quantile(age, probs = c(0.025), na.rm = T),
              `25%` = stats::quantile(age, probs = c(0.25), na.rm = T),
              `50%` = stats::quantile(age, probs = c(0.50), na.rm = T),
              `75%` = stats::quantile(age, probs = c(0.75), na.rm = T),
              `97.5%` = stats::quantile(age, probs = c(0.975), na.rm = T))

  return(new_ages_sum)

}


#' Summarise Posterior Age Models
#'
#' @inheritParams plot_hamstr
#' @description Extracts the summary statistics of posterior age models and attaches the depths
#' @return data.frame / tibble
#' @importFrom readr parse_number
#' @keywords internal
#' @examples
#' \dontrun{
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
#' summarise_age_models(fit)
#' }
summarise_age_models <- function(hamstr_fit){

  if (is_hamstr_interpolated_ages(hamstr_fit)){
    age_summary <- summarise_new_ages(hamstr_fit)
  } else {
    age_summary <- rstan::summary(hamstr_fit$fit, par = "c_ages")[["summary"]] %>%
      tibble::as_tibble(., rownames = "par")

    depths <- tibble::tibble(depth = hamstr_fit$data$modelled_depths,
                     idx = 1:length(hamstr_fit$data$modelled_depths))
    age_summary <- age_summary %>%
      dplyr::mutate(idx = readr::parse_number(par)) %>%
      dplyr::left_join(depths, ., by = "idx")

  }
  return(age_summary)
}

# Methods --------

#' Interpolate Age Models at Given Depths
#' @description Method for generic function predict. Returns the posterior age
#' models interpolated to new depths given in depth.
#' @param object hamstr_fit object
#' @param depth Defaults to "modelled", which returns the modelled depths.
#' "data" returns age models at the depths of the observations, or a numerical
#'  vector to specify depths. Accumulation rates are only returned at the 
#'  modelled depths.
#' @return
#'
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error
#'   )
#'
#' predict(fit, depth = seq(1, 100, by = 10))
#' predict(fit, depth = "data")
#' predict(fit)
#' predict(fit, type = "acc_rates")
#' }
#'
#' @export
#' @method predict hamstr_fit
predict.hamstr_fit <- function(object, type = c("age_models", "acc_rates"),
                               depth = c("modelled", "data"),
                               #tau = 0,
                               #kern = c("U", "G", "BH"),
                               ...){
  
  type <- match.arg(type)
  
  switch(type,
         age_models = interpolate_age_models(object, depth),
         acc_rates = get_posterior_acc_rates(object, #tau = tau, kern = kern,
                                             ...)
  )
}



#' Title
#'
#' @param object hamstr_fit object
#' @param type age models "age_models" or accumulation rates "acc_rates"
#'
#' @return
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error,
#'   K = c(10, 10))
#'
#' summary(fit)
#' summary(fit, type = c("acc_rates"))
#' }
#' @export
#' @method summary hamstr_fit
summary.hamstr_fit <- function(object, type = c("age_models", "acc_rates"),
                               #tau = 0, kern = c("U", "G", "BH"),
                               ...){

  type <- match.arg(type)

  switch(type,
         age_models = summarise_age_models(object),
         acc_rates = summarise_hamstr_acc_rates(object,
                                                #tau = tau, kern = kern,
                                                ...)
  )
}


#' Title
#'
#' @param object hamstr_interpolated_ages object
#' @return
#'
#' @export
#' @method summary hamstr_interpolated_ages
summary.hamstr_interpolated_ages <- function(object){
    summarise_new_ages(object)
  }



## Accumulation rates -----

#' Get Posterior Accumulation Rates
#'
#' @inheritParams plot_hamstr
#' @inheritParams filter_hamstr_acc_rates
#'
#' @return a dataframe/tibble with posterior ages for all iterations after warmup
#' @export
#' @importFrom readr parse_number
#' @importFrom rstan extract
#' @examples
#' \dontrun{
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
#' get_posterior_acc_rates(fit)
#' }
get_posterior_acc_rates <- function(hamstr_fit, tau = 0, kern = c("U", "G", "BH")){

  depths <- tibble::as_tibble(
    hamstr_fit$data[c("c", "c_depth_top", "c_depth_bottom")]
    ) %>%
    dplyr::mutate(depth = c_depth_top) %>%
    dplyr::rename(idx = c)

  out <- as.data.frame(hamstr_fit$fit, pars = "x") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(iter = 1:nrow(.)) %>%
    tidyr::gather(par, time_per_depth, -iter) %>%
    dplyr::mutate(idx = readr::parse_number(par),
                  par = "x") %>%
    dplyr::left_join(depths, .data$.) %>%
    dplyr::mutate(depth_per_time = 1000/time_per_depth) %>%
    dplyr::arrange(.data$par, .data$iter, .data$idx, .data$depth) %>%
    dplyr::select(iter, idx, depth, c_depth_top, c_depth_bottom, time_per_depth, depth_per_time)
  
  class(out) <- append("hamstr_acc_rates", class(out))

  #if (tau > 0){
    out <- filter_hamstr_acc_rates(out, tau = tau, kern = kern)
  #}

  return(out)

}

#' Smooth accumulation rates
#'
#' @param hamstr_acc_rates Accumulation rates from get_posterior_acc_rates
#' @param tau Scale of the smoothing kernel in depth units. If tau > 0,
#'  accumulation rates are smoothed (filtered) before summary statistics are 
#'  calculated, so that the accumulation rate at a given depth corresponds to the
#'  average rate over the depth interval tau. Default to 0. 
#' @param kern Choice of smoothin kernal. U for uniform (moving average), G for
#'  Gaussian, BH for Berger and Heath (exponential). Defaults to U
#' @return
#' @keywords internal
filter_hamstr_acc_rates <- function(hamstr_acc_rates, tau = 0, kern = c("U", "G", "BH")){
  
  kern <- match.arg(kern)
  
  # scale tau by delta_c
  tau_scl <- tau / diff(hamstr_acc_rates$depth[1:2])
  
  
  if (tau > 0){
    
    if (kern == "G"){
      
      fl <- ceiling(2*tau_scl)+1
      z <- (-fl):fl
      f <- dnorm(z, 0, tau_scl)
      f <- f / sum(f)  
      
    } else if (kern == "BH"){
      
      fl <- (3 * tau_scl+1)
      z  <- (-3*tau_scl):(3*tau_scl)
      
      # PDF of Berger and Heath
      f  <- dexp(z+tau_scl, 1/tau_scl)
      f <- f / sum(f)
    }  else if (kern == "U"){
      
      fl <- ceiling(tau_scl/2)
      z  <- (-tau_scl/2):(tau_scl/2)
      
      f <- dunif(z, (-tau_scl/2), tau_scl/2)
      f <- f / sum(f)
    }
    
    #if (tau_scl > 0) plot(z, f)
    
  }
  
  x_sum <- hamstr_acc_rates %>% 
    dplyr::group_by(iter) %>% 
    dplyr::mutate(time_per_depth = if (tau_scl == 0) {time_per_depth} else {
      # pad the vector with reversed head and tail
      # fl is half the filter length
      stats::filter(
        c(rev(head(time_per_depth, fl)), time_per_depth, rev(tail(time_per_depth, fl))),
        filter = f)[(fl+1):(dplyr::n()+fl)]
    }) %>% 
    dplyr::mutate(depth_per_time = 1000 * 1/time_per_depth, 
                  tau = tau)
  
  x_sum
  
}

#' Summarise accumulation rates
#'
#' @param hamstr_fit A hamstr_fit object
#' @inheritParams filter_hamstr_acc_rates
#' @return
#' @keywords internal
summarise_hamstr_acc_rates <- function(hamstr_fit,
                                       tau = 0,
                                       kern =  c("U", "G", "BH")
                                       ){

  kern <- match.arg(kern)
  
  x <- predict(hamstr_fit, type = "acc_rates", tau = tau, kern = kern) %>%
    tidyr::pivot_longer(cols = c(time_per_depth, depth_per_time),
                        names_to = "acc_rate_unit") 
  
  x_sum <- x %>%
    dplyr::group_by(depth, c_depth_top, c_depth_bottom, acc_rate_unit, idx, tau) %>%
    dplyr::summarise(mean = mean(value),
                     #se_mean = NA,
                     sd = stats::sd(value),
                     `2.5%` = stats::quantile(value, probs = c(0.025), na.rm = T),
                     `25%` = stats::quantile(value, probs = c(0.25), na.rm = T),
                     `50%` = stats::quantile(value, probs = c(0.50), na.rm = T),
                     `75%` = stats::quantile(value, probs = c(0.75), na.rm = T),
                     `97.5%` = stats::quantile(value, probs = c(0.975), na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(acc_rate_unit, depth)
  
  return(x_sum)

}





