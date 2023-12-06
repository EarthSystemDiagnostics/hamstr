# Extractor functions -------

#' Get Posterior Age Models
#' @inheritParams plot_hamstr
#' @return a dataframe/tibble with posterior ages for all iterations after warmup

#' @importFrom rstan extract
#' @keywords internal
#' @examples 
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error)
#'
#' get_posterior_ages(fit)
#' }
get_posterior_ages <- function(hamstr_fit){
  depths <- tibble::tibble(depth = hamstr_fit$data$modelled_depths,
                   idx = 1:length(hamstr_fit$data$modelled_depths))

  posterior_ages <- as.data.frame(hamstr_fit$fit, pars = "c_ages") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(iter = 1:dplyr::n()) %>%
    tidyr::pivot_longer(tidyr::starts_with("c_ages"),
                        names_to = "par", values_to = "age") %>% 
    dplyr::mutate(idx = get_par_idx(.data$par),
           par = "c_ages") %>%
    dplyr::left_join(depths, .data$., by = "idx") %>%
    dplyr::select("iter", "depth", "age") %>%
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
#'   obs_err = MSB2K$error)
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


#' Summarise to Quantiles and Moments
#'
#' @param dat a dataframe or tibble
#' @param var the variable to summarise (unquoted)
#' @param probs the quantiles at which to summarise
#'
#' @return a tibble containing the summarised posterior
#' @keywords internal
summarise_q <- function(dat,
                        var,
                        probs = c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975)){
  dat %>% 
    dplyr::reframe(mean = mean({{ var }}, na.rm = TRUE),
              sd = stats::sd({{ var }}, na.rm = TRUE),
              x = stats::quantile({{ var }}, probs, na.rm = TRUE),
              q = paste0(round(100*probs, 1), "%")) %>% 
    tidyr::pivot_wider(names_from = "q", values_from = "x") %>% 
    dplyr::as_tibble()
}



#' Summarise Interpolated Posterior Age Models
#'
#' @param new_ages Object of class "hamstr_interpolated_ages"
#'
#' @return data.frame / tibble
#' @keywords internal
summarise_new_ages <- function(new_ages, probs = c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975)){

  new_ages_sum <- new_ages %>%
    dplyr::group_by(.data$depth) %>%
    summarise_q(var = .data$age, probs = probs)

  return(new_ages_sum)

}


#' Summarise Posterior Age Models
#'
#' @inheritParams plot_hamstr
#' @description Extracts the summary statistics of posterior age models and attaches the depths
#' @return data.frame / tibble

#' @keywords internal
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error)
#'
#' summarise_age_models(fit)
#' }  
summarise_age_models <- function(hamstr_fit, probs = c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975)){

  if (is_hamstr_interpolated_ages(hamstr_fit)){
    age_summary <- summarise_new_ages(hamstr_fit, probs = probs)
  } else {
    age_summary <- rstan::summary(hamstr_fit$fit, par = "c_ages", probs = probs)[["summary"]] %>%
      tibble::as_tibble(rownames = "par")

    depths <- tibble::tibble(depth = hamstr_fit$data$modelled_depths,
                     idx = 1:length(hamstr_fit$data$modelled_depths))
    age_summary <- age_summary %>%
      dplyr::mutate(idx = get_par_idx(.data$par)) %>%
      dplyr::left_join(depths, .data$., by = "idx") %>% 
      dplyr::select( "depth", "idx", tidyr::everything())

  }
  return(age_summary)
}

#' Summarise Hamstr Parameters
#'
#' @param object A hamstr_fit object
#' @param par Character vector of parameters to include
#'
#' @return a tibble of summarised posterior of hamstr parameters
#' @keywords internal
summarise_hamstr_parameters <- function(object,
                                        pars = c("alpha[1]", "R", "w", "L", "D",
                                                 "H_depth", "H_length"),
                                        probs = c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975)) {
  rstan::summary(object$fit,
                 pars = pars, probs = probs)$summary %>%
    dplyr::as_tibble(rownames = "Parameter") %>%
    dplyr::select(-"se_mean")
}


# Methods --------

#' Interpolate Age Models at Arbitrary Depths
#' @description predict method for class "hamstr_fit". Returns the posterior age
#' models interpolated to arbitrary depths.
#' @param object hamstr_fit object
#' @param type age models "age_models" or accumulation rates "acc_rates"
#' @param depth defaults to "modelled", which returns the modelled depths.
#' "data" returns age models at the depths of the observations, or a numerical
#'  vector to specify depths. Accumulation rates are only returned at the 
#'  modelled depths.
#' @param ... additional arguments to hamstr predict methods
#' @return a tibble of hamstr age-depth model realisations
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
predict.hamstr_fit <- function(object,
                               type = c("age_models", "acc_rates"),
                               depth = c("modelled", "data"),
                               ...){
  
  type <- match.arg(type)
  
  switch(type,
         age_models = interpolate_age_models(object, depth),
         acc_rates = get_posterior_acc_rates(object, #tau = tau, kern = kern,
                                             ...)
  )
}



#' Summarise hamstr Models
#' @description summary method for class "hamstr_fit"
#' @param object hamstr_fit object
#' @param type age models "age_models" or accumulation rates "acc_rates"
#' @param ... additional arguments to hamstr summary methods
#' @return a tibble
#' @examples 
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
summary.hamstr_fit <- function(object, type = c("age_models", "acc_rates", "pars"),
                               #probs = c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975), 
                               #tau = 0, kern = c("U", "G", "BH"),
                               ...){

  type <- match.arg(type)

  switch(type,
         age_models = summarise_age_models(object, ...),
         acc_rates = summarise_hamstr_acc_rates(object,
                                                #probs = probs,
                                                #tau = tau, kern = kern,
                                                ...),
         pars = summarise_hamstr_parameters(object, 
                                            #probs = probs,
                                            ...)
         
  )
}



#' Summarise Age Models That Have Been Interpolated
#' @description summary method for class "hamstr_interpolated_ages"
#' @param object hamstr_interpolated_ages object
#' @param ... additional arguments to hamstr summary methods
#' @return a tibble containing the summarised posterior

#'
#' @export
#' @method summary hamstr_interpolated_ages
summary.hamstr_interpolated_ages <- function(object, ...){
    summarise_new_ages(object, ...)
  }



## Accumulation rates -----

#' Get Posterior Accumulation Rates
#'
#' @inheritParams plot_hamstr
#' @inheritParams filter_hamstr_acc_rates
#'
#' @return a dataframe/tibble with posterior ages for all iterations after warmup
#' @keywords internal

#' @importFrom rstan extract
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error)
#'
#' get_posterior_acc_rates(fit)
#' }
get_posterior_acc_rates <- function(hamstr_fit, tau = 0, kern = c("U", "G", "BH")){

  depths <- tibble::as_tibble(
    hamstr_fit$data[c("c", "c_depth_top", "c_depth_bottom")]
    ) %>%
    dplyr::mutate(depth = .data$c_depth_top) %>%
    dplyr::rename(idx = "c")

  out <- as.data.frame(hamstr_fit$fit, pars = "x") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(iter = 1:dplyr::n()) %>%
    tidyr::pivot_longer(cols = -tidyr::contains("iter"),
                        names_to = "par", values_to = "time_per_depth") %>% 
    dplyr::mutate(idx = get_par_idx(.data$par),
                  par = "x") %>%
    dplyr::right_join(depths) %>%
    dplyr::mutate(depth_per_time = 1000/.data$time_per_depth) %>%
    dplyr::arrange(.data$par, .data$iter, .data$idx, .data$depth) %>%
    dplyr::select("iter", "idx",  "depth",
                  "c_depth_top", "c_depth_bottom", 
                  "time_per_depth", "depth_per_time")
  
  class(out) <- append("hamstr_acc_rates", class(out))
  
  out <- filter_hamstr_acc_rates(out, tau = tau, kern = kern)
  
  return(out)

}

#' Smooth accumulation rates
#'
#' @param hamstr_acc_rates accumulation rates from get_posterior_acc_rates
#' @param tau scale of the smoothing kernel in depth units. If tau > 0,
#'  accumulation rates are smoothed (filtered) before summary statistics are 
#'  calculated, so that the accumulation rate at a given depth corresponds to the
#'  average rate over the depth interval tau. Default to 0. 
#' @param kern choice of smoothin kernal. U for uniform (moving average), G for
#'  Gaussian, BH for Berger and Heath (exponential). Defaults to U
#' @return a tibble
#' @keywords internal
filter_hamstr_acc_rates <- function(hamstr_acc_rates, tau = 0, kern = c("U", "G", "BH")){
  
  kern <- match.arg(kern)
  
  # scale tau by delta_c
  tau_scl <- tau / diff(hamstr_acc_rates$depth[1:2])
  
  
  if (tau > 0){
    
    if (kern == "G"){
      
      fl <- ceiling(2*tau_scl)+1
      z <- (-fl):fl
      f <- stats::dnorm(z, 0, tau_scl)
      f <- f / sum(f)  
      
    } else if (kern == "BH"){
      
      fl <- (3 * tau_scl+1)
      z  <- (-3*tau_scl):(3*tau_scl)
      
      # PDF of Berger and Heath
      f  <- stats::dexp(z+tau_scl, 1/tau_scl)
      f <- f / sum(f)
    }  else if (kern == "U"){
      
      fl <- ceiling(tau_scl/2)
      z  <- (-tau_scl/2):(tau_scl/2)
      
      f <- stats::dunif(z, (-tau_scl/2), tau_scl/2)
      f <- f / sum(f)
    }
    
    #if (tau_scl > 0) plot(z, f)
    
  }
  
  x_sum <- hamstr_acc_rates %>% 
    dplyr::group_by(.data$iter) %>% 
    dplyr::mutate(time_per_depth = if (tau_scl == 0) {.data$time_per_depth} else {
      # pad the vector with reversed head and tail
      # fl is half the filter length
      stats::filter(
        c(rev(utils::head(.data$time_per_depth, fl)), .data$time_per_depth, rev(utils::tail(.data$time_per_depth, fl))),
        filter = f)[(fl+1):(dplyr::n()+fl)]
    }) %>% 
    dplyr::mutate(depth_per_time = 1000 * 1/.data$time_per_depth, 
                  tau = tau)
  
  x_sum
  
}

#' Summarise accumulation rates
#'
#' @param hamstr_fit a hamstr_fit object
#' @inheritParams filter_hamstr_acc_rates
#' @return a tibble containing the summarised posterior accumulation rates
#' @keywords internal
summarise_hamstr_acc_rates <- function(hamstr_fit,
                                       tau = 0,
                                       kern =  c("U", "G", "BH"),
                                       probs = c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975)
                                       ){

  kern <- match.arg(kern)
  
  x <- stats::predict(hamstr_fit, type = "acc_rates", tau = tau, kern = kern) %>%
    tidyr::pivot_longer(cols = c("time_per_depth", "depth_per_time"),
                        names_to = "acc_rate_unit") 
  
  x_sum <- x %>%
    dplyr::group_by(.data$depth, 
                    .data$c_depth_top, .data$c_depth_bottom,
                    .data$acc_rate_unit, .data$idx, .data$tau) %>%
    summarise_q(var = .data$value, probs = probs) %>% 
    dplyr::ungroup() %>%
    dplyr::arrange(.data$acc_rate_unit, .data$depth)
  
  return(x_sum)

}





