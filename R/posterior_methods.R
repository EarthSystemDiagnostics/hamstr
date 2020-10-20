# Extractor functions

#' Get Posterior Age Models
#'
#' @inheritParams plot_hamstr
#'
#' @return a dataframe/tibble with posterior ages for all iterations after warmup
#' @export
#' @importFrom readr parse_number
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error,
#'   K = c(10, 10), nu = 6,
#'   acc_mean_prior = 20,
#'   mem_mean = 0.7, mem_strength = 4,
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
    dplyr::left_join(depths, .data$.) %>% 
    dplyr::arrange(.data$par, .data$iter, .data$idx, .data$depth)
  
  return(posterior_ages)
  
}


#' Interpolate Posterior Age Model At New Depths
#'
#' @inheritParams plot_hamstr
#' @param new_depth a vector of depths at which to interpolate the age models
#'
#' @return hamstr_interpolated_ages object
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error,
#'   K = c(10, 10), nu = 6,
#'   acc_mean_prior = 20,
#'   mem_mean = 0.7, mem_strength = 4,
#'   inflate_errors = 0,
#'   iter = 2000, chains = 3)
#'   
#' interpolate.age.models(fit, new_depth = seq(1000, 15000, by = 1000))
#' }
#' 
interpolate_age_models <- function(hamstr_fit, new_depth){
  
  # get posterior age models
  pst_age <- get_posterior_ages(hamstr_fit)
  
  new_age <- pst_age %>% 
    dplyr::group_by(iter) %>% 
    dplyr::do({
      tibble::tibble(
        iter = .$iter[1],
        depth = new_depth,
        age = stats::approx(.$depth, .$age, new_depth)$y
      )
    }) %>% 
    dplyr::ungroup()
  
  class(new_age) <- append("hamstr_interpolated_ages", class(new_age))
  
  return(new_age)
}


#' Title
#'
#' @param object 
#' @param type 
#'
#' @return
#'
#' @examples
#' @export
#' @method summary hamstr_fit
summary.hamstr_fit <- function(object, type = "age_models"){
  if (type == "age_models"){
    summarise_age_models(object)
  }
}


#' Title
#'
#' @param object 
#' @param type 
#'
#' @return
#'
#' @examples
#' @export
#' @method summary hamstr_interpolated_ages
summary.hamstr_interpolated_ages <- function(object, type = "age_models"){
 
   if (type == "age_models"){
    summarise_new_ages(object)
     }
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
#' @description Extracts the summary statistics of posterior age models and attached the depths 
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
#'   mem_mean = 0.7, mem_strength = 4,
#'   inflate_errors = 0,
#'   iter = 2000, chains = 3)
#'   
#' summarise_age_models(fit)
#' }
summarise_age_models <- function(hamstr_fit){
  
  if (is_hamstr_interpolated_ages(hamstr_fit)){
    age_summary <-summarise_new_ages(hamstr_fit)
  } else {
    age_summary <- rstan::summary(hamstr_fit$fit, par = "c_ages")[["summary"]] %>% 
      tibble::as_tibble(., rownames = "par")
    
    depths <- tibble::tibble(depth = hamstr_fit$data$modelled_depths,
                     idx = 1:length(hamstr_fit$data$modelled_depths))
    age_summary <- age_summary %>% 
      dplyr::mutate(idx = readr::parse_number(par)) %>% 
      dplyr::left_join(depths, .)
    
  }
  return(age_summary)
}




