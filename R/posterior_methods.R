# Extractor functions

#' Get Posterior Age Models
#'
#' @param hamstr_fit fitted hamstr model, output from hamstr()
#'
#' @return a dataframe/tibble with posterior ages for all iterations after warmup
#' @export
#'
#' @examples
#' \dontrun{
#' get_posterior_ages(fitted.hamstr.object)
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
    dplyr::left_join(depths, .) %>% 
    dplyr::arrange(par, iter, idx, depth)
  
  return(posterior_ages)
  
}


#' Interpolate Posterior Age Model At New Depths
#'
#' @param hamstr_fit fitted hamstr model, output from hamstr()
#' @param new_depth a vector of depths at which to interpolate the age models
#'
#' @return 
#' @export
#'
#' @examples
#' \dontrun{
#' interpolate.age.models(fitted.hamstr.object, new_depth = seq(1000, 15000, by = 1000))
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
    })
  
  class(new_age) <- append(class(new_age), "hamstr_interpolated_ages")
  
  return(new_age)
}


#' Summarise Interpolated Posterior Age Models
#'
#' @param new_ages 
#'
#' @return
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
#' @param hamstr_fit an hamstr_fit object or hamstr_interpolated_ages object
#' @description Extracts the summary statistics of posterior age models and attached the depths 
#' @return
#' @export
#'
#' @examples
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



#' Plot Summary of Posterior Age Models
#'
#' @param hamstr_fit 
#'
#' @return
#' @export
#'
#' @examples
plot_summary_age_models <- function(hamstr_fit){
  
  age_summary <- summarise_age_models(hamstr_fit)

  obs_ages <- data.frame(
    depth = hamstr_fit$data$depth,
    age = hamstr_fit$data$obs_age,
    err = hamstr_fit$data$obs_err)
  
  obs_ages <- dplyr::mutate(obs_ages,
                            age_upr = age + 2*err,
                            age_lwr = age - 2*err)
  
  
  infl_errs <- rstan::summary(hamstr_fit$fit, par = "obs_err_infl")$summary %>% 
    tibble::as_tibble(., rownames = "par") %>% 
    dplyr::mutate(dat_idx = readr::parse_number(par))
  
  p.age.sum <- age_summary %>% 
    ggplot2::ggplot(ggplot2::aes(x = depth, y = mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymax = `2.5%`, ymin = `97.5%`), fill = "Lightgrey") +
    ggplot2::geom_ribbon(ggplot2::aes(ymax = `75%`, ymin = `25%`), fill = "Darkgrey") +
    ggplot2::geom_line() +
    ggplot2::geom_line(ggplot2::aes(y = `50%`), colour = "Green") +
    ggplot2::labs(x = "Depth", y = "Age") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
   
   
  if (hamstr_fit$data$inflate_errors == 1){
    obs_ages <- obs_ages %>% 
      dplyr::mutate(infl_err = infl_errs$mean,
             age_lwr_infl = age + 2*infl_err,
             age_upr_infl = age - 2*infl_err)
    
   p.age.sum <- p.age.sum +
    ggplot2::geom_linerange(
      data = obs_ages,
      ggplot2::aes(x = depth, ymax = age_upr_infl, ymin = age_lwr_infl),
      group = NA,
      colour = "Red",
      alpha = 0.5, inherit.aes = F)
  }
  
  p.age.sum <- p.age.sum +
    ggplot2::geom_linerange(data = obs_ages,
                   ggplot2::aes(x = depth, 
                       ymax = age_upr, ymin = age_lwr), inherit.aes = FALSE,
                    colour = "Blue", size = 1.25) +
    ggplot2::geom_point(data = obs_ages, ggplot2::aes(y = age),
               colour = "Blue")
  
  
  p.age.sum <- add_subdivisions(p.age.sum, hamstr_fit)
  
  p.age.sum
}


