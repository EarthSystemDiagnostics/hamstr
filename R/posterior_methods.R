# Extractor functions

#' Get Posterior Age Models
#'
#' @param adam_fit fitted adam model, output from adam()
#'
#' @return a dataframe/tibble with posterior ages for all iterations after warmup
#' @export
#'
#' @examples
#' \dontrun{
#' get_posterior_ages(fitted.adam.object)
#' }
get_posterior_ages <- function(adam_fit){
  
  depths <- tibble(depth = adam_fit$data$modelled_depths,
                   idx = 1:length(adam_fit$data$modelled_depths))
  
  posterior_ages <- as.data.frame(adam_fit$fit, pars = "c_ages") %>% 
    as_tibble() %>% 
    mutate(iter = 1:nrow(.)) %>% 
    gather(par, age, -iter) %>% 
    mutate(idx = readr::parse_number(par),
           par = "c_ages") %>% 
    left_join(depths, .) %>% 
    arrange(par, iter, idx, depth)
  
  return(posterior_ages)
  
}


#' Interpolate Posterior Age Model At New Depths
#'
#' @param adam_fit fitted adam model, output from adam()
#' @param new_depth a vector of depths at which to interpolate the age models
#'
#' @return 
#' @export
#'
#' @examples
#' \dontrun{
#' interpolate.age.models(fitted.adam.object, new_depth = seq(1000, 15000, by = 1000))
#' }
#' 
interpolate_age_models <- function(adam_fit, new_depth){
  
  # get posterior age models
  pst_age <- get_posterior_ages(adam_fit)
  
  new_age <- pst_age %>% 
    group_by(iter) %>% 
    do({
      tibble(
        iter = .$iter[1],
        depth = new_depth,
        age = approx(.$depth, .$age, new_depth)$y
      )
    })
  
  class(new_age) <- append(class(new_age), "adam_new_ages")
  
  return(new_age)
}



#' Summarise Interpolated Posterior Age Models
#'
#' @param new_ages 
#'
#' @return
#' @export
#'
#' @examples
summarise_new_ages <- function(new_ages){
  
  new_ages_sum <- new_ages %>% 
    group_by(depth) %>% 
    summarise(mean = mean(age),
              #se_mean = NA,
              sd = sd(age),
              `2.5%` = quantile(age, probs = c(0.025), na.rm = T),
              `25%` = quantile(age, probs = c(0.25), na.rm = T),
              `50%` = quantile(age, probs = c(0.50), na.rm = T),
              `75%` = quantile(age, probs = c(0.75), na.rm = T),
              `97.5%` = quantile(age, probs = c(0.975), na.rm = T))
  
  return(new_ages_sum)
  
}



#' Summarise Posterior Age Models
#'
#' @param adam_fit 
#' @description Extracts the summary statistics of posterior age modells and attached the depths 
#' @return
#' @export
#'
#' @examples
summarise_age_models <- function(adam_fit){
  
  age_summary <- summary(adam_fit$fit, par = "c_ages")$summary %>% 
    as_tibble(., rownames = "par")
  
  depths <- tibble(depth = adam_fit$data$modelled_depths,
                   idx = 1:length(adam_fit$data$modelled_depths))
  age_summary <- age_summary %>% 
    mutate(idx = readr::parse_number(par)) %>% 
    left_join(depths, .)
  
  return(age_summary)
}



#' Plot Summary of Posterior Age Models
#'
#' @param adam_fit 
#'
#' @return
#' @export
#'
#' @examples
plot_summary_age_models <- function(adam_fit){
  
  age_summary <- summarise_age_models(adam_fit)

  obs_ages <- data.frame(
    depth = adam_fit$data$depth,
    age = adam_fit$data$obs_age,
    err = adam_fit$data$obs_err)
  
  obs_ages <- dplyr::mutate(obs_ages,
                            age_upr = age + 2*err,
                            age_lwr = age - 2*err)
  
  
  infl_errs <- summary(adam_fit$fit, par = "obs_err_infl")$summary %>% 
    as_tibble(., rownames = "par") %>% 
    mutate(dat_idx = readr::parse_number(par))
  
  p.age.sum <- age_summary %>% 
    ggplot(aes(x = depth, y = mean)) +
    geom_ribbon(aes(ymax = `2.5%`, ymin = `97.5%`), fill = "Lightgrey") +
    geom_ribbon(aes(ymax = `75%`, ymin = `25%`), fill = "Darkgrey") +
    geom_line() +
    geom_line(aes(y = `50%`), colour = "Green") +
    labs(x = "Depth", y = "Age") +
    theme_bw() 
   
   
  if (adam_fit$data$inflate_errors == 1){
    obs_ages <- obs_ages %>% 
      mutate(infl_err = infl_errs$mean,
             age_lwr_infl = age + 2*infl_err,
             age_upr_infl = age - 2*infl_err)
    
   p.age.sum <- p.age.sum +
    ggplot2::geom_linerange(
      data = obs_ages,
      aes(x = depth, ymax = age_upr_infl, ymin = age_lwr_infl),
      group = NA,
      colour = "Red",
      alpha = 0.5, inherit.aes = F)
  }
  
  p.age.sum <- p.age.sum +
    geom_linerange(data = obs_ages,
                   aes(y = age, ymax = age_upr, ymin = age_lwr),
                    colour = "Blue", size = 1.5) +
    geom_point(data = obs_ages, aes(y = age),
               colour = "Yellow")
  
  p.age.sum
}


