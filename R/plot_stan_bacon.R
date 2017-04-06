#' Title
#'
#' @param stan_bacon_fit The object returned from \code{stan_bacon}.
#' @param n.iter The number of iterations of the model to plot, defaults to 100.
#' 
#' @description Plots the Bacon modelled Age~Depth relationship together with
#'   the depths, ages, and age uncertainties in the observed data. A sample of
#'   size \code{n.iter} of the interactions of the posterior distribution are
#'   plotted as grey lines. The observed data are plotted as points with +- 2*se error bars.
#'
#' @return A ggplot2 object
#' @export
#' @examples
plot_stan_bacon <- function(stan_bacon_fit, n.iter = 100) {
  
  fit_data <- stan_bacon_fit$data
  
  c_depths <- c(fit_data$c_depth_top[1], fit_data$c_depth_bottom) %>%
    dplyr::tibble(Depth = ., Section = paste0("V", 1:length(.)))
  
  obs_ages <- dplyr::tibble(
    Depth = fit_data$depth,
    Age = fit_data$obs_age,
    Err = fit_data$obs_err
  ) %>%
    dplyr::mutate(Age_upr = Age + 2*Err,
           Age_lwr = Age - 2*Err)
  
  posterior <- rstan::extract(stan_bacon_fit$fit)
  post_depth_age <- posterior$c_ages %>%
    dplyr::tbl_df() %>%
    dplyr::mutate(Iter = 1:nrow(posterior$c_ages)) %>%
    tidyr::gather(Section, Age,-Iter) %>%
    dplyr::left_join(., c_depths, by = "Section")
  
  p <- post_depth_age %>%
    dplyr::filter(Iter %in% sample(unique(.$Iter), n.iter, replace = FALSE)) %>%
    ggplot2::ggplot(aes(x = Depth, y = Age, group = Iter)) +
    ggplot2::geom_line(alpha = 1 / sqrt(n.iter))  +
    ggplot2::geom_pointrange(
      data = obs_ages,
      aes(y = Age, ymax = Age_upr, ymin = Age_lwr),
      group = NA,
      colour = "Red",
      alpha = 0.5) +
    ggplot2::theme_bw()
  
  p
}

