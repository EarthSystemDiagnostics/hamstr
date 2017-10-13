#' Title
#'
#' @param stan_bacon_fit The object returned from \code{stan_bacon}.
#' @param n.iter The number of iterations of the model to plot, defaults to 1000.
#' 
#' @description Plots the Bacon modelled Age~Depth relationship together with
#'   the depths, ages, and age uncertainties in the observed data. A sample of
#'   size \code{n.iter} of the interactions of the posterior distribution are
#'   plotted as grey lines. The observed data are plotted as points with +- 2*se error bars.
#'
#' @return A ggplot2 object
#' @export
#' @importFrom ggpubr ggarrange
#' @importFrom rstan extract
#' @examples
plot_stan_bacon <- function(stan_bacon_fit, n.iter = 1000) {
  
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
  
  p.fit <- post_depth_age %>%
    dplyr::filter(Iter %in% sample(unique(.$Iter), n.iter, replace = FALSE)) %>%
    ggplot2::ggplot(aes(x = Depth, y = Age, group = Iter)) +
    ggplot2::geom_line(alpha = 0.5 / sqrt(n.iter))  +
    ggplot2::geom_pointrange(
      data = obs_ages,
      aes(y = Age, ymax = Age_upr, ymin = Age_lwr),
      group = NA,
      colour = "Red",
      alpha = 0.5) +
    ggplot2::theme_bw()
  
  
  ## Prior and posterior figures
  acc.rng <- qgamma(c(0.000001, 0.999), shape = stan_bacon_fit$data$acc_alpha,  rate = stan_bacon_fit$data$acc_beta)
  
  acc.prior <- tibble(acc.rate = seq(acc.rng[1], acc.rng[2], length.out = 1000)) %>% 
    mutate(acc.dens = dgamma(acc.rate, shape = stan_bacon_fit$data$acc_alpha,  rate = stan_bacon_fit$data$acc_beta))
  
  acc.post <- tibble(alpha = as.vector(rstan::extract(stan_bacon_fit$fit, "alpha")$alpha))
  
  p.acc <- acc.prior %>% 
    ggplot(aes(x = acc.rate, y = acc.dens)) +
    geom_density(data = acc.post, aes(x = alpha), fill = "Grey", inherit.aes = FALSE) +
    geom_line(colour = "Red") + 
    scale_x_continuous("Acc. rate [yr/cm]", limits = acc.rng) +
    scale_y_continuous("Density") +
    theme_bw()
  
  # memory prior
  mem.prior <- tibble(mem = seq(0, 1, length.out = 1000)) %>% 
    mutate(mem.dens = dbeta(mem, shape1 = stan_bacon_fit$data$mem_alpha,  shape2 = stan_bacon_fit$data$mem_beta))
  
  mem.post <- tibble(R = as.vector(rstan::extract(stan_bacon_fit$fit, "R")$R),
                     w = as.vector(rstan::extract(stan_bacon_fit$fit, "w")$w))
  
  p.mem <- mem.prior %>% 
    ggplot(aes(x = mem, y = mem.dens)) +
    geom_density(data = mem.post, aes(x = R), fill = "Blue", inherit.aes = FALSE) +
    geom_density(data = mem.post, aes(x = w), fill = "Grey", inherit.aes = FALSE) +
    geom_line(colour = "Red") + 
    scale_x_continuous("Memory [correlation at 1 cm]", limits = c(0, 1)) + 
    scale_y_continuous("") +
    theme_bw()
  
  ggpubr::ggarrange(
    ggpubr::ggarrange(p.acc, p.mem, ncol = 2),
    p.fit,
    nrow = 2, heights = c(1, 2)) 
}

