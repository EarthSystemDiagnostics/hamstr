#' Title
#'
#' @param stan_fit The object returned from \code{stan_bacon}.
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
#' @importFrom magrittr %>%
#' @examples
plot_stan_bacon <- function(stan_fit, n.iter = 1000, plot_priors = TRUE) {

  p.fit <- plot_age_models(stan_fit, n.iter = n.iter)

  if (plot_priors == FALSE) return(p.fit)
  
  if (plot_priors){
    p.mem <- plot_memory_prior_posterior(stan_fit)
    p.acc <- plot_acc_rate_prior_posterior(stan_fit)
    }    
  
  t.lp <- rstan::traceplot(stan_fit$fit, pars = c("lp__"), include = TRUE) +
    theme(legend.position = "top")

  ggpubr::ggarrange(
    ggpubr::ggarrange(t.lp, p.acc, p.mem, ncol = 3, widths = c(3,2,3)),
    p.fit,
    nrow = 2, heights = c(1, 2))
}


plot_acc_rate_prior_posterior <- function(stan_fit){
  
  if (exists("K1", where = stan_fit$data)){
    p.acc <- ggplot() + theme_void()
  } else {
    acc.rng <- qgamma(c(0.000001, 0.9999), shape = stan_fit$data$acc_alpha,
                      rate = stan_fit$data$acc_beta)
    
    acc.prior <- tibble(acc.rate = seq(acc.rng[1], acc.rng[2], length.out = 1000)) %>%
      mutate(acc.dens = dgamma(acc.rate, shape = stan_fit$data$acc_alpha,
                               rate = stan_fit$data$acc_beta))
    
    acc.post <- tibble(alpha = as.vector(rstan::extract(stan_fit$fit, "alpha")$alpha))
    
    if (any(acc.post$alpha > acc.rng[2])){
      acc.outside <- sum(acc.post$alpha > acc.rng[2])
      p.acc.o <- round(100 * acc.outside / length(acc.post$alpha), 2)
      warning(p.acc.o, paste0("% of accumulation rate samples were outside plotted region of the acc.rate prior"))
    }
    
    acc.post <- acc.post[acc.post$alpha <= acc.rng[2], ]
    
    p.acc <- acc.prior %>%
      ggplot(aes(x = acc.rate, y = acc.dens)) +
      geom_density(data = acc.post, aes(x = alpha), fill = "Grey", inherit.aes = FALSE) +
      geom_line(colour = "Red") +
      scale_x_continuous("Acc. rate [yr/cm]", limits = acc.rng) +
      scale_y_continuous("Density") +
      theme_bw()
  }

return(p.acc)
  
}



plot_memory_prior_posterior <- function(stan_fit){
  # memory prior
  mem.prior <- tibble(mem = seq(0, 1, length.out = 1000)) %>%
    mutate(mem.dens = dbeta(mem, shape1 = stan_fit$data$mem_alpha,
                            shape2 = stan_fit$data$mem_beta))
  
  w <- rstan::extract(stan_fit$fit, "w")$w
  ifelse(is.matrix(w), w <- apply(w, 1, median),  w <- as.vector(w))
  
  mem.post <- tibble(w = w,
                     R = as.vector(rstan::extract(stan_fit$fit, "R")$R))
  
  
  p.mem <- mem.prior %>%
    ggplot(aes(x = mem, y = mem.dens)) +
    geom_density(data = mem.post, aes(x = R, fill = "at 1 cm 'R'"),
                 inherit.aes = FALSE, show.legend = TRUE) +
    geom_density(data = mem.post, aes(x = w, fill = "between\nsections 'w'"),
                 inherit.aes = FALSE, show.legend = TRUE) +
    geom_line(colour = "Red") +
    scale_x_continuous("Memory [correlation]", limits = c(0, 1)) +
    scale_y_continuous("") +
    scale_fill_discrete("") +
    theme_bw() +
    theme(legend.position = "top")
  
  return(p.mem)
}




plot_age_models <- function(stan_fit, n.iter = 1000){
  
  fit_data <- stan_fit$data
  
  c_depths <- c(fit_data$c_depth_top[1], fit_data$c_depth_bottom) %>%
    dplyr::tibble(Depth = ., Section = paste0("V", 1:length(.)))
  
  obs_ages <- dplyr::tibble(
    Depth = fit_data$depth,
    Age = fit_data$obs_age,
    Err = fit_data$obs_err)
  
  obs_ages <- dplyr::mutate(obs_ages,
                            Age_upr = Age + 2*Err,
                            Age_lwr = Age - 2*Err)
  
  posterior <- rstan::extract(stan_fit$fit)
  
 
  s1 <- rstan::summary(stan_fit$fit)
  s1 <- as_tibble(s1$summary, rownames = "par")
  
  infl <- s1 %>% 
    filter(grepl("infl[", par, fixed = T),
           grepl("err_infl", par, fixed = T) == FALSE)
  
  if (nrow(infl) > 0){
    obs_ages <- obs_ages %>% 
    mutate(infl_fac = infl$mean,
           infl_err = Err + Err * infl_fac,
           age_lwr_infl = Age + 2*infl_err,
           age_upr_infl = Age - 2*infl_err)
  }
  
  post_depth_age <- posterior$c_ages %>%
    dplyr::as_tibble(.) %>%
    dplyr::mutate(Iter = 1:nrow(posterior$c_ages)) %>%
    tidyr::gather(Section, Age,-Iter) %>%
    dplyr::left_join(., c_depths, by = "Section")
  
  p.fit <- post_depth_age %>%
    dplyr::filter(Iter %in% sample(unique(.$Iter), n.iter, replace = FALSE)) %>%
    ggplot2::ggplot(aes(x = Depth, y = Age, group = Iter))
  
  if (exists("K1", where = fit_data)){
    
    bnds <- tapply(fit_data$c_depth_top, fit_data$whichK1, min)
    
    p.fit <- p.fit + 
      ggplot2::geom_vline(xintercept = bnds[-1], colour = "grey", linetype = 2) 
  }
  
  p.fit <- p.fit +
    ggplot2::geom_line(alpha = 0.5 / sqrt(n.iter)) 
  
  if (nrow(infl) > 0){
    p.fit <- p.fit +
    ggplot2::geom_linerange(
      data = obs_ages,
      aes(ymax = age_upr_infl, ymin = age_lwr_infl),
      group = NA,
      colour = "Blue",
      alpha = 0.5)
  }
  p.fit <- p.fit +
    ggplot2::geom_linerange(
      data = obs_ages,
      aes(ymax = Age_upr, ymin = Age_lwr),
      group = NA,
      colour = "Red",
      size = 1.2,
      alpha = 1) +
    ggplot2::geom_point(
      data = obs_ages,
      aes(y = Age),
      group = NA,
      colour = "Red",
      #size = 1.01,
      alpha = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  
  return(p.fit)
  
}





