#' Plot an adam_fit object
#'
#' @param adam_fit The object returned from \code{stan_bacon}.
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
plot_adam <- function(adam_fit, type = c("ribbon", "spaghetti"), n.iter = 1000, plot_diagnostics = TRUE) {

  type <- match.arg(type)

  if (type == "ribbon"){
    p.fit <- plot_summary_age_models(adam_fit)
  } else if (type == "spaghetti"){
    p.fit <- plot_age_models(adam_fit, n.iter = n.iter)
  }

  if (plot_diagnostics == FALSE) return(p.fit)

  if (plot_diagnostics){
    p.mem <- plot_memory_prior_posterior(adam_fit)
    p.acc <- plot_hierarchical_acc_rate(adam_fit)
    }

  t.lp <- rstan::traceplot(adam_fit$fit, pars = c("lp__"), include = TRUE) +
    theme(legend.position = "top") +
    labs(x = "Iteration")

  ggpubr::ggarrange(
    ggpubr::ggarrange(t.lp, p.acc, p.mem, ncol = 3, widths = c(3,3,2)),
    p.fit,
    nrow = 2, heights = c(1, 2))
}


#' Title
#'
#' @param prior 
#' @param posterior 
#'
#' @return
#' @keywords internal
#'
#' @examples
plot_prior_posterior_hist <- function(prior, posterior){
  clrs <- c("Posterior" = "Blue", "Prior" = "Red")
  ggplot() +
    geom_histogram(data = posterior,
                   aes(x = x, after_stat(density),
                       fill = "Posterior"),
                   alpha = 0.5, bins = 100) +
    geom_line(data = prior, aes(x = x, y = d, colour = "Prior")) +
    facet_wrap(~par, scales = "free") +
    scale_fill_manual(values = clrs) +
    scale_colour_manual(values = clrs) +
    guides(fill = guide_legend(override.aes = list(alpha = c(0.5))))+
    labs(
      x = "Value",
      y = "Density",
      colour = "",
      fill = ""
    ) +
    theme_bw() 
}


#' Title
#'
#' @param adam_fit 
#'
#' @return
#' @export
#' @import rstan 
#'
#' @examples
plot_infl_prior_posterior <- function(adam_fit){
  
  clrs <- c("Posterior" = "Blue", "Prior" = "Red")
  
  adam_dat <- adam_fit$data
  
  infl_mean_shape_post <-
    tibble(infl_mean = as.vector(rstan::extract(adam_fit$fit, "infl_mean")[[1]]),
           infl_shape = as.vector(rstan::extract(adam_fit$fit, "infl_shape")[[1]])) %>% 
    mutate(iter = 1:n())
  
  
  max_x_shape <- with(adam_dat, {
    infl_shape_prior_upr <- qgamma(c(0.99), shape = infl_shape_shape, rate =  infl_shape_shape / infl_shape_mean)
  
     max(c(infl_shape_prior_upr, infl_mean_shape_post$infl_shape))
  }) 
  
  max_x_mean <- with(adam_dat, {
    infl_mean_prior_upr <- qnorm(c(0.99), 0, infl_sigma_sd)
    max(c(infl_mean_prior_upr, infl_mean_shape_post$infl_mean))
  })
  
  
  infl_fac <- rstan::extract(adam_fit$fit, "infl")[[1]] %>% 
    as_tibble() %>% 
    gather() %>% 
    mutate(key = readr::parse_number(key))
  
 
  p.infl.fac <- rstan::stan_plot(adam_fit$fit, pars = "infl")
  
 
  
  infl_prior_shape <-
    tibble(x = seq(0, max_x_shape, length.out = 1000)) %>%
    mutate(
      #infl_mean = 2*dnorm(x, 0,  sd = adam_dat$infl_sigma_sd),
      d = dgamma(x-1, adam_dat$infl_shape_shape,  rate = adam_dat$infl_shape_shape / adam_dat$infl_shape_mean),
      par = "infl_shape")
  
  infl_prior_mean <-
    tibble(x = seq(0, max_x_mean, length.out = 1000)) %>%
    mutate(
      d = 2*dnorm(x, 0,  sd = adam_dat$infl_sigma_sd),
      par = "infl_mean")
  
  infl_priors <- bind_rows(infl_prior_mean, infl_prior_shape)
  
  infl_mean_shape_post_long <- infl_mean_shape_post %>% 
    gather(par, x, -iter)
 
  p.pars <- plot_prior_posterior_hist(infl_priors, infl_mean_shape_post_long)
  
  
  infl_mean_shape_post <- infl_mean_shape_post %>% 
    mutate(q99 = qgamma(0.75, infl_shape, rate = infl_shape/infl_mean))
  
  
  infl_pars_prior_dist <- infl_mean_shape_post %>% 
    filter(iter %in% sample.int(n(), 10)) %>% 
    crossing(., tibble(x = exp(seq(log(0.01), log(quantile(infl_mean_shape_post$q99, prob = 0.95)), length.out = 100)))) %>% 
    mutate(d = dgamma(x, shape = infl_shape, rate = infl_shape / infl_mean),
           #d = dgamma(x, shape = infl_shape, rate = infl_shape / 1),
           par = "Modelled prior for infl_fac")
  
  
  p.priors <- infl_pars_prior_dist  %>% 
    ggplot(aes(x = x, y = d, group = iter)) +
    geom_line(alpha = 1,#/sqrt(100),
              colour = "Red") +
    theme_bw() +
    #facet_wrap(~par+iter, scales = "free") +
    labs(y = "Density", x = "Value") 
  
  p <- egg::ggarrange(plots = list(p.pars, p.priors,
                                   #p.sigma, 
                                   p.infl.fac), ncol = 2)
  
  return(p)
  
}





#' Plot Mean Accumulation Rate Prior and Posterior Distributions
#'
#' @param adam_fit 
#'
#' @return
#' @export
plot_acc_mean_prior_posterior <- function(adam_fit) {
  clrs <- c("Posterior" = "Blue", "Prior" = "Red")
  
  adam_dat <- adam_fit$data
  
  prior_mean <- adam_dat$acc_mean_prior
  
  acc_prior_rng <- qnorm(c(0.99), mean = 0, sd = 10 * prior_mean)
  
  acc_prior <-
    tibble(acc_rate = seq(0, acc_prior_rng[1], length.out = 1000)) %>%
    mutate(
      density = 2 * dnorm(acc_rate, 0, 10 * prior_mean),
      density = ifelse(acc_rate <= 0, 0, density)
    )
  
  acc_post <-
    tibble(alpha = as.vector(rstan::extract(adam_fit$fit, "alpha[1]")[[1]]))
  
  p <- acc_prior %>%
    ggplot(aes(x = acc_rate, y = density)) +
    # plot the posterior first
    geom_histogram(
      data = acc_post,
      aes(x = alpha, after_stat(density), fill = "Posterior"),
      inherit.aes = FALSE,
      alpha = 0.5,
      # set the colour for the outline of the bins but don't include in colour 
      # legend
      colour = clrs["Posterior"],
      bins = 100
    ) +
    geom_line(aes(colour = "Prior")) +
    labs(
      x = "Mean accumulation rate",
      y = "Density",
      colour = "",
      fill = ""
    ) +
    scale_fill_manual(values = clrs) +
    scale_colour_manual(values = clrs) +
    guides(fill = guide_legend(override.aes = list(alpha = c(0.5)))) +
    
    theme_bw()
  
  return(p)
  
}


plot_memory_prior_posterior <- function(adam_fit){
  # memory prior
  mem.prior <- tibble(mem = seq(0, 1, length.out = 1000)) %>%
    mutate(mem.dens = dbeta(mem, shape1 = adam_fit$data$mem_alpha,
                            shape2 = adam_fit$data$mem_beta))

  w <- rstan::extract(adam_fit$fit, "w")$w
  ifelse(is.matrix(w), w <- apply(w, 1, median),  w <- as.vector(w))

  mem.post <- tibble(w = w,
                     R = as.vector(rstan::extract(adam_fit$fit, "R")$R))


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

add_subdivisions <- function(gg, adam_fit){

  tick_dat <- hierarchical_depths(adam_fit$data)

  for (x in seq_along(tick_dat)){

    df <- data.frame(x = tick_dat[[x]])

    lnth <- length(tick_dat) - (x-1)

    gg <- gg + geom_rug(data = df, aes(x = x),
                        inherit.aes = F, sides = "top",
                        length = unit(0.01*lnth, "npc"))

  }

  return(gg)
}


plot_age_models <- function(adam_fit, n.iter = 1000){


  posterior_ages <- get_posterior_ages(adam_fit)

  obs_ages <- dplyr::tibble(
    depth = adam_fit$data$depth,
    age = adam_fit$data$obs_age,
    err = adam_fit$data$obs_err)

  obs_ages <- dplyr::mutate(obs_ages,
                            age_upr = age + 2*err,
                            age_lwr = age - 2*err)
  
  infl_errs <- rstan::summary(adam_fit$fit, par = "obs_err_infl")$summary %>% 
    as_tibble(., rownames = "par") %>% 
    mutate(dat_idx = readr::parse_number(par))
  
  p.fit <- posterior_ages %>%
    dplyr::filter(iter %in% sample(unique(.$iter), n.iter, replace = FALSE)) %>%
    ggplot2::ggplot(aes(x = depth, y = age, group = iter))


  p.fit <- p.fit +
    ggplot2::geom_line(alpha = 0.5 / sqrt(n.iter))

  if (adam_fit$data$inflate_errors == 1){
    obs_ages <- obs_ages %>% 
      mutate(infl_err = infl_errs$mean,
             age_lwr_infl = age + 2*infl_err,
             age_upr_infl = age - 2*infl_err)
    
    p.fit <- p.fit +
      ggplot2::geom_linerange(
        data = obs_ages,
        aes(x = depth, ymax = age_upr_infl, ymin = age_lwr_infl),
        group = NA,
        colour = "Red",
        alpha = 0.5, inherit.aes = F)
  }
  
  p.fit <- p.fit +
    ggplot2::geom_linerange(
      data = obs_ages,
      aes(ymax = age_upr, ymin = age_lwr),
      group = NA,
      colour = "Blue",
      size = 1.2,
      alpha = 1) +
    ggplot2::geom_point(
      data = obs_ages,
      aes(y = age),
      group = NA,
      colour = "Blue",
      #size = 1.01,
      alpha = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    labs(x = "Depth", y = "Age")


  # add subdivisions
  p.fit <- add_subdivisions(p.fit, adam_fit)

  return(p.fit)

}


#' Plot the hierarchical accumulation rate parameters
#'
#' @param adam_fit
#'
#' @return
#' @export
#'
#' @examples
plot_hierarchical_acc_rate <- function(adam_fit){

  idx <- as_tibble(alpha_indices(adam_fit$data$K)[1:3]) %>%
    mutate(alpha_idx = (alpha_idx))

  a3 <- rstan::summary(adam_fit$fit, pars = "alpha")$summary
  
  alph <- as_tibble(a3, rownames = "par") %>%
    mutate(alpha_idx = readr::parse_number(par)) %>%
    left_join(idx, .) %>%
    mutate(lvl = factor(lvl))

  # for each unit at each level in hierarchy get max and min depth 
  alph$depth1 <- c(min(adam_fit$data$modelled_depths),
                   unlist(sapply((
                     hierarchical_depths(adam_fit$data)
                   ),
                   function(x) {
                    head(x, -1)
                   })))
  
  alph$depth2 <- c(max(adam_fit$data$modelled_depths),
                   unlist(sapply((
                     hierarchical_depths(adam_fit$data)
                   ),
                   function(x) {
                     tail(x, -1)
                   })))
  
  alph2 <- alph %>% 
    select(lvl, alpha_idx, depth1, depth2, mean) %>% 
    group_by(lvl) %>% 
    gather(type, depth, -mean, -lvl, -alpha_idx) %>% 
    select(lvl, alpha_idx, depth, mean) %>% 
    arrange(lvl, alpha_idx, depth, mean)
  

  gg <- alph2 %>%
    ggplot(aes(x = depth, y = mean, colour = lvl)) +
    geom_path() +
    expand_limits(y = 0) +
    labs(y = "Accummulation rate [age/depth]", x = "Depth",
         colour = "Hierarchical\nlevel") +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "top")

  return(gg)
}
