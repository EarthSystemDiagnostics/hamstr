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
plot_adam <- function(adam_fit, type = c("ribbon", "spaghetti"), n.iter = 1000, plot_diagnostics = FALSE) {

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


plot_acc_rate_pars <- function(adam_fit){


  K1.df <- rstan::summary(adam_fit$fit, par = "section_acc_mean")
  K1.df <- as_tibble(K1.df$summary, rownames = "par") %>%
    separate(par, into = c("par", "K1"), sep = "\\[") %>%
    separate(K1, into = c("K1", "h"), sep = "\\]") %>%
    mutate(K1 = factor(as.numeric(K1), ordered = T))


  K.df <- rstan::summary(adam_fit$fit, par = "alpha")
  K.df <- as_tibble(K.df$summary, rownames = "par")%>%
    separate(par, into = c("par", "K"), sep = "\\[") %>%
    separate(K, into = c("K", "h"), sep = "\\]")

  ind <- tibble(K1 = factor(adam_fit$data$whichK1, ordered = T),
                K = factor(adam_fit$data$c, ordered = T))

  K.df <- left_join(ind, K.df) %>%
    mutate(K1 = factor(K1, ordered = T, levels = 1:20))

  K.df %>%
    ggplot(aes(x = mean)) +
    #geom_histogram(position = "identity", alpha = 0.5)+
    geom_histogram(aes(fill = K1, colour = K1), # position = "identity",
                   alpha = 1) +
    geom_vline(data = K1.df, aes(xintercept = mean, colour = K1)) +
    # scale_fill_brewer(type = "qual", palette = "Set3",
    #                   guide = guide_legend(ncol = 2)) +
    # scale_color_brewer(type = "qual", palette = "Set3",
    #                    guide = guide_legend(ncol = 2)) +
    labs(x = "Acc. rate") +
    theme_bw()

}

#plot_acc_rate_pars(sarn.fit3)

plot_acc_rate_prior_posterior <- function(adam_fit){

  if (exists("K1", where = adam_fit$data)){
    p.acc <- plot_acc_rate_pars(adam_fit)
  } else {
    acc.rng <- qgamma(c(0.000001, 0.9999), shape = adam_fit$data$acc_alpha,
                      rate = adam_fit$data$acc_beta)

    acc.prior <- tibble(acc.rate = seq(acc.rng[1], acc.rng[2], length.out = 1000)) %>%
      mutate(acc.dens = dgamma(acc.rate, shape = adam_fit$data$acc_alpha,
                               rate = adam_fit$data$acc_beta))

    acc.post <- tibble(alpha = as.vector(rstan::extract(adam_fit$fit, "alpha")$alpha))

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


  s1 <- rstan::summary(adam_fit$fit)
  s1 <- as_tibble(s1$summary, rownames = "par")

  infl <- s1 %>%
    filter(grepl("infl[", par, fixed = T),
           grepl("err_infl", par, fixed = T) == FALSE)

  if (nrow(infl) > 0){
    obs_ages <- obs_ages %>%
    mutate(infl_fac = infl$mean,
           infl_err = err + err * infl_fac,
           age_lwr_infl = age + 2*infl_err,
           age_upr_infl = age - 2*infl_err)
  }


  p.fit <- posterior_ages %>%
    dplyr::filter(iter %in% sample(unique(.$iter), n.iter, replace = FALSE)) %>%
    ggplot2::ggplot(aes(x = depth, y = age, group = iter))


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
      aes(ymax = age_upr, ymin = age_lwr),
      group = NA,
      colour = "Red",
      size = 1.2,
      alpha = 1) +
    ggplot2::geom_point(
      data = obs_ages,
      aes(y = age),
      group = NA,
      colour = "Red",
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

  alph$depth <- c(mean(adam_fit$data$modelled_depths),
                  unlist(sapply((
                    hierarchical_depths(adam_fit$data)
                  ),
                  function(x) {
                    y <- stats::filter(x, rep(1 / 2, 2))
                    y[is.na(y) == FALSE]
                  })))

  gg <- alph %>%
    mutate(alpha_idx = as.numeric(alpha_idx)) %>%
    ggplot(aes(x = depth, y = mean, colour = lvl)) +
    #geom_point() +
    geom_segment(data = filter(alph, lvl == 1),
                 aes(y = mean, yend = mean,
                     x = min(alph$depth), xend = max(alph$depth))) +
    geom_line() +

    expand_limits(y = 0) +
    labs(y = "Accummulation rate [age/depth]", x = "Depth",
         colour = "Hierarchical\nlevel") +
    #scale_colour_brewer(palette = "YlOrRd") +
    #scale_colour_viridis_d() +
    #theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "top")

  return(gg)
}


