# Methods -------

#' Plot hamstr Objects
#' 
#' @description plot method for class "hamstr_fit".
#'
#' @param x a hamstr_fit object
#' @param type one of "default", "age_models","acc_rates", "hier_acc_rates",
#' "acc_mean_prior_post", "mem_prior_post", "L_prior_post", "D_prior_post", "PDF_14C"
#' @param ... additional arguments to hamstr plotting methods
#' @inheritParams plot_hamstr
#' @return a ggplot object
#'
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error)
#'
#' plot(fit)
#' plot(fit, type = "acc_rates", tau = 5, kern = "U")
#' }
#' @export
#' @method plot hamstr_fit
plot.hamstr_fit <- function(x,
                            type = c("default",
                              "age_models",
                              "acc_rates",
                              "hier_acc_rates",
                              "acc_mean_prior_post",
                              "mem_prior_post",
                              "L_prior_post",
                              "D_prior_post",
                              "PDF_14C"
                              ),
                            summarise = TRUE,
                            ...){

  type <- match.arg(type)

  switch(type,
         default = plot_hamstr(x, summarise = summarise, ...),
         age_models = plot_hamstr(x, summarise = summarise,
                                  plot_diagnostics  = FALSE, ...),
         acc_rates = plot_hamstr_acc_rates(x, ...),
         hier_acc_rates = plot_hierarchical_acc_rate(x),
         acc_mean_prior_post = plot_acc_mean_prior_posterior(x),
         mem_prior_post = plot_memory_prior_posterior(x),
         L_prior_post = plot_L_prior_posterior(x),
         D_prior_post = plot_D_prior_posterior(x),
         PDF_14C = plot_14C_PDF(x, ...)
  )
  }



# Palettes ------

hamstr_pal <- c("Age models" = "darkgrey", "Median" = "black",
                "68%" = "darkgrey", "95%" = "lightgrey",
                "Mean" = "#A1D76A", 
                "Age point" = "#5785C1", "Obs age" = "#5785C1", "Latent age" = "#FBA72A", "Infl err" = "#FBA72A") 

add_colour_scale <- function(gg, 
                             clrs = hamstr_pal,
                             lbls = names(hamstr_pal)){
  gg +
    ggplot2::scale_fill_manual(name = "Interval",
                               values = clrs,
                               breaks = lbls,
                               labels = lbls,
                               guide = "legend") +
    ggplot2::scale_colour_manual(name = "",
                                 values = clrs,
                                 breaks = lbls,
                                 labels = lbls,
                                 guide = "legend")+ 
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  
}



# Functions ------

#' Plot an hamstr_fit object
#'
#' @param hamstr_fit the object returned from \code{stan_hamstr}.
#'
#' @param n.iter the number of iterations of the model to plot, defaults to
#'   1000.
#' @param summarise logical TRUE or FALSE. Plot the realisations as a summarised
#'  "ribbon" showing 50% and 95% intervals (faster), or as a spaghetti plot
#'  showing individual realisations. Defaults to TRUE (ribbon).
#' @param plot_diagnostics logical, include diagnostic plots: traceplot of
#'   log-posterior, hierarchical accumulations rates, memory parameter. Defaults
#'   to TRUE.
#' @description Plots the HAMStR modelled age ~ depth relationship together with
#'   the depths, ages, and age uncertainties in the observed data. A random
#'   sample of size \code{n.iter} of the iterations of the posterior
#'   distribution are plotted as grey lines. The observed data are plotted as
#'   points with +- 2*se error bars.
#'
#' @return A ggplot2 object
#' @keywords internal
#' @import patchwork
#' @importFrom rstan extract
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error)
#'
#' # With age models summarised as a ribbon. Faster than spaghetti plots.
#' plot_hamstr(fit)
#'
#' # With age models as spaghetti plots. Can see individual realisations, but slower to plot.
#' plot_hamstr(fit, summarise = FALSE)
#' }
plot_hamstr <- function(hamstr_fit, summarise = TRUE,
                        n.iter = 1000, plot_diagnostics = TRUE) {

  #summarise <- match.arg(summarise)



  if (summarise == TRUE){
    p.fit <- plot_summary_age_models(hamstr_fit)
  } else if (summarise == FALSE){
    p.fit <- plot_age_models(hamstr_fit, n.iter = n.iter)
  }

  ## Add latent variable estimated ages accounting for bioturbation error

  if (hamstr_fit$data$model_bioturbation == 1){

    post <- as.data.frame(hamstr_fit$fit, pars = c("bt_age")) %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(iter = 1:dplyr::n()) %>%
      tidyr::pivot_longer(cols = -"iter") %>%
      dplyr::mutate(dpt = get_par_idx(.data$name))

    tmp <- dplyr::tibble(depth = hamstr_fit$data$depth,
                  dpt = 1:hamstr_fit$data$N) %>%
    dplyr::left_join(post, by = "dpt")

    p.fit <- p.fit +
      ggplot2::geom_violin(data = tmp,
                           ggplot2::aes(x = .data$depth, y = .data$value,
                                        group = as.factor(.data$dpt),
                                        colour = "Latent age"), fill = NA,
                           scale = "area",
                           position = ggplot2::position_identity(), alpha = 0.75,
                           show.legend = TRUE, key_glyph = draw_key_linerange)

  }



  if (plot_diagnostics == FALSE) return(p.fit)

  if (plot_diagnostics){


    p.mem <- plot_memory_prior_posterior(hamstr_fit) +
      ggplot2::labs(x = "AR1 coefficient") +
      ggplot2::theme(panel.grid = ggplot2::element_blank())

    # only plot if model has been sampled
    if (hamstr_fit$data$sample_posterior){
     p.acc <- plot_hierarchical_acc_rate(hamstr_fit)

    t.lp <- rstan::traceplot(hamstr_fit$fit, pars = c("lp__"), include = TRUE, alpha = 0.75) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top") +
      ggplot2::labs(x = "Iteration", colour = "Chain") +
      ggplot2::scale_colour_brewer(name = "Set1", type = "qual", direction = -1) +
      ggplot2::guides(color = guide_legend(override.aes = list(alpha = 1)))
    
    diag.list <- list(t.lp, p.acc, p.mem)
    
    diag.nest <- t.lp | p.acc | p.mem
    } else {

      diag.nest <- patchwork::plot_spacer() | p.mem | patchwork::plot_spacer()

    }

   

    if (hamstr_fit$data$model_bioturbation == 1){
      p.L <- plot_L_prior_posterior(hamstr_fit) +
        theme(legend.position = "top")

        diag.nest <- diag.nest | p.L
      
    }

    if (hamstr_fit$data$model_displacement == 1){
      p.D <- plot_D_prior_posterior(hamstr_fit) +
        theme(legend.position = "top")
       
      diag.nest <- diag.nest | p.D
    } 

    p.fit <- p.fit / (diag.nest) + patchwork::plot_layout(heights = c(3, 1))
    
    return(p.fit)

  }
}


#' Title
#'
#' @param hamstr_fit 
#'
#' @return dataframe
#' @keywords internal
get_obs_age_err <- function(hamstr_fit){
  obs_ages <- data.frame(
    depth = hamstr_fit$data$depth,
    age = hamstr_fit$data$obs_age,
    err = hamstr_fit$data$obs_err)
  
  obs_ages <- dplyr::mutate(obs_ages,
                            age_upr1 = .data$age + 1*.data$err,
                            age_lwr1 = .data$age - 1*.data$err,
                            age_upr2 = .data$age + 2*.data$err,
                            age_lwr2 = .data$age - 2*.data$err)
  
  return(obs_ages)
}



#' Plot Summary of Posterior Age Models
#'
#' @inheritParams plot_hamstr
#'
#' @return A ggplot2 object
#' @keywords internal

#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error)
#'
#' plot_summary_age_models(fit)
#' }
plot_summary_age_models <- function(hamstr_fit){

 
  obs_ages <- get_obs_age_err(hamstr_fit)

  if (hamstr_fit$data$sample_posterior == FALSE){

    gg <- ggplot2::ggplot(obs_ages, aes(x = .data$depth, y = .data$age, fill = NA)) +
      geom_blank() +
      theme_bw() 

    gg <- add_datapoints(gg, obs_ages)
    gg <- add_subdivisions(gg, hamstr_fit)
    gg <- add_colour_scale(gg)

    return(gg)

  }

  age_summary <- summarise_age_models(hamstr_fit)

  p.age.sum <- plot_downcore_summary(age_summary)
  p.age.sum <- add_colour_scale(p.age.sum)

  if (hamstr_fit$data$inflate_errors == 1 | hamstr_fit$data$model_displacement == 1){

    p.age.sum <- add_infl_err(p.age.sum, hamstr_fit, obs_ages)
  }

  p.age.sum <- add_datapoints(p.age.sum, obs_ages)

  p.age.sum <- add_subdivisions(p.age.sum, hamstr_fit)

  p.age.sum
}


#' Add Age Control Points to Plot
#'
#' @param gg a ggplot2 object 
#' @param dat a dataframe containing the datapoints to add
#'
#' @return modifies a ggplot2 object
#' @keywords internal
add_datapoints <- function(gg, dat){
  gg +
    ggplot2::geom_linerange(data = dat,
                            ggplot2::aes(x = .data$depth,
                                         ymax = .data$age_upr2, ymin = .data$age_lwr2,
                                         colour = "Obs age"),
                            show.legend = FALSE,
                            group = NA,
                            inherit.aes = FALSE,
                            linewidth = 0.5) +
    ggplot2::geom_linerange(data = dat,
                            ggplot2::aes(x = .data$depth,
                                         ymax = .data$age_upr1, ymin = .data$age_lwr1,
                                         colour = "Obs age"),
                            show.legend = FALSE,
                            group = NA,
                            inherit.aes = FALSE,
                            linewidth = 1) +
    ggplot2::geom_point(
      data = dat, ggplot2::aes(
        y = .data$age,
        group = NA,
        colour = "Obs age")
    ) +
    ggplot2::labs(x = "Depth", y = "Age")
}

add_infl_err <- function(gg, hamstr_fit, obs_ages){
  
  infl_errs <- rstan::summary(hamstr_fit$fit, par = "obs_err_infl")$summary %>%
      tibble::as_tibble(rownames = "par") %>%
      dplyr::mutate(dat_idx = get_par_idx(.data$par))
    
  obs_ages <- obs_ages %>%
      dplyr::mutate(infl_err = infl_errs$mean,
                    age_lwr_infl1 = .data$age + 1*.data$infl_err,
                    age_upr_infl1 = .data$age - 1*.data$infl_err,
                    age_lwr_infl2 = .data$age + 2*.data$infl_err,
                    age_upr_infl2 = .data$age - 2*.data$infl_err)
    
  gg <- gg +
      ggplot2::geom_linerange(
        data = obs_ages,
        ggplot2::aes(x = .data$depth, ymax = .data$age_upr_infl2, ymin = .data$age_lwr_infl2,
                     colour = "Infl err"), show.legend = TRUE,
        group = NA,
        lwd = 0.5, inherit.aes = F)+
    ggplot2::geom_linerange(
      data = obs_ages,
      ggplot2::aes(x = .data$depth, ymax = .data$age_upr_infl1, ymin = .data$age_lwr_infl1,
                   colour = "Infl err"), show.legend = TRUE,
      group = NA,
      lwd = 1,      inherit.aes = F)
 
  return(gg)
}


#' Plot Age Models as Spaghetti Plot
#'
#' @inheritParams plot_hamstr
#'
#' @return A ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom rlang .data

#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error)
#'
#' plot_age_models(fit)
#' }
plot_age_models <- function(hamstr_fit, n.iter = 1000){
  obs_ages <- get_obs_age_err(hamstr_fit)

  if (hamstr_fit$data$sample_posterior == FALSE){

    gg <- ggplot2::ggplot(
      obs_ages,
      aes(x = .data$depth, y = .data$age, fill = NA)) +
      geom_blank() +
      theme_bw()

    gg <- add_datapoints(gg, obs_ages)
    gg <- add_subdivisions(gg, hamstr_fit)
    gg <- add_colour_scale(gg)

    return(gg)

  }

  posterior_ages <- get_posterior_ages(hamstr_fit)


  p.fit <- posterior_ages %>%
    dplyr::filter(.data$iter %in% sample(unique(.data$iter), n.iter, replace = FALSE)) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data$depth, y = .data$age, group = .data$iter, fill = NA
        )
      )

  p.fit <- p.fit +
    ggplot2::geom_line(alpha = 0.5 / sqrt(n.iter),
                       aes(colour = "Age models"),
                       key_glyph = "abline"
                       )


  if (hamstr_fit$data$inflate_errors == 1){

    p.fit <- add_infl_err(p.fit, hamstr_fit, obs_ages)
    
  }


  p.fit <- add_datapoints(p.fit, obs_ages)+
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::theme_bw()

  # add subdivisions
  p.fit <- add_subdivisions(p.fit, hamstr_fit)

  # add colour scales
  p.fit <- add_colour_scale(p.fit)

  return(p.fit)

}



## Accumulation rates ----



#' Plot Downcore Summary
#' @param ds a downcore summary of age or accumulation rate
#' @param axis units for the x axis, depth or age
#' @return a ggplot2 object
#' @keywords internal
plot_downcore_summary <- function(ds, axis = c("depth", "age")){

  axis <- match.arg(axis)

  p <- ds %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[axis]], y = .data[["mean"]]), show.legend = FALSE) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymax = .data$`2.5%`, ymin = .data$`97.5%`,
                   fill = "95%")) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymax = .data$`15.9%`, ymin = .data$`84.1%`,
                   fill = "68%")) +
    ggplot2::geom_line(aes(colour = "Mean"),
                       key_glyph = "abline") +
    ggplot2::geom_line(ggplot2::aes(y = .data$`50%`, colour = "Median"),
                       key_glyph = "abline") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())


  return(p)
}


#' Plot accumulation rates
#' @inheritParams plot_hamstr
#' @param axis plot accumulation rate against depth or age
#' @param units plot accumulation rate in depth per time, or time per depth
#' (or both)
#' @inheritParams filter_hamstr_acc_rates
#' @return a ggplot2 object
#' @keywords internal
plot_hamstr_acc_rates <- function(hamstr_fit,
                                  axis = c("depth", "age"),
                                  units = c("depth_per_time", "time_per_depth"),
                                  tau = 0, kern = c("U", "G", "BH"), ...){

  units <- match.arg(units,
                     #choices = c("depth_per_time", "time_per_depth"),
                     several.ok = TRUE)

  axis <- match.arg(axis,
                    #choices = c("depth_per_time", "time_per_depth"),
                    several.ok = TRUE)

  kern <- match.arg(kern)

  acc_rates <- summarise_hamstr_acc_rates(hamstr_fit, tau = tau, kern=kern,
                                          probs = c(0.025, 0.159, 0.25, 0.5, 0.75, 0.841, 0.975))

  if ("depth" %in% axis){
    acc_rates_long <- acc_rates %>%
      dplyr::select(-"depth") %>%
      tidyr::pivot_longer(cols = c("c_depth_top", "c_depth_bottom"),
                          names_to = "depth_type", values_to = "depth")

    rug_dat <- data.frame(d = hamstr_fit$data$depth)

    p.depth <- acc_rates_long %>%
      dplyr::filter(.data$acc_rate_unit %in% units) %>%
      plot_downcore_summary() +
      ggplot2::labs(x = "Depth", y = "Accumulation rate") +
      ggplot2::facet_wrap(~.data$acc_rate_unit, scales = "free_y") +
      ggplot2::geom_rug(data = rug_dat, aes(x = .data$d, colour = "Age point"),
                        inherit.aes = FALSE, key_glyph = draw_key_vline)
    
    p.depth <- add_colour_scale(p.depth) 

    if ("hamstr_fit" %in% class(hamstr_fit)){
      p.depth <- add_subdivisions(p.depth, hamstr_fit = hamstr_fit)
    }


  }

  if ("age" %in% axis){

    median_age <- summary(hamstr_fit) %>%
      #mutate(unit = "age") %>%
      dplyr::rename(age = "50%") %>%
      dplyr::select("depth", "age")

    jnt <- dplyr::left_join(median_age, acc_rates) %>%
      dplyr::filter(stats::complete.cases(.data$mean))



    p.age <- jnt %>%
      dplyr::filter(.data$acc_rate_unit %in% units) %>%
      plot_downcore_summary(axis = "age") +
      ggplot2::labs(x = "Age", y = "Accumulation rate") +
      ggplot2::facet_wrap(~.data$acc_rate_unit, scales = "free_y") #+
    
    p.age <- add_colour_scale(p.age) 
    
    
    if ("hamstr_fit" %in% class(hamstr_fit)){
      rug_dat <- data.frame(a = hamstr_fit$data$obs_age)
      p.age <- p.age +
        ggplot2::geom_rug(data = rug_dat, aes(x = .data$a, colour = "Age point"),
                          inherit.aes = FALSE, key_glyph = draw_key_vline)
    }
  }

  if (length(axis) == 2){
    #p <- ggpubr::ggarrange(p.depth, p.age, nrow = 2)
    p <- p.depth / p.age
  } else if (axis == "depth"){
    p <- p.depth
  } else {
    p <- p.age
  }

  p

}

seg_lvls_2_brk_lvls <- function(seg_lvls){
  
  lvls <- unique(seg_lvls)
  cnts <- table(seg_lvls)
  
  rep(lvls, times = cnts+1)

}



#' Plot the hierarchical accumulation rate parameters
#'
#' @inheritParams plot_hamstr
#'
#' @return ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error)
#'
#' plot_hierarchical_acc_rate(fit)
#' 
#' plot_hierarchical_acc_rate(hamstr_fit = fit) +
#'  facet_wrap(~lvl)
#' }
plot_hierarchical_acc_rate <- function(hamstr_fit){
  
  if (is.null(hamstr_fit$data$brks)) {
    
    gg <- plot_hierarchical_acc_rate(hamstr_fit)
    
    return(gg)
  } 
  
  brks <- hamstr_fit$data$brks
  dpth_rng <- hamstr_fit$data$bottom_depth - hamstr_fit$data$top_depth
  
  idx <- dplyr::tibble(
    brks = unlist(brks)
  ) %>% 
    dplyr::mutate(
      depth = brks * dpth_rng + hamstr_fit$data$top_depth
    )
  
  a3 <- rstan::summary(hamstr_fit$fit, pars = "alpha")$summary
  
  alph <- tibble::as_tibble(a3, rownames = "par") %>%
    dplyr::mutate(alpha_idx = get_par_idx(.data$par),
                  lvl = hamstr_fit$data$lvl) %>%
    dplyr::select(lvl, mean) %>% 
    dplyr::group_by(.data$lvl) %>% 
    dplyr::reframe(mean = c(mean, utils::tail(mean, 1)))
  
  out <- dplyr::bind_cols(idx, alph) %>% 
    dplyr::group_by(.data$lvl) %>% 
    dplyr::mutate(depth = dplyr::case_when(.data$depth < hamstr_fit$data$top_depth ~ hamstr_fit$data$top_depth,
                                           .data$depth > hamstr_fit$data$bottom_depth ~ hamstr_fit$data$bottom_depth,
                                           TRUE ~ .data$depth)) %>%
    dplyr::ungroup()
  
  out <- out %>% 
    ggplot2::ggplot(aes(x=.data$depth, y=mean, colour = factor(.data$lvl-1))) +
    ggplot2::geom_step() +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(y = "Acc. rate [age/depth]", x = "Depth",
                  colour = "Hierarchical\nlevel") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "top") +
    ggplot2::guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
    ggplot2::scale_colour_brewer(palette = "Dark2")
  
  return(out)
}

## Prior and posteriors --------

#' Plot a Prior and Posterior
#'
#' @param prior Numerical PDF of prior
#' @param posterior Sample from posterior distribution
#' @return a ggplot2 object
#' @keywords internal
#' @import ggplot2
#' @importFrom stats density
plot_prior_posterior_hist <- function(prior, posterior, bins = 50){
  clrs <- c("Posterior" = "#5785C1", "Prior" = "#CB7A5C")

  gg <- ggplot2::ggplot(data = NULL, aes(fill = "Posterior"))

  if (is.null(posterior) == FALSE){

    bw <- stats::sd(posterior$x) / 2

    gg <- gg + ggplot2::geom_histogram(data = posterior,
                   ggplot2::aes(x = .data$x, ggplot2::after_stat(density),
                   fill = "Posterior"),
                   colour = NA,
                   alpha = 0.75, binwidth = bw)
  }
    gg +
    ggplot2::geom_line(data = prior, ggplot2::aes(x = .data$x, y = .data$d,
                                                  colour = "Prior")) +
    ggplot2::facet_wrap(~.data$par, scales = "free_y", ncol = 1) +
    ggplot2::scale_colour_manual("", values = clrs[2]) +
    ggplot2::scale_fill_manual("", values = clrs[1]) +
    ggplot2::labs(
      x = "Value",
      y = "Density",
      colour = "",
      fill = ""
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 2), 
                      shape = ggplot2::guide_legend(order = 1)) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::theme_bw()
}


#' Plot the Prior and Posterior Distributions of the Inflation Factor Parameters
#'
#' @return a ggplot2 object
#' @import rstan
#' @import ggplot2

#' @importFrom rlang .data
#' @inheritParams plot_hamstr
#' @keywords internal
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error,
#'   inflate_errors = 1)
#'
#' plot_infl_prior_posterior(fit)
#' }
plot_infl_prior_posterior <- function(hamstr_fit) {
  clrs <- c("Posterior" = "Blue", "Prior" = "Red")

  hamstr_dat <- hamstr_fit$data

  infl_mean_shape_post <-
    tibble::tibble(
      infl_mean = as.vector(rstan::extract(hamstr_fit$fit, "infl_mean")[[1]]),
      infl_shape = as.vector(rstan::extract(hamstr_fit$fit, "infl_shape")[[1]])
    ) 
  
  infl_mean_shape_post <- infl_mean_shape_post %>%
    dplyr::mutate(iter = 1:nrow(infl_mean_shape_post))


  max_x_shape <- with(hamstr_dat, {
    infl_shape_prior_upr <- stats::qgamma(c(0.99),
                                          shape = infl_shape_shape,
                                          rate = infl_shape_shape / infl_shape_mean)

    max(c(infl_shape_prior_upr, infl_mean_shape_post$infl_shape))
  })

  max_x_mean <- with(hamstr_dat, {
    infl_mean_prior_upr <- stats::qnorm(c(0.99), 0, infl_sigma_sd)
    max(c(infl_mean_prior_upr, infl_mean_shape_post$infl_mean))
  })


  infl_fac <- rstan::extract(hamstr_fit$fit, "infl")[[1]] %>%
    tibble::as_tibble(.data$., .name_repair = c("unique")) %>%
    tidyr::gather() %>%
    dplyr::mutate(key = get_par_idx(.data$key))


  p.infl.fac <- rstan::stan_plot(hamstr_fit$fit, pars = "infl")



  infl_prior_shape <-
    tibble::tibble(x = seq(0, max_x_shape, length.out = 1000)) %>%
    dplyr::mutate(
      d = stats::dgamma(.data$x - 1,
                        hamstr_dat$infl_shape_shape,
                        rate = hamstr_dat$infl_shape_shape / hamstr_dat$infl_shape_mean),
      par = "infl_shape"
    )

  infl_prior_mean <-
    tibble::tibble(x = seq(0, max_x_mean, length.out = 1000)) %>%
    dplyr::mutate(
      d = 2 * stats::dnorm(.data$x, 0, sd = hamstr_dat$infl_sigma_sd),
      par = "infl_mean"
    )

  infl_priors <- dplyr::bind_rows(infl_prior_mean, infl_prior_shape)

  infl_mean_shape_post_long <- infl_mean_shape_post %>%
    tidyr::pivot_longer(cols = c("infl_mean", "infl_shape"), names_to = "par", values_to = "x")

  p.pars <- plot_prior_posterior_hist(infl_priors, infl_mean_shape_post_long, bins = 1000) +
    facet_wrap(~.data$par, scales = "free")


  infl_mean_shape_post <- infl_mean_shape_post %>%
    dplyr::mutate(q99 = stats::qgamma(0.75, .data$infl_shape,
      rate = .data$infl_shape / .data$infl_mean
    ))


  infl_pars_prior_dist <- infl_mean_shape_post %>%
    dplyr::filter(.data$iter %in% sample.int(dplyr::n(), 10)) %>%
    tidyr::crossing(tibble::tibble(x = exp(seq(
      log(0.01), log(stats::quantile(infl_mean_shape_post$q99, prob = 0.95)),  length.out = 100
    )))) %>%
    dplyr::mutate(
      d = stats::dgamma(.data$x, shape = .data$infl_shape, rate = .data$infl_shape / .data$infl_mean),
      par = "Modelled prior for infl_fac"
    )


  p.priors <- infl_pars_prior_dist %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$d,
                                 group = .data$iter)) +
    ggplot2::geom_line(
      alpha = 1, # /sqrt(100),
      aes(colour = "Infl err")
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(y = "Density", x = "Value")

 
  p <- p.pars / (p.priors | p.infl.fac)

  return(p)
}


#' Plot Mean Accumulation Rate Prior and Posterior Distributions
#' @inheritParams plot_hamstr
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @return a ggplot2 object
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
#'   inflate_errors = 0)
#'
#' plot_acc_mean_prior_posterior(fit)
#' }
plot_acc_mean_prior_posterior <- function(hamstr_fit) {

  prior_mean <- hamstr_fit$data$acc_mean_prior

  acc_prior_rng <- stats::qnorm(c(0.999), mean = 0, sd = 10 * prior_mean)

  prior <-  tibble::tibble(
    x = seq(0, acc_prior_rng[1], length.out = 1000)
    ) %>%

    dplyr::mutate(
      par = "acc_mean",
      d = 2 * stats::dnorm(.data$x, 0, 10 * prior_mean)
    )

  prior$d[prior$x <= 0] <- 0


  if (hamstr_fit$data$sample_posterior == TRUE){
    post <-
      tibble::tibble(x = as.vector(rstan::extract(hamstr_fit$fit, "alpha[1]")[[1]]))
  } else {
    post <- NULL
  }

  p <- plot_prior_posterior_hist(prior, post) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    labs(x = "Mean accumulation rate [age/depth]")


  return(p)

}


#' Plot Memory Prior and Posterior
#'
#' @inheritParams plot_hamstr
#'
#' @return a ggplot2 object
#' @import ggplot2
#' @importFrom rlang .data
#' @keywords internal
#' @examples
#' \dontrun{
#' fit <- hamstr(
#'   depth = MSB2K$depth,
#'   obs_age = MSB2K$age,
#'   obs_err = MSB2K$error)
#'
#' plot_memory_prior_posterior(fit)
#' }
plot_memory_prior_posterior <- function(hamstr_fit){
  # memory prior
  mem.prior <- tibble::tibble(x = seq(0, 1, length.out = 1000)) %>%
    dplyr::mutate(d = stats::dbeta(.data$x,
                                   shape1 = hamstr_fit$data$mem_alpha,
                                   shape2 = hamstr_fit$data$mem_beta),
                  par = "Memory at 1 depth unit")


  if (hamstr_fit$data$sample_posterior == TRUE){

    w <- rstan::extract(hamstr_fit$fit, "w")$w

    ifelse(is.matrix(w), w <- apply(w, 1, stats::median),  w <- as.vector(w))

     mem.post <- tibble::tibble(
       w = w,
       R = as.vector(rstan::extract(hamstr_fit$fit, "R")$R)
       ) %>%
    dplyr::rename(`Memory between sections` = "w",
                  `Memory at 1 depth unit` = "R") %>%
    tidyr::pivot_longer(cols = c("Memory between sections", "Memory at 1 depth unit"),
                        names_to = "par", values_to = "x")
  } else {
    mem.post <- NULL
  }

  p.mem <- plot_prior_posterior_hist(mem.prior, mem.post) +
    expand_limits(x = c(0, 1)) +
    theme(legend.position = "top")


  return(p.mem)
}


#' Plot Mixing Depth Prior and Posterior
#'
#' @inheritParams plot_hamstr
#'
#' @return A ggplot2 object
#' @import ggplot2
#' @importFrom rlang .data
#' @keywords internal
plot_L_prior_posterior <- function(hamstr_fit){

  if (hamstr_fit$data$L_prior_shape > 0){
    post <- as.data.frame(hamstr_fit$fit, pars = c("L")) %>%
    dplyr::as_tibble() %>%
    tidyr::pivot_longer(cols = dplyr::everything(),
                        names_to = "par", values_to = "x") %>%
    dplyr::mutate(dpt = get_par_idx(.data$par),
           par = "L")


  L_shp <- hamstr_fit$data$L_prior_shape
  L_mean <- hamstr_fit$data$L_prior_mean

  L_prior_rng <- stats::qgamma(c(0, 0.999),
                               shape = L_shp,
                               scale = L_mean / L_shp)

  L_prior_rng[1] <- min(c(L_prior_rng[1]), post$x)
  L_prior_rng[2] <- max(c(L_prior_rng[2]), post$x)

  L_prior <-  tibble::tibble(
    x = seq(L_prior_rng[1], L_prior_rng[2], length.out = 1000)
    ) %>%
    dplyr::mutate(
      par = "L",
      d = stats::dgamma(
        .data$x,
        shape = L_shp,
        scale = L_mean / L_shp)
      )

  plot_prior_posterior_hist(L_prior, post)+
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "Mixing depth [L]")
  } else {

    ggplot2::ggplot(data = tibble::tibble(x = hamstr_fit$data$L_prior_mean * c(1, 1), y = c(0, 1))) +
      ggplot2::geom_line( aes(x = .data$x , y = .data$y, colour = "Fixed")) +
      ggplot2::expand_limits(x = c(0, 2*hamstr_fit$data$L_prior_mean))+
      ggplot2::labs(x = "Mixing depth [L]", y = "Density") +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_manual("", values = c(Fixed = "Red"))

  }

  }


#' Plot Displacement Depth Prior and Posterior
#'
#' @inheritParams plot_hamstr
#'
#' @return A ggplot2 object
#' @import ggplot2
#' @importFrom rlang .data
#' @keywords internal
plot_D_prior_posterior <- function(hamstr_fit) {

  prior_mean <- hamstr_fit$data$D_prior_scale

  prior_rng <- stats::qnorm(c(0.999), mean = 0, sd = prior_mean)

  prior <-  tibble::tibble(
    x = seq(0, prior_rng[1], length.out = 1000)
  ) %>%

    dplyr::mutate(
      par = "D",
      d = 2 * stats::dnorm(.data$x, 0, prior_mean)
    )

  prior$d[prior$x <= 0] <- 0

  if (hamstr_fit$data$sample_posterior == TRUE){
    post <-
      tibble::tibble(x = as.vector(rstan::extract(hamstr_fit$fit, "D[1]")[[1]]))
  } else {
    post <- NULL
  }

  p <- plot_prior_posterior_hist(prior, post) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    labs(x = "Displacement [D]", y = "Density")

  return(p)

}




plot_14C_PDF <- function(hamstr_fit, nu = 6, cal_curve) {

  compare_14C_PDF(age.14C = hamstr_fit$data$obs_age,
                  age.14C.se = hamstr_fit$data$obs_err,
                  cal_curve = cal_curve,
                  nu = nu)

}


#' Add subdivision tickmarks
#'
#' @param gg A ggplot2 object
#' @inheritParams plot_hamstr
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @keywords internal
add_subdivisions <- function(gg, hamstr_fit) {

  tick_dat <- hierarchical_depths(hamstr_fit$data)

  for (i in seq_along(tick_dat)) {

    df <- data.frame(x = tick_dat[[i]])

    lnth <- length(tick_dat) - (i - 1)

    gg <- gg + ggplot2::geom_rug(
      data = df,
      ggplot2::aes(x = .data$x),
      inherit.aes = F,
      sides = "top",
      length = ggplot2::unit(0.01 * lnth, "npc")
      )

  }

  return(gg)

}


