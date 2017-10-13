plot_prior_posterior <- function(fit){
  fit$data
  
  acc.rng <- qgamma(c(0.000001, 0.999), shape = fit$data$acc_alpha,  rate = fit$data$acc_beta)
  
  acc.prior <- tibble(acc.rate = seq(acc.rng[1], acc.rng[2], length.out = 1000)) %>% 
    mutate(acc.dens = dgamma(acc.rate, shape = fit$data$acc_alpha,  rate = fit$data$acc_beta))
  
  acc.post <- tibble(alpha = as.vector(rstan::extract(fit$fit, "alpha")$alpha))
  
  acc.prior %>% 
    ggplot(aes(x = acc.rate, y = acc.dens)) +
    geom_density(data = acc.post, aes(x = alpha), fill = "Grey", inherit.aes = FALSE) +
    geom_line(colour = "Red") + 
    scale_x_continuous(limits = acc.rng)
    
  
}
