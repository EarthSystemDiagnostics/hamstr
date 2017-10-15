#' #' @title Plot
#' #' @description FUNCTION_DESCRIPTION
#' #' @param fit PARAM_DESCRIPTION
#' #' @return OUTPUT_DESCRIPTION
#' #' @details DETAILS
#' #' @examples 
#' #' \dontrun{
#' #' if(interactive()){
#' #'  #EXAMPLE1
#' #'  }
#' #' }
#' #' @seealso 
#' #'  \code{\link[ggpubr]{ggarrange}}
#' 
#' #'  \code{\link[rstan]{extract}}
#' #' @rdname plot_prior_posterior
#' #' @export 
#' #' @importFrom ggpubr ggarrange
#' #' @importFrom rstan extract
#' plot_prior_posterior <- function(fit){
#'   acc.rng <- qgamma(c(0.000001, 0.999), shape = fit$data$acc_alpha,  rate = fit$data$acc_beta)
#'   
#'   acc.prior <- tibble(acc.rate = seq(acc.rng[1], acc.rng[2], length.out = 1000)) %>% 
#'     mutate(acc.dens = dgamma(acc.rate, shape = fit$data$acc_alpha,  rate = fit$data$acc_beta))
#'   
#'   acc.post <- tibble(alpha = as.vector(rstan::extract(fit$fit, "alpha")$alpha))
#'   
#'   p.acc <- acc.prior %>% 
#'     ggplot(aes(x = acc.rate, y = acc.dens)) +
#'     geom_density(data = acc.post, aes(x = alpha), fill = "Grey", inherit.aes = FALSE) +
#'     geom_line(colour = "Red") + 
#'     scale_x_continuous("Acc. rate [yr/cm]", limits = acc.rng) +
#'     scale_y_continuous("Density") +
#'     theme_bw()
#'   
#'   # memory prior
#'   mem.prior <- tibble(mem = seq(0, 1, length.out = 1000)) %>% 
#'     mutate(mem.dens = dbeta(mem, shape1 = fit$data$mem_alpha,  shape2 = fit$data$mem_beta))
#'   
#'   mem.post <- tibble(R = as.vector(rstan::extract(fit$fit, "R")$R),
#'                      w = as.vector(rstan::extract(fit$fit, "w")$w))
#'   
#'   p.mem <- mem.prior %>% 
#'     ggplot(aes(x = mem, y = mem.dens)) +
#'     geom_density(data = mem.post, aes(x = R), fill = "Grey", inherit.aes = FALSE) +
#'     #geom_density(data = mem.post, aes(x = w), fill = "Grey", inherit.aes = FALSE) +
#'     geom_line(colour = "Red") + 
#'     scale_x_continuous("Memory [correlation at 1 cm]", limits = c(0, 1)) + 
#'     scale_y_continuous("") +
#'     theme_bw()
#'   
#'   p.fit <- plot_stan_bacon(fit)
#'   
#'     ggpubr::ggarrange(
#'     ggpubr::ggarrange(p.acc, p.mem, ncol = 2),
#'     p.fit,
#'     nrow = 2, heights = c(1, 2)) 
#' }



