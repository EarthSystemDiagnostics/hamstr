#' The 'hamstr' package.
#'
#' @description hamstr implements a "Bacon-like" Bayesian sediment accumulation
#'   or age-depth model with hierarchically structured multi-resolution discrete
#'   sediment sections. The Bayesian model is implemented in the Stan
#'   probabilistic programming language <https://mc-stan.org/>. "Blaauw,
#'   Maarten, and J. Andrés Christen. 2011. 'Flexible Paleoclimate Age-Depth
#'   Models Using an Autoregressive Gamma Process.'  Bayesian Analysis 6 (3):
#'   457-74. \doi{doi:10.1214/ba/1339616472}
#'
#' @docType package
#' @name hamstr-package
#' @useDynLib hamstr, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @import rstantools
#' @importFrom RcppParallel RcppParallelLibs CxxFlags
#'
#' @references
#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. <https://mc-stan.org>
#'
NULL
