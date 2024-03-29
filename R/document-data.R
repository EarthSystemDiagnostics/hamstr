#' @title MSB2K
#' @description Dataset for core MSB2K as provided in the package \link[rbacon]{rbacon}
#' @format A data frame with 40 rows and 4 variables:
#' \describe{
#'   \item{\code{labID}}{factor Sample ID number from the lab.}
#'   \item{\code{age}}{integer Radiocarbon age.}
#'   \item{\code{error}}{integer Age uncertainty in years.}
#'   \item{\code{depth}}{double Depth in core in cm.}
#'}
#' @details The example dataset MSB2K as provided by the \link[rbacon]{rbacon} package.
#' @references
#'
#' Blaauw, Maarten, and J. Andrés Christen. 2011. Flexible
#'   Paleoclimate Age-Depth Models Using an Autoregressive Gamma Process.
#'   Bayesian Analysis 6 (3): 457-74. \doi{10.1214/ba/1339616472}
#'
#' Maarten Blaauw, J. Andres Christen and Marco A. Aquino L. (2020). rbacon:
#'   Age-Depth Modelling using Bayesian Statistics. R package version 2.4.3.
#'   <https://CRAN.R-project.org/package=rbacon>
"MSB2K"
