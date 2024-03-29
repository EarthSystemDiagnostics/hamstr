% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_stan_hamstr.R
\name{plot_hamstr_acc_rates}
\alias{plot_hamstr_acc_rates}
\title{Plot accumulation rates}
\usage{
plot_hamstr_acc_rates(
  hamstr_fit,
  axis = c("depth", "age"),
  units = c("depth_per_time", "time_per_depth"),
  tau = 0,
  kern = c("U", "G", "BH"),
  ...
)
}
\arguments{
\item{hamstr_fit}{the object returned from \code{stan_hamstr}.}

\item{axis}{plot accumulation rate against depth or age}

\item{units}{plot accumulation rate in depth per time, or time per depth
(or both)}

\item{tau}{scale of the smoothing kernel in depth units. If tau > 0,
accumulation rates are smoothed (filtered) before summary statistics are
calculated, so that the accumulation rate at a given depth corresponds to the
average rate over the depth interval tau. Default to 0.}

\item{kern}{choice of smoothin kernal. U for uniform (moving average), G for
Gaussian, BH for Berger and Heath (exponential). Defaults to U}
}
\value{
a ggplot2 object
}
\description{
Plot accumulation rates
}
\keyword{internal}
