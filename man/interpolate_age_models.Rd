% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_methods.R
\name{interpolate_age_models}
\alias{interpolate_age_models}
\title{Interpolate Posterior Age Model At New Depths}
\usage{
interpolate_age_models(hamstr_fit, depth)
}
\arguments{
\item{hamstr_fit}{the object returned from \code{stan_hamstr}.}

\item{depth}{a vector of depths at which to interpolate the age models.
If left NULL, the depths of the age control points are used.}
}
\value{
hamstr_interpolated_ages object
}
\description{
Interpolate Posterior Age Model At New Depths
}
\examples{
\dontrun{
fit <- hamstr(
  depth = MSB2K$depth,
  obs_age = MSB2K$age,
  obs_err = MSB2K$error)

interpolate.age.models(fit, depth = seq(1000, 15000, by = 1000))
}

}
\keyword{internal}
