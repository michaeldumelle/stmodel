% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initial.R
\name{initial}
\alias{initial}
\title{Create Initial Values}
\usage{
initial(
  s_de,
  s_ie,
  t_de,
  t_ie,
  st_de,
  st_ie,
  v_s,
  v_t,
  s_range,
  t_range,
  estmethod,
  stcov
)
}
\arguments{
\item{s_de}{The spatial dependent variance (spatial partial sill).}

\item{s_ie}{The spatial independent variance (spatial nugget).}

\item{t_de}{The temporal dependent variance (temporal partial sill).}

\item{t_ie}{The temporal independent variance (temporal nugget).}

\item{st_de}{The spatio-temporal dependent variance (spatio-temporal partial sill).}

\item{st_ie}{The spatio-temporal independent variance (spatio-temporal nugget).}

\item{v_s}{The proportion of spatial dependent variance
(if \code{estmethod = "product"}).}

\item{v_t}{The proportion of temporal dependent variance
(if \code{estmethod = "product"}).}

\item{s_range}{The spatial effective range (the spatial distance at which
the correlation equals 0.05 (for non-compact) or 0 (for compact correlations)}

\item{t_range}{The spatial effective range (the spatial distance at which
the correlation equals 0.05 (for non-compact) or 0 (for compact correlations)}

\item{estmethod}{The estimation method
\describe{
\item{\code{reml}}{Restricted Maximum Likelihood}
\item{\code{svwls}}{Semivariogram Weighted Least Squares}
}}

\item{stcov}{The spatio-temporal covariance type
\describe{
\item{\code{product}}{The product LMM}
\item{\code{sum_with_error}}{The sum-with-error LMM}
\item{\code{productsum}}{The product sum LMM}
}}
}
\value{
A named vector with covariance parameters having class equal to
the \code{estmethod} argument and the \code{stcov} argument.
}
\description{
A wrapper around \code{make_covparam_object} that requires specification
of \code{estmethod}.
}
