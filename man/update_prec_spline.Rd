% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/semi2.R
\name{update_prec_spline}
\alias{update_prec_spline}
\title{Update precision of the spline terms}
\usage{
update_prec_spline(coefs, prior)
}
\arguments{
\item{coefs}{Matrix, each column the coef of (a subject|a level in a
factor). `attr(coefs, 'penalty')`: a full-row rank penalty matrix (for
ordinary splines). Kmat %*% coefs are penalised (i.e. are
random). `attr(coefs, 'block_dim')`: the dimension corresponding to the
block cov matrix (for subject splines)}
}
\value{
precision matrix for the coefs
}
\description{
Update precision of the spline terms
}
