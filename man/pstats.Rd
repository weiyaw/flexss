% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/semi2.R
\name{pstats}
\alias{pstats}
\alias{pstats_v4}
\alias{pmean_v4}
\title{Get posterior statistics}
\usage{
pstats_v4(samples, fun)

pmean_v4(samples)
}
\arguments{
\item{samples}{A vector, matrix or array of samples}

\item{fun}{A function for statitics calculation. Support `tidy` style of
function specification.}
}
\value{
A numeric value, vector or matrix
}
\description{
Calculate posterior statistics from an array or vector of samples. The samples are
populated along the last dimension, e.g. columns of a matrix are the samples;
the depth of a 3D array are the samples.
}
\details{
`pstats` calculates posterior statistics given by `fun`.

`pmean` calculates posterior means, and faster than pstats(samples, mean).
}
