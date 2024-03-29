% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main functions.R
\name{supclus}
\alias{supclus}
\title{SupClus explanations}
\usage{
supclus(k, x, y, h, dmax)
}
\arguments{
\item{k}{Number of interpretable clusters.}

\item{x}{Matrix \eqn{n \times d} with covariates for each sample unit.}

\item{y}{\eqn{n}-dimensional vector with predictions from the original model.}

\item{h}{Bandwidth of gaussian kernel.}

\item{dmax}{Maximum number of covariates (\eqn{\le d}).}
}
\value{
Matrix \eqn{n \times dmax} with local coefficients and \eqn{n}-dimensional vector with cluster of each instance.
}
\description{
This function calculates the local coefficients using the SubClus method for a sample.
}
\examples{
library(interpretability)

data("artificial_data")

supcluscoef <- supclus(k = 2,
                       x = artificial_data$x_test,
                       y = artificial_data$y_pred,
                       h = 0.1,
                       dmax = ncol(artificial_data$x_test))
}
