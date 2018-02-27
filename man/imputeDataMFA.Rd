\name{imputeDataMFA}

\encoding{latin1}
\alias{imputeDataMFA}

\title{
Impute Missing Rows of Data Tables
}

\description{
Impute the missing rows of data tables using the alternating least 
squares algorithm used in PCA. This function is called internally by \code{\link{MIMFA}} and is not usually called directly by a user.
}

\usage{
imputeDataMFA(datasets, U, missRows, comp, max.iter = 500,
              tol = 1e-10)
}

\arguments{

  \item{datasets}{a list containing the data tables with missing
        rows. Tables in the list should be arranged in samples
        \eqn{\times} variables, with samples order matching in all
        data tables.}

  \item{U}{the compromise configuration, a matrix with the individuals
        coordinates as returned by \code{\link{wrapperSTATIS}} function.}

  \item{missRows}{a character vector containing the name of
        the missing rows.}

  \item{comp}{number of components kept for imputation.}

  \item{max.iter}{integer, maximum number of iterations for the 
        iterative algorithm.}

 \item{tol}{positive value, the threshold for assessing convergence.}
}

\details{
Since the core of MFA is a PCA of the weighted data table \eqn{K}, 
the algorithm suggested to estimate MFA axes and impute missing 
values is inspired from the alternating least squares algorithm 
used in PCA. This consists in finding matrices \eqn{F} and \eqn{U} 
which minimize the following criterion:

\eqn{||K-M-FU||^2 = \sum_{i}\sum_{k}\left( K_{ik} - M_{ik} - \sum_{d=1}^D F_{id} U_{kd} \right)^2,}

where \eqn{M} is a matrix with each row equal to a vector of the 
mean of each variable and \eqn{D} is the kept dimensions in PCA. 
The solution is obtained by alternating two multiple regressions 
until convergence, one for estimating axes (loadings \eqn{\hat{U}}) 
and one for components (scores \eqn{\hat{F}}):

  \eqn{
    \hat{U}' = (\hat{F}'\hat{F})^{-1}\hat{F}'(K-\hat{M})\\
    \hat{F} = (K-\hat{M})\hat{U}(\hat{U}'\hat{U})^{-1}.
  }
}

\value{
An list containing components with the imputed rows for
each data table.
}

\author{
Ignacio González
}

\keyword{multivariate}
