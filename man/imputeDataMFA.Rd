\name{imputeDataMFA}

\encoding{latin1}
\alias{imputeDataMFA}

\title{
Impute Missing Rows and Estimates MFA Axes
}

\description{
Impute the missing rows of data tables using the alternating least 
squares algorithm used in PCA. This function is internally called by 
\code{\link{MIMFA}} and is not usually called directly by a user.
}

\usage{
imputeDataMFA(datasets, U, missRows, comp, maxIter=500, tol=1e-10)
}

\arguments{

    \item{datasets}{a list containing the data tables with missing
        rows. Tables in the list should be arranged in samples
        \eqn{\times}{x} variables, with samples order matching in all
        data tables.}

    \item{U}{the compromise configuration, a matrix with the individuals
        coordinates as returned by \code{\link{STATIS}} function.}

    \item{missRows}{a list containing character vectors with 
        the name of the missing individuals (rows) per table.}

    \item{comp}{a number of components kept for imputation.}

    \item{maxIter}{integer, maximum number of iterations for the 
        iterative algorithm.}

    \item{tol}{positive value, the threshold for assessing convergence.}
}

\details{
Since the core of MFA is a PCA of the merged data tables \eqn{K}, 
the algorithm suggested to estimate MFA axes and impute missing 
values is inspired from the alternating least squares algorithm 
used in PCA. This consists in finding matrices \eqn{F} and \eqn{U} 
which minimize the following criterion:

\deqn{
||K-M-FU||^2 = \sum_{i}\sum_{k}\left( K_{ik} - M_{ik} - 
                    \sum_{d=1}^D F_{id} U_{kd} \right)^2,}{
||K-M-FU||^2 = \sum_{i}\sum_{k}( K_{ik} - M_{ik} - 
                    \sum_{d=1}^D F_{id} U_{kd})^2,}

where \eqn{M} is a matrix with each row equal to a vector of the 
mean of each variable and \eqn{D} is the kept dimensions in PCA. 
The solution is obtained by alternating two multiple regressions 
until convergence, one for estimating axes (loadings \eqn{\hat{U}}) 
and one for components (scores \eqn{\hat{F}}):

    \deqn{
        \hat{U}' = (\hat{F}'\hat{F})^{-1}\hat{F}'(K - \hat{M})
    }
    \deqn{
        \hat{F} = (K - \hat{M})\hat{U}(\hat{U}'\hat{U})^{-1}.
    }

The \code{imputeDataMFA} algorithm first consists in imputing missing 
values in \eqn{K} with initial values (the column means on the 
non-missing entries), then \eqn{\hat{M}} is computed. The second step 
of the iterative algorithm is to calculate 
\eqn{\hat{F} = (K - \hat{M})U(U'U)^{-1}} on the completed dataset by 
using \eqn{D} components of \eqn{U}. Missing values are estimated as 
\eqn{\hat{K} = \hat{M} + \hat{F}U'}. The new imputed data set \eqn{K} 
is obtained by replacing the missing values of the original \eqn{K} 
matrix with the corresponding elements of \eqn{\hat{K}}, whilst keeping 
the observed values unaltered. These steps of estimation of the 
parameters and imputation of the missing values are iterate until 
convergence. The number \eqn{D} of components used in the algorithm can 
be estimated setting the \code{estim.ncp} argument to \code{TRUE} in the 
function \code{\link{MIMFA}}.
}

\value{
A list containing components with the imputed rows for
each data table.
}

\author{
Ignacio González
}

\keyword{multivariate}
