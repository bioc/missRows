\name{wrapperSVD}

\alias{wrapperSVD}

\title{
Singular Value Decomposition of a Matrix
}

\description{
Wrapper function from the \pkg{FactoMineR} package to perform
Singular Value Decomposition of a matrix. This function
is called internally by \code{\link{MFA}} and is not usually 
called directly by a user.
}

\usage{
wrapperSVD(X, rWeights=NULL, cWeights=NULL, ncp=Inf)
}

\arguments{
    \item{X}{a data matrix.}

    \item{rWeights}{vector with the weights of each row 
        (\code{NULL} by default and the weights are uniform).}

    \item{cWeights}{vector with the weights of each column 
        (\code{NULL} by default and the weights are uniform).}

    \item{ncp}{the number of components kept for the outputs.}
}

\value{
A list that contains the following components:

    \item{vs}{a vector containing the singular values of \code{X}.}

    \item{U}{a matrix whose columns contain the left singular 
        vectors of \code{X}.}

    \item{V}{a matrix whose columns contain the right singular 
        vectors of \code{X}.}
}

\seealso{
\code{\link{svd}}, \code{\link[FactoMineR]{svd.triplet}}
}

\examples{
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
X <- hilbert(9)[, 1:6]

wrapperSVD(X)
}

\keyword{algebra}
