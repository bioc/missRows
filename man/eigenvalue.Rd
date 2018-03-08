\name{eigenvalue}

\alias{eigenvalue}

\title{
Scaled Firts Singular Value of the SVD of a Matrix
}

\description{
Obtain the scaled first singular value of the singular-value decomposition
of a rectangular matrix \code{x} as computed by \code{\link{svd}}. Scaling 
is done by dividing the firts singular value by the root square of the number 
of rows in \code{x}. This function is internally called by \code{\link{MFA}} 
and is not usually called directly by a user.
}

\usage{
eigenvalue(x)
}

\arguments{
    \item{x}{a numeric or complex matrix whose SVD decomposition
            can be computed.}
}

\value{
The scaled first singular value of \code{svd(x)}.
}

\examples{
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
x <- hilbert(9)[, 1:6]
eigenvalue(x)
}

\keyword{algebra}
