\name{eigenvalue}

\alias{eigenvalue}

\title{
Scaled Firts Singular Value of the SVD of a Matrix
}

\description{
Obtain the scaled first singular value of the singular-value decomposition
of a rectangular matrix \code{X} as computed by \code{\link{svd}}. Scaling 
is done by dividing the firts singular value by the root square of the number 
of rows in \code{X}. This function is internally called by \code{\link{MFA}} 
and is not usually called directly by a user.
}

\usage{
eigenvalue(X)
}

\arguments{
    \item{X}{a numeric matrix whose SVD decomposition
            can be computed.}
}

\value{
The scaled first singular value of \code{svd(X)}.
}

\keyword{algebra}
