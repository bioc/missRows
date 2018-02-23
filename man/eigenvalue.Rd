\name{eigenvalue}

\alias{eigenvalue}

\title{
Firts Singular Value of the SVD of a Matrix
}
\description{
Obtain the firts singular value of the singular-value decomposition
of a rectangular matrix as computed by \code{\link{svd}}. This function
is called internally by \code{\link{MFA}} and is not usually called
directly by a user.
}

\usage{
eigenvalue(x)
}

\arguments{
  \item{x}{a numeric or complex matrix whose SVD decomposition
        can be computed.}
}

\value{
The firts singular value of \code{svd(x)}.
}

\keyword{ kwd1 }
\keyword{ kwd2 }
