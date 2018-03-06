\name{RVcoeff}

\alias{RVcoeff}

\title{
Calculate the RV Coefficient
}

\description{
Calculate the RV coefficient between two matrices 
\code{X} and \code{Y}. This function is internally called by \code{\link{tuneM}} and is not usually 
called directly by a user.
}

\usage{
RVcoeff(X, Y)
}

\arguments{
  \item{X, Y}{a matrix with \eqn{n} rows (individuals) and 
        \eqn{p} columns (variables).}
}

\value{
The RV coefficient between the two matrices.
}

\references{
Robert, P., Escoufier, Y. (1976). 
A unifying tool for linear multivariate statistical methods: 
The RV coefficient. 
\emph{Journal of the Royal Statistical Society}.
25(3), 257--265.
}

\keyword{multivariate}
