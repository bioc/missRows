\name{RVcoeff}

\alias{RVcoeff}

\title{
Calculate the RV Coefficient
}

\description{
Calculates the RV coefficient between two matrice 
\code{X} and \code{Y}. This function is called 
internally by \code{\link{tuneM}} and is not usually 
called directly by a user.
}

\usage{
RVcoeff(X, Y)
}

\arguments{
  \item{X}{a matrix with \eqn{n} rows (individuals) and 
        \eqn{p} columns (variables).}
        
  \item{Y}{a matrix with \eqn{n} rows (individuals) and 
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
