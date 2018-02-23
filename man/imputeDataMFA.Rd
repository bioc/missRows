\name{imputeDataMFA}

\encoding{latin1}
\alias{imputeDataMFA}

\title{
imputeDataMFA
}

\description{
This function is called internally by \code{\link{MIMFA}} and
is not usually called directly by a user.
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

  \item{U}{the compromise configuration, a data
        frame or matrix?? with the individuals coordinates as returned by
        \code{\link{wrapperSTATIS}} function.}

  \item{missRows}{a character vector containing the name of
        the missing rows.}

  \item{comp}{number of components kept for imputation.}

  \item{max.iter}{Describe \code{max.iter} here}

  \item{tol}{Describe \code{tol} here}
}

\details{
 If necessary, more details than the description above
}

\value{
An list containing components with the imputed rows for
each data table.
}

\author{
Ignacio González
}

\keyword{ kwd1 }
\keyword{ kwd2 }
