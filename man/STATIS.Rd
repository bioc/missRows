\name{STATIS}

\alias{STATIS}

\title{
STATIS Analysis of Multiple Data Tables
}

\description{
Partially borrowed from the \code{statis} 
function implemented within the \pkg{ade4} package, \code{STATIS} 
performs a STATIS analysis of multiple data tables. This function is 
internally called by \code{\link{MIMFA}} and is not usually 
called directly by a user.
}

\usage{
STATIS(Ktab, nf = 3, tol = 1e-07)
}

\arguments{
    \item{Ktab}{a list containing the data tables. Tables in the 
        list should be arranged in samples \eqn{\times}{x} variables, 
        with samples order matching in all data tables.}

    \item{nf}{an integer indicating the number of kept axes for the 
        compromise configuration.}

    \item{tol}{a tolerance threshold to test whether an eigenvalue 
        is positive.}
}

\value{
A data frame with the row coordinates.
}

\references{
Lavit C., Escoufier Y., Sabatier R., Traissac P. (1994)
The ACT (STATIS method).
\emph{Computational Statistics & Data Analysis},
18(1), 97--119.
}

\seealso{
\code{\link[ade4]{statis}}
}

\examples{
#-- load data 
data(NCI60)

#-- STATIS analysis of two data tables
result <- STATIS(NCI60$completeData, nf = 3)
}

\keyword{multivariate}
