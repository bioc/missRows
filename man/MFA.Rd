\name{MFA}

\encoding{latin1}
\alias{MFA}

\title{
Multiple Factor Analysis (MFA)
}

\description{
Perform Multiple Factor Analysis on quantitative variables.
This function is internally called by \code{\link{MIMFA}} and
is not usually called directly by a user.
}

\usage{
MFA(dataTables, ncomp, nbRows, nbTables, ncTables)
}

\arguments{
    \item{dataTables}{a list containing data tables without missing
        data. Tables in the list should be arranged in samples
        \eqn{\times}{x} variables, with samples order matching in all
        data tables.}

    \item{ncomp}{a number of components to include in MFA.}

    \item{nbRows}{a number of rows, equal for all data tables.}

    \item{nbTables}{a number of data tables in \code{dataTables}.}

    \item{ncTables}{a vector containing the number of columns per data
        table.}
}

\details{
\code{MFA} function performs Multiple Factor Analysis in the sense of
Escofier-Pages on data tables of quantitative variables.
}

\value{
\code{MFA} returns a matrix of individuals coordinates.
}

\references{
Escofier, B. and Pages, J. (1994) Multiple Factor Analysis (AFMULT package).
\emph{Computational Statistics and Data Analysis}, \bold{18}, 121-140.
}

\author{
Ignacio González, Valentin Voillet
}

\seealso{
\code{\link[FactoMineR]{MFA}} in \code{FactoMineR}
}

\examples{
#-- load complete data
data(NCI60)
dataTables <- NCI60$completeData
nc <- sapply(dataTables, ncol)

mfa <- MFA(dataTables, ncomp=3, nbRows=nrow(dataTables[[1]]), 
            nbTables=length(dataTables), ncTable=nc)
}

\keyword{multivariate}
