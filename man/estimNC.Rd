\name{estimNC}

\alias{estimNC}

\title{
Estimate the Number of Components for Data Imputation
}

\description{
\code{estimNC} estimates the number of MFA components for 
data imputation. This function is called internally 
by \code{\link{MIMFA}} and is not usually called
directly by a user.
}

\usage{
estimNC(X, minNC=0, maxNC)
}

\arguments{
    \item{X}{a numeric matrix.}

    \item{minNC}{minimum number of components to consider, by default 0.}

    \item{maxNC}{maximum number of components to test.}
}

\details{
Partially borrowed from the \code{estim_npc} 
function in the \pkg{FactoMineR} package, \code{estimNC} estimates 
the number of MFA components for data imputation using the 
generalized cross-validation approximation method. 
}

\value{
Return the number of MFA components to use in data imputation. 
}

\references{
Josse, J. and Husson, F. (2012).
Selecting the number of components in PCA using cross-validation 
approximations. 
\emph{Computational Statistics and Data Analysis}, \bold{56}, 1869-1879.
}

\examples{
#-- load data and create MIDTList object
data(NCI60)
midt <- newMIDTList(NCI60$incompleteData, strata=NCI60$cell.line)

#-- estimates the number of MFA components for data imputation
#-- ncomp is chosen to being enough large
\dontrun{
midt <- MIMFA(midt, ncomp=50, M=10, estimNC=TRUE)
midt}
}

\keyword{multivariate}
