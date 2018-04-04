\name{imputedData}

\encoding{latin1}
\alias{imputedData}

\title{
Get the Imputed Data Tables
}

\description{
\code{imputedData} function exports the list of data tables 
with imputed data.
}

\usage{
imputedData(object)
}

\arguments{
    \item{object}{an object of class \code{MIDTList} as returned by
        \code{\link{MIMFA}} function.}
}

\value{
A list of length the data tables number, each component
containing a completed data table.
}

\author{
Ignacio González
}

\seealso{
\code{\link{MIMFA}}
}

\examples{
#-- load data and create MIDTList object
data(NCI60)
midt <- newMIDTList(NCI60$mae)

#-- performs MIMFA
midt <- MIMFA(midt, ncomp=2, M=5)

#-- exports the imputed data tables
completeData <- imputedData(midt)
}

\keyword{utilities}
