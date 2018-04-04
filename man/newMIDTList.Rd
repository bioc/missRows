\name{newMIDTList}

\encoding{latin1}
\docType{methods}
\alias{newMIDTList}
\alias{newMIDTList.default}
\alias{newMIDTList.MultiAssayExperiment}

\title{
Creates an MIDTList Object
}

\description{
Creates a \code{MIDTList} object from data tables, stratum indicator
for each individual and table names (optional).
}

\usage{
\method{newMIDTList}{default}(..., colData=NULL, strata=NULL, assayNames=NULL)

\method{newMIDTList}{MultiAssayExperiment}(object, strata=NULL, ...)
}

\arguments{
    \item{\dots}{data tables with missing individuals. Two or more objects
        which can be interpreted as matrices (or data frames) or a list
        containing these. Data tables passed as arguments in \code{\dots} 
        (or \code{list(\dots)}) must be arranged in variables (rows) 
        \eqn{\times}{x} individuals (columns), with individual names matching 
        row names of \code{colData}. 
        
        For \code{'MultiAssayExperiment'}, \code{\dots} is not currently used.}
        
    \item{object}{an object of class 
        \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}.}

    \item{colData}{a data frame giving the stratum for each individual. 
        The row names of \code{colData} must contain individual identifiers.}

    \item{assayNames}{optional. A character vector giving the name for
        each table.}

    \item{strata}{a character ...}
}

\details{
To facilitate programming pipelines, \code{NULL} values are input for
\code{compromise}, \code{configurations}, \code{imputedIndv} and
\code{MIparam} slots, in which case the default value is used as if the
argument had been missing. These slots will be updated after multiple
imputation (\code{\link{MIMFA}}) approach.
}

\value{
A \code{\link{MIDTList}} object.
}

\author{
Ignacio González
}

\seealso{
\code{\linkS4class{MIDTList}}, \code{\link{MIMFA}}
}

\examples{
#-- load data
data(NCI60)

#-- MIDTList object directly from a 'MultiAssayExperiment' object
midt <- newMIDTList(NCI60$mae)
midt

#-- MIDTList object from separate data tables
table1 <- NCI60$dataTables$trans
table2 <- NCI60$dataTables$prote
colData <- NCI60$dataTables$cell.line

midt <- newMIDTList(table1, table2, colData=colData,
                    assayNames=c("transcrip", "proteome"))
midt
}

\keyword{classes}
