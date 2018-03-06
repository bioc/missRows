\name{newMIDTList}

\encoding{latin1}
\docType{methods}
\alias{newMIDTList}
\alias{tableNames<-}
\alias{tableNames<-,MIDTList-method}

\title{
Creates a MIDTList Object
}

\description{
Creates a \code{MIDTList} object from data tables, stratum indicator
for each individual and table names (optional).
}

\usage{
newMIDTList(..., strata = NULL, tableNames = NULL)

\S4method{tableNames}{MIDTList}(object) <- value
}

\arguments{
  \item{\dots}{data tables with missing rows. Two or more objects
        which can be interpreted as matrices (or data frames).
        Data tables passed as arguments in \code{\dots} must be
        arranged in samples \eqn{\times} variables, with samples
        order and row names matching in all data tables.}

  \item{strata}{a named vector or data frame giving the stratum for
        each individual. Names (in vector) or row names (in data frame)
        must be equal to the row names among data tables.}

  \item{tableNames}{optional. A character vector giving the name for
        each table.}

  \item{object}{an object of class \code{MIDTList}.}

  \item{value}{a character vector of length the total number of
        data tables.}
}

\details{
To facilitate programming pipelines, \code{NULL} values are input for
\code{compromise}, \code{configurations}, \code{imputedRows} and
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

#-- MIDTList object directly from a list of data tables
midt <- newMIDTList(NCI60$incompleteData, strata = NCI60$cell.line)
midt

#-- MIDTList object from separate data tables
table1 <- NCI60$incompleteData$trans
table2 <- NCI60$incompleteData$prote

midt <- newMIDTList(table1, table2, strata = NCI60$cell.line,
                    tableNames = c("transcrip", "proteome"))
midt
}

\keyword{classes}
