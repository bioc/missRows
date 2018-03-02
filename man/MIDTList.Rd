\name{MIDTList-class}

\encoding{latin1}
\docType{class}
\alias{MIDTList}
\alias{MIDTList-class}
\alias{incompleteData}
\alias{incompleteData,MIDTList-method}
\alias{tableNames}
\alias{tableNames,MIDTList-method}
\alias{strata}
\alias{strata,MIDTList-method}
\alias{missingRows}
\alias{missingRows,MIDTList-method}
\alias{compromise}
\alias{compromise,MIDTList-method}
\alias{configurations}
\alias{configurations,MIDTList-method}
\alias{imputedRows}
\alias{imputedRows,MIDTList-method}
\alias{MIparam}
\alias{MIparam,MIDTList-method}

\title{Infrastructure for Multiple Imputation and Data Tables List}

\description{
\code{MIDTList} is an S4 class providing the infrastructure to store the
input data tables, intermediate calculations and results of a multiple
imputation approach.
}

\section{Objects from the Class}{
Objects of this class are typically created by \code{\link{newMIDTList}}
constructor and returned by \code{\link{MIMFA}} function. Usually there
should be no need to create them directly.
}

\section{Slots}{
Initial slots:

  \describe{
    \item{\code{incompleteData}:}{a list containing the data tables 
          with missing rows. Tables in the list should be arranged 
          in samples \eqn{\times}{x} variables, with samples order 
          matching in all data tables.}

    \item{\code{strata}:}{named vector with the individual names
          and the corresponding condition (stratum).}

    \item{\code{tableNames}:}{a character vector of the same length as
          \code{incompleteData} with the name of the each table.}

    \item{\code{missingRows}:}{a list containing character vectors with 
          the name of the missing individuals (rows) per table.}\cr

  Slots returned after \code{\link{MIMFA}}:

    \item{\code{compromise}:}{the compromise configuration, a 
          matrix with the individuals coordinates as returned by
          \code{\link{STATIS}} function.}

    \item{\code{configurations}:}{a list containing the individuals
          coordinates for each imputed dataset as returned by
          \code{\link{MIMFA}} function.}

    \item{\code{imputedRows}:}{a list containing the imputed rows
          (individuals) for each data table as returned by
          \code{\link{imputeDataMFA}} function.}

    \item{\code{MIparam}:}{a list containing the parameters used
          in the \code{\link{MIMFA}} function.}
  }
}

\section{Extends}{
Class \code{\linkS4class{MIDTList}}, directly.
}

\section{Methods}{
Class-specific methods return the corresponding objects:

  \describe{

    \item{incompleteData}{\code{signature(object = "MIDTList")}:
          Returns list containing the initial data tables.}

    \item{tableNames}{\code{signature(object = "MIDTList")}:
          Access table names in \code{incompleteData}.}

    \item{strata}{\code{signature(object = "MIDTList")}:
          Returns factor of strata giving the stratum for
          each individual.}

    \item{missingRows}{\code{signature(object = "MIDTList")}:
          Returns list containing character vectors with 
          the name of the missing individuals (rows) per table.}

    \item{compromise}{\code{signature(object = "MIDTList")}:
          Returns matrix with the individuals coordinates as returned 
          by \code{\link{STATIS}} function.}

    \item{configurations}{\code{signature(object = "MIDTList", M = "all")}:
          Gets all configurations. If \code{M} is a positive integer, the
          \code{M}th configuration is returned.}

    \item{imputedRows}{\code{signature(object = "MIDTList")}:
          Returns list containing the imputed rows
          (individuals) for each data table.}

    \item{MIparam}{\code{signature(object = "MIDTList")}:
          Returns list containing the parameters used
          in the \code{\link{MIMFA}} function.}
  }

Standard generic methods:

  \describe{

    \item{show}{\code{signature(object = "MIDTList")}: Informatively
          display object contents.}
	}
}

\author{
Ignacio González
}

\seealso{
\code{\link{newMIDTList}}
}

\keyword{classes}
