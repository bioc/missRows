\name{MIDTList-class}

\encoding{latin1}
\docType{class}
\alias{MIDTList}
\alias{MIDTList-class}
\alias{initialize,MIDTList-method}
\alias{strata}
\alias{strata,MIDTList-method}
\alias{missingIndv}
\alias{missingIndv,MIDTList-method}
\alias{compromise}
\alias{compromise,MIDTList-method}
\alias{configurations}
\alias{configurations,MIDTList-method}
\alias{imputedIndv}
\alias{imputedIndv,MIDTList-method}
\alias{MIparam}
\alias{MIparam,MIDTList-method}

\title{Infrastructure for Multiple Imputation and Multi-omics Experiments}

\description{
\code{MIDTList} is an S4 class that extends the class 
\code{\link{MultiAssayExperiment}} by providing the infrastructure (slots) 
to store the input data, intermediate calculations and results of 
a multiple imputation approach.
}

\section{Objects from the Class}{
Objects can be created by calls of the form:

\code{MIDTList(..., colData=NULL, strata=NULL, assayNames=NULL)}

\code{new("MIDTList", ..., colData=NULL, strata=NULL, assayNames=NULL)}
}

\section{Arguments}{
\describe{
    \item{\code{\dots}}{arguments passed to '\code{\dots}' can be:
    
        1) data tables with missing individuals. Two or more objects
        which can be interpreted as matrices (or data frames). Data tables
        passed as arguments in \code{\dots} must be arranged in variables 
        (rows) \eqn{\times}{x} individuals (columns), with individual names
        matching row names of \code{colData}.
        
        2) a \code{list} containing two or more data tables with
        missing individuals. Data tables (matrices or data frames) within 
        list must be arranged in variables (rows) \eqn{\times}{x} individuals
        (columns), with individual names matching row names of \code{colData}.
        
        3) an object of class 
        \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}. In this case
        \code{colData} and \code{assayNames} arguments are ignored.
        }

    \item{\code{colData}}{a \code{DataFrame} giving the characteristics for 
        all individuals (biological units). The row names of \code{colData}
        must contain individual identifiers.}

    \item{\code{assayNames}}{optional. A character vector giving the name for
        each table.}

    \item{\code{strata}}{a character indicating the column of \code{colData} 
        to be used as \code{strata} in the construction of \code{MIDTList}.}
    }
}

\section{Details}{
To facilitate programming pipelines, \code{NULL} values are input for
\code{compromise}, \code{configurations}, \code{imputedIndv} and
\code{MIparam} slots, in which case the default value is used as if the
argument had been missing. These slots will be updated after multiple
imputation (\code{\link{MIMFA}}) approach.
}

\section{Slots}{

\describe{
    
    \item{\code{ExperimentList}}{an \code{\link{ExperimentList}} class 
        object for each assay dataset.}

    \item{\code{colData}}{a \code{DataFrame} of all clinical/specimen data
        available across experiments.}

    \item{\code{sampleMap}}{a \code{DataFrame} of translatable identifiers
        of samples (individuals) and participants.}

    \item{\code{metadata}}{additional data describing the 
        \code{MultiAssayExperiment} object.}

    \item{\code{drops}}{a metadata \code{list} of dropped information.}

    \item{\code{strata}:}{a numeric value or character. The column 
        of \code{colData} to be used as strata in \code{MIMFA}.}

    \item{\code{missingIndv}:}{a list containing character vectors with 
        the name of the missing individuals per table.}

    \item{\code{compromise}:}{the compromise configuration, a 
        matrix with the individuals coordinates as returned by
        \code{\link{STATIS}} function.}

    \item{\code{configurations}:}{a list containing the individuals
        coordinates for each imputed dataset as returned by
        \code{\link{MIMFA}} function.}

    \item{\code{imputedIndv}:}{a list containing the imputed 
        individuals for each data table as returned by
        \code{\link{imputeDataMFA}} function.}

    \item{\code{MIparam}:}{a list containing the parameters used
        in the \code{\link{MIMFA}} function.}
    }
}

\section{Extends}{
Class \code{\linkS4class{MultiAssayExperiment}}, directly.
}

\section{Methods}{

\describe{
    \item{initialize}{\code{signature(.Object = "MIDTList")}: See 'Objects 
        from the Class' section for description.}
}

Class-specific methods return the corresponding objects:

\describe{

    \item{strata}{\code{signature(object="MIDTList")}:
        Return factor of strata giving the stratum for
        each individual.}

    \item{missingIndv}{\code{signature(object="MIDTList")}:
        Return list containing character vectors with 
        the name of the missing individuals per table.}

    \item{compromise}{\code{signature(object="MIDTList")}:
        Return matrix with the individuals coordinates as returned 
        by \code{\link{STATIS}} function.}

    \item{configurations}{\code{signature(object="MIDTList", M="all")}:
        Get all configurations. If \code{M} is a positive integer, the
        \code{M}th configuration is returned.}

    \item{imputedIndv}{\code{signature(object="MIDTList")}:
        Return list containing the imputed individuals for each data table.}

    \item{MIparam}{\code{signature(object="MIDTList")}:
        Return list containing the parameters used
        in the \code{\link{MIMFA}} function.}
}

Standard generic methods:

\describe{

    \item{show}{\code{signature(object="MIDTList")}: Informatively
        display object contents.}
}
    
See \link[MultiAssayExperiment]{MultiAssayExperiment-class} for generic
methods associated to the \code{MultiAssayExperiment} class.
}

\author{
Ignacio González
}

\seealso{
\link{MultiAssayExperiment-class},
\link{MultiAssayExperiment-methods}
}

\examples{
#-- load data
data(NCI60)

#-- MIDTList object from separate data tables
table1 <- NCI60$dataTables$trans
table2 <- NCI60$dataTables$prote
colData <- NCI60$dataTables$cell.line

midt <- MIDTList(table1, table2, colData=colData,
                assayNames=c("transcrip", "proteome"))
midt

#-- MIDTList object from a list
tablesList <- NCI60$dataTables[1:2]
colData <- NCI60$dataTables$cell.line

midt <- MIDTList(tablesList, colData=colData)
midt

#-- MIDTList object directly from a 'MultiAssayExperiment'
midt <- MIDTList(NCI60$mae)
midt
}

\keyword{classes}
