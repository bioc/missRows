\name{MIDTList-class}

\encoding{latin1}
\docType{class}
\alias{MIDTList}
\alias{MIDTList-class}
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
to store the data, intermediate calculations and results of 
a multiple imputation approach.
}

\section{Objects from the Class}{
Objects of this class are typically created by \code{\link{newMIDTList}}
constructor and returned by \code{\link{MIMFA}} function. Usually there
should be no need to create them directly.
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
\code{\link{newMIDTList}}, \link{MultiAssayExperiment-class},
\link{MultiAssayExperiment-methods}
}

\examples{
#-- load data
data(NCI60)

#-- MIDTList object from a list of data tables
midt <- newMIDTList("trans" = NCI60$dataTables$trans, 
                    "prote" = NCI60$dataTables$prote,
                    colData = NCI60$dataTables$cell.line)
midt
}

\keyword{classes}
