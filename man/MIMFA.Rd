\name{MIMFA}

\encoding{latin1}
\alias{MIMFA}

\title{
Handling Missing Individuals in MFA
}

\description{
The \code{MIMFA} function estimates coordinates of
individuals and variables on the MFA components by implementing 
a multiple imputation (MI) approach in order to deal
with multiple tables in presence of missing individuals
(rows).
}

\usage{
MIMFA(object, ncomp = 2, M = NULL, estim.ncp = FALSE,
      max.iter = 500, tol = 1e-10)
}

\arguments{
  \item{object}{object of class \code{MIDTList} as returned by
        \code{\link{newMIDTList}}.}

  \item{ncomp}{number of components to include in MFA when 
        \code{estim.ncp = FALSE} (default to 2). If 
        \code{estim.ncp = TRUE}, then \code{ncomp} correspond 
        to the maximum number of components to test.}

  \item{M}{integer, number of imputations. Default to 
        \code{min(30, M.total)}, where \code{M.total} is the 
        total number of possible imputations.}

  \item{estim.ncp}{logical. If \code{TRUE} the number of MFA 
        components for data imputation is estimated. Default 
        is \code{FALSE}.}

  \item{max.iter}{integer, maximum number of iterations for the 
        \code{\link{imputeDataMFA}} function.}

  \item{tol}{positive value, the threshold for assessing convergence
        in the \code{\link{imputeDataMFA}} algorithm.}
}

\details{
According to the MI methodology, missing rows are filled
in by several sets of plausible values, resulting in \code{M}
completed data. MFA is then applied to each completed data
leading to \code{M} different configurations.
Finally, the \code{M} configurations are combined using the
STATIS method to yield one consensus solution.

If \code{estim.ncp = TRUE}, the number of MFA components for 
data imputation is estimated using the generalized cross-validation 
approximation method. In this case, \code{ncomp} correspond 
to the maximum number of components to test.
}

\value{
A \code{\link{MIDTList}} object containing additional slots for:

  \describe{
    \item{\code{compromise}}{}
    \item{\code{configurations}}{}
    \item{\code{imputedRows}}{}
    \item{\code{MIparam}}{}
  }
  
See \code{\link{MIDTList}} for description.
}

\references{
Voillet V., Besse P., Liaubet L., San Cristobal M., González I. (2016).
Handling missing rows in multi-omics data integration:
Multiple Imputation in Multiple Factor Analysis framework.
\emph{BMC Bioinformatics}, 17(40).

Lavit C., Escoufier Y., Sabatier R., Traissac P. (1994).
The ACT (STATIS method).
\emph{Computational Statistics & Data Analysis},
18(1), 97--119.

Josse J., Husson F. (2012).
Selecting the number of components in PCA using cross-validation
approximations.
\emph{Computational Statistics and Data Analysis}, 56, 1869--1879.
}

\author{
Ignacio González and Valentin Voillet
}

\seealso{
\code{\link{plotInd}}, \code{\link{plotVar}}, \code{\link{tuneM}}
}

\examples{
#-- load data and create MIDTList object
data(NCI60)
midt <- newMIDTList(NCI60$incompleteData, strata = NCI60$cell.line)
midt

#-- performs MIMFA
midt <- MIMFA(midt, ncomp = 3, M = 10)
midt

#-- estimates the number of MFA components for data imputation
#-- ncomp is chosen to being enough large
\dontrun{
midt <- MIMFA(midt, ncomp = 50, M = 10, estim.ncp = TRUE)
midt}
}

\keyword{multivariate}
