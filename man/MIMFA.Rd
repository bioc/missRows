\name{MIMFA}

\encoding{latin1}
\alias{MIMFA}

\title{
Handling Missing Individuals in MFA
}

\description{
The \code{MIMFA} function estimates coordinates of
individuals on the MFA components by implementing a
multiple imputation (MI) approach in order to deal
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

  \item{ncomp}{number of components to include in MFA. Default to 2.}

  \item{M}{integer, number of imputations (default \code{min(30, M.total)},
        where \code{M.total} is the total number of possible imputations.}

  \item{estim.ncp}{Describe \code{estim.ncp} here}

  \item{max.iter}{Describe \code{max.iter} here}

  \item{tol}{Describe \code{tol} here}
}

\details{
According to the MI methodology, missing rows are filled
in by several sets of plausible values, resulting in \code{M}
completed data. MFA is then applied to each completed data
leading to \code{M} different configurations.
Finally, the \code{M} configurations are combined using the
STATIS method to yield one consensus solution.

\code{estim.ncp} ... Estimate the number of components in ...
method used to estimate the number of components, "GCV" for the
generalized cross-validation approximation.
Returns ncp the best number of dimensions to use (find the minimum
or the first local minimum) and the mean error for each dimension
tested
}

\value{
A \code{\link{MIDTList}} object containing additional slots for:

  \describe{
    \item{compromise}{}
    \item{configurations}{}
    \item{imputedRows}{}
    \item{MIparam}{}
  }
See \code{\link{MIDTList}} for description.
}

\references{
Voillet V., Besse P., Liaubet L., San Cristobal M., González I. (2016).
Handling missing rows in multi-omics data integration:
Multiple Imputation in Multiple Factor Analysis framework.
\emph{BMC Bioinformatics}, 17(40).

Lavit C., Escoufier Y., Sabatier R., Traissac P. (1994)
The ACT (STATIS method).
\emph{Computational Statistics & Data Analysis},
18(1), 97-119.

Josse J., Husson F. (2012).
Selecting the number of components in PCA using cross-validation
approximations.
\emph{Computational Statistics and Data Analysis}, 56, 1869-1879.
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
}

\keyword{ kwd1 }
\keyword{ kwd2 }
