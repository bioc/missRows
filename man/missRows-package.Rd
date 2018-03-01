\name{missRows-package}

\encoding{latin1}
\alias{missRows-package}
\alias{missRows}

\docType{package}

\title{
Handling Missing Rows in Multi-Omics Data Integration
}

\description{
More about what it does (maybe more than one line)
 A concise (1-5 lines) description of the package
}

\details{
\tabular{ll}{
Package: \tab missRows\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2017-04-22\cr
License: \tab What license is it under?\cr
Depends: \tab methods\cr
}
 An overview of how to use the package, including the most important functions
}

\author{
Ignacio González and Valentin Voillet

Maintainer: Who to complain to <yourfault@somewhere.net>
}

\references{
Voillet V., Besse P., Liaubet L., San Cristobal M., González I. (2016).
Handling missing rows in multi-omics data integration:
Multiple Imputation in Multiple Factor Analysis framework.
\emph{BMC Bioinformatics}, 17(40).
}

\keyword{package}
\seealso{
 Optional links to other man pages, e.g.
% \code{\link[<pkg>:<pkg>-package]{<pkg>}}
}

\examples{
## A typical MI-MFA session might look like the following. 
## Here we assume there are two data tables with missing rows, 
## "table1" and "table2", and the stratum for each individual 
## is stored in a data frame df.

\dontrun{

#-- Data preparation
midt <- newMIDTList(table1, table2, strata = df)

#-- Performing MI
midt <- MIMFA(midt, ncomp = 2, M = 30)

#-- Analysis of the results
plotInd(midt)
plotVar(midt)
}
}
