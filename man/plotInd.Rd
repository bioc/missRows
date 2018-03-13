\name{plotInd}

\encoding{latin1}
\alias{plotInd}

\title{
Plot of Individuals (Experimental Units)
}

\description{
This function provides scatter plots for individuals (experimental units)
representation from \code{\link{MIMFA}} results.
}

\usage{
plotInd(object, comp=1:2, colStrata=NULL, colMissing="white",
        confAreas=c("none", "ellipse", "convex.hull"), confLevel=0.95,
        ellipseType=c("norm", "t"), alpha=0.1, lwd=0.3, cex=3,
        legTitle="Strata")
}

\arguments{
    \item{object}{an object of class \code{MIDTList} as returned by
        \code{\link{MIMFA}} function.}

    \item{comp}{an integer vector of length two. The components
        that will be used on the horizontal and the vertical axis
        respectively to project the individuals.}

    \item{colStrata}{a character vector of the same length than the
        number of strata containing the color names to be used to
        annotate the individuals per stratum.}

    \item{colMissing}{the fill color for imputated individuals.}

    \item{confAreas}{a character string indicating whether to plot
        \code{"none"}, \code{"ellipse"} or \code{"convex.hull"}
        confidence areas.}

    \item{confLevel}{a numerical value indicating the confidence level
        of ellipses being plotted when \code{confAreas = "ellipse"}.
        The default is set to 0.95, for a 95\% confidence region.}

    \item{ellipseType}{the type of ellipse. The default \code{"norm"}
        assumes a multivariate normal distribution, and \code{"t"}
        assumes a multivariate t-distribution.}

    \item{alpha}{the \code{alpha} transparency for filled color of the
        confidence areas, values are any numbers from 0 (transparent)
        to 1 (opaque).}

    \item{lwd}{a positive number. The border line width of the
        confidence areas.}

    \item{cex}{a positive number. The amount by which plotting
        symbols should be magnified.}

    \item{legTitle}{character. The legend title.}
}

\details{
\code{plotInd} function makes scatter plot for individuals
representation from \code{\link{MIMFA}} results. Each point
corresponds to an individual. The individuals are colored with
rapport to their stratum, whereas imputed individuals are colored 
according to the \code{colMissing} argument.

Multiple imputation generates \code{M} imputed datasets and 
the variance between-imputations reflects the uncertainty associated 
to the estimation of the missing row values. The \code{plotInd} 
function propose two approaches to visualize the uncertainty 
due to missing data: confidence ellipses and convex hulls. The idea 
is to project all the multiple imputed datasets onto the compromise 
configuration. Each individual is represented by \code{M} points, each
corresponding to one of the \code{M} configurations. Confidence ellipses 
and convex hulls can then be constructed for the \code{M} configurations 
for each individual. For ease of understanding, not all individuals for 
the \code{M} configurations obtained are plotted.

Confidence ellipses can be created by setting the \code{confAreas} 
argument to \code{"ellipse"}. The 95\% confidence ellipses are showed by 
default. Convex hulls are plotted by setting the \code{confAreas} argument 
to \code{"convex.hull"}. The computed convex hull results in a polygon 
containing all \code{M} solutions.
}

\value{
An object of class \code{ggplot}.
}

\author{
Ignacio González and Valentin Voillet
}

\seealso{
\code{\link{plotVar}}
}

\examples{
#-- load data and create MIDTList object
data(NCI60)
midt <- newMIDTList(NCI60$incompleteData, strata=NCI60$cell.line)

#-- performs MIMFA
midt <- MIMFA(midt, ncomp=2, M=5)

#-- default plot
plotInd(midt)

#-- with confidence ellipses
plotInd(midt, confAreas="ellipse")

#-- with convex hull areas
plotInd(midt, confAreas="convex.hull")
}

\keyword{hplot}
\keyword{dplot}
