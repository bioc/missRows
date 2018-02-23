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
plotInd(object, comp = 1:2, col.strata = NULL, col.missing = "white",
        conf.areas = c("none", "ellipse", "convex.hull"), conf.level = 0.95,
        ellipse.type = c("norm", "t"), alpha = 0.1, lwd = 0.3, cex = 3,
        legend.title = "Strata")
}

\arguments{
  \item{object}{object of class \code{MIDTList} as returned by
        \code{\link{MIMFA}} function.}

  \item{comp}{integer vector of length two. The components
        that will be used on the horizontal and the vertical axis
        respectively to project the individuals.}

  \item{col.strata}{character vector of the same length than the
        number of strata containing the color names to be used to
        annotate the individuals per stratum.}

  \item{col.missing}{the fill color for imputated individuals.}

  \item{conf.areas}{character string indicating whether to plot
        \code{'none'}, \code{'ellipse'} or \code{'convex.hull'}
        confidence areas.}

  \item{conf.level}{numerical value indicating the confidence level
        of ellipses being plotted when \code{conf.areas = 'ellipse'}.
        The default is set to 0.95, for a 95\% region.}

  \item{ellipse.type}{the type of ellipse. The default \code{"norm"}
        assumes a multivariate normal distribution, and \code{"t"}
        assumes a multivariate t-distribution.}

  \item{alpha}{the \code{alpha} transparency for filled color of the
        confidence areas, values are any numbers from 0 (transparent)
        to 1 (opaque).}

  \item{lwd}{a positive number. The border line width of the
        confidence areas.}

  \item{cex}{a positive number. The amount by which plotting
        symbols should be magnified.}

  \item{legend.title}{character. The legend title.}
}

\details{
\code{plotInd} function makes scatter plot for individuals
representation from \code{\link{MIMFA}} results. Each point
corresponds to an individual. The individuals are colored with
rapport to their stratum. Imputated individuals are colored with
\code{col.missing} color.

Confidences areas ...
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
midt <- newMIDTList(NCI60$incompleteData, strata = NCI60$cell.line)

#-- performs MIMFA
midt <- MIMFA(midt, ncomp = 2, M = 5)

#-- default plot
plotInd(midt)

#-- with confidence ellipses
plotInd(midt, conf.areas = "ellipse")

#-- with convex hull areas
plotInd(midt, conf.areas = "convex.hull")
}

\keyword{hplot}
\keyword{dplot}
