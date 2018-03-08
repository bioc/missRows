\name{plotVar}

\encoding{latin1}
\alias{plotVar}

\title{
Plot of Variables: Correlation Circle
}

\description{
This function provides "correlation circle", scatter plots for variables
representation from \code{\link{MIMFA}} results.
}

\usage{
plotVar(object, comp = 1:2, col = NULL, var.names = FALSE, cex = 3,
        pch = 19, alpha = 0.7, spty = TRUE, cutoff = 0, rad.in = 0.5,
        overlap = TRUE, ncols = 2, legend.title = "Tables")
}

\arguments{

    \item{object}{an object of class \code{MIDTList} as returned by
        \code{\link{MIMFA}} function.}

    \item{comp}{an integer vector of length two. The components that
        will be used on the horizontal and the vertical axis
        respectively to project the variables.}

    \item{col}{a character or integer vector of colors for plotted
        character and symbols, it must be of length the total number 
        of data tables (see Details).}

    \item{var.names}{either a character vector with data table 
        names, or \code{FALSE} for no plotting variable
        labels. If \code{TRUE}, the variable names in the data
        tables are used as variable labels in the plot. See Details.}

    \item{cex}{a numeric vector of character expansion sizes for the
        plotted character and symbols, can be of length one
        or of length the total number of data tables. See Details.}

    \item{pch}{plotting 'character'. A vector of single characters or
        integers, can be of length one or of length the total number
        of data tables (see Details). See \code{\link{points}} for
        all alternatives.}

    \item{alpha}{the \code{alpha} transparency for plotting color of
        the symbols, values are any numbers from 0 (transparent)
        to 1 (opaque).}

    \item{spty}{logical, specifying the type of plot region to be used.
        If \code{TRUE} (the default), a square plotting region is
        generated. If not, a maximal plotting region is produced.}

    \item{cutoff}{a numeric between 0 and 1. Variables with correlations
        below this cutoff in absolute value are not plotted (see Details).}

    \item{rad.in}{a numeric between 0 and 1, the radius of the inner circle.
        Defaults to 0.5.}

    \item{overlap}{logical. Whether the variables in data tables should be
        plotted in one single panel. Default is \code{TRUE}.}

    \item{ncols}{numeric. When \code{overlap = FALSE} subsequent figures will
        be drawn in a multipanel on the device with \code{ncols} columns.}

    \item{legend.title}{character. The legend title.}
}

\details{
\code{plotVar} produce a "correlation circle", \emph{i.e.} the 
correlations between each variable and the selected components 
are plotted as scatter plot, with concentric circles of radius 
one and radius given by \code{rad.in}. Each point corresponds 
to a variable.

The \code{var.names} argument can be used in order to select a 
part of the data table variable labels that are drawn. For example 
if you have data tables named \code{"table1"} and \code{"table2"}, 
you can use \code{var.names = "table1"} and then the variable 
names in the \code{"table1"} are drawn.

The arguments \code{cex} and \code{pch} can be either vectors 
of length one or of length the total number of data tables. In the 
first case, the single value determine the graphics attributes for 
all data tables variables. Otherwise, multiple argument values can 
be specified so that each data table variable can be given its own 
graphic attributes. In this case, each component of the vector 
correspond to the attributes of the each data table variables.
}

\value{
A list containing the following components:

    \item{df}{a data frame used to generate the \code{ggplot}.}

    \item{ggp}{an object of class \code{ggplot}.}
}

\author{
Ignacio González
}

\seealso{
\code{\link{plotInd}}
}

\examples{
#-- load data and create MIDTList object
data(NCI60)
midt <- newMIDTList(NCI60$incompleteData, strata = NCI60$cell.line,
                    tableNames = c("trans", "prote"))

#-- performs MIMFA
midt <- MIMFA(midt, ncomp = 2, M = 5)

#-- default plot
plotVar(midt)

#-- select data table variables to draw and cutoff
plotVar(midt, var.names = "trans", cutoff = 0.55)
plotVar(midt, var.names = TRUE, cutoff = 0.55)
}

\keyword{hplot}
\keyword{dplot}
