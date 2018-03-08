\name{missPattern}

\encoding{latin1}
\alias{missPattern}
\alias{print.missPattern}

\title{
Inspects Pattern of Missingness
}

\description{
This function inspects and plots the structure of missing rows
in data tables.
}

\usage{
missPattern(object, col.strata = NULL, col.missing = "grey70", cex.ttitle = 12,
            legend.title = "Strata", miss.lab = "miss", show.plot = TRUE)

\method{print}{missPattern}(x, \ldots)
}

\arguments{
    \item{object}{an object of class \code{MIDTList} as returned by
        \code{\link{newMIDTList}}.}

    \item{x}{an object of class inheriting from \code{missPattern}}

    \item{col.strata}{a character vector of the same length than the
        number of strata, containing the color names to be used to
        annotate the individuals per stratum.}

    \item{col.missing}{the fill color for missing individuals (rows).}

    \item{cex.ttitle}{a positive number. The amount by which table
        titles should be magnified.}

    \item{legend.title}{character. The legend title.}

    \item{miss.lab}{character. The label legend for missing rows.}

    \item{show.plot}{logical. Whether the plot will be displayed. 
        Default is \code{TRUE}.}

    \item{\ldots}{not used currently.}
}

\details{
\code{missPattern} calculates and plots the amount of
missing/available rows in each stratum per data table.
Data tables are plotted separately on a same device showing
the pattern of missingness. The individuals are colored according
to their stratum whereas missing rows are specific
colored (see \code{col.missing}).
}

\value{
A list with the following components:

    \item{nbMissing}{\code{data.frame} containing the amount of
        missing/available rows in each stratum per data table.}

    \item{isMissing}{\code{data.frame} containing the indicator
        matrix for the missing rows.}

    \item{ggp}{an object of class \code{ggplot}.}
}

\author{
Ignacio González
}

\examples{
#-- load data and create MIDTList object
data(NCI60)
midt <- newMIDTList(NCI60$incompleteData, strata = NCI60$cell.line)

#-- inspects pattern of missingness
patt <- missPattern(midt)
patt
}

\keyword{hplot}
\keyword{dplot}
