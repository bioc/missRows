\name{tuneM}

\encoding{latin1}
\alias{tuneM}
\alias{print.tuneM}

\title{
Tune the Number of Imputations in MI-MFA
}

\description{
\code{tuneM} can be used to determine the appropriate number 
of imputed datasets needed to obtain satisfactory results 
in MI-MFA.
}

\usage{
tuneM(object, ncomp = 2, M.max = 30, inc = 5, N = 10, 
        tol = 1e-06, show.plot = TRUE)

\method{print}{tuneM}(x, \ldots)
}

\arguments{
    \item{object}{an object of class \code{MIDTList} as returned by
        \code{\link{newMIDTList}}.}

    \item{x}{an object of class inheriting from \code{tuneM}.}

    \item{ncomp}{a number of components to include in MFA.}

    \item{M.max}{an integer corresponding to the maximum number of 
        imputed datasets. See Details.}

    \item{inc}{integer. The increment of the sequence for the
        number \code{M} of imputations considered. See Details.}

    \item{N}{integer. Collections of size \code{N} are generated for 
        each number of imputations \code{M}. See Details.}

    \item{tol}{a positive value, the tolerance used for assessing 
        stabilization.}

    \item{show.plot}{logical. If \code{TRUE} (the default) a plot 
        showing the stability of the estimated MFA configurations 
        is displayed.}

    \item{\ldots}{not currently used.}
}

\details{
The appropriate number of imputations can be informally determined 
by carrying out MI-MFA on \eqn{N} replicate sets of \eqn{M_l} 
imputations for \eqn{l=0,1,2,\ldots ,}{l = 0, 1, 2, . . . ,} with 
\eqn{M_0 < M_1 < M_2 < \cdots < M_{max}}{M_0 < M_1 < M_2 < . . . < M_max},
until the estimate compromise configurations are stabilized.

\code{tuneM} function implements such a procedure. Collections 
of size \code{N} are generated for each number of imputations 
\code{M}, with \code{M = seq(inc, M.max, by = inc)}. The stability 
of the estimated MI-MFA configurations is then determined by 
calculating the RV coefficient between the configurations obtained 
using \code{M}\eqn{_{\,l}}{_l} and \code{M}\eqn{_{\,l+1}}{_{l+1}} 
imputations.

If \code{show.plot = TRUE} a plot showing the stability of the 
estimated MFA configurations is displayed. The values shown are 
the mean RV coefficients for the \code{N} configurations as a 
function of the number of imputations. Error bars represent the 
standard deviation of the RV coefficients.
}

\value{
A list with the following components:

    \item{stats}{\code{data.frame} containing the information
        used to generate the plot.}

    \item{ggp}{an object of class \code{ggplot}.}
}

\references{
Voillet V., Besse P., Liaubet L., Cristobal M.S., González I.
(2016). Handling missing rows in multi-omics data integration:
Multiple Imputation in Multiple Factor Analysis framework.
\emph{BMC Bioinformatics}, 17(40).
}

\author{
Ignacio González, Valentin Voillet
}

\seealso{
\code{\link{MIMFA}}
}

\examples{
#-- load data and create MIDTList object
data(NCI60)
midt <- newMIDTList(NCI60$incompleteData, strata = NCI60$cell.line)

#-- tune the number of imputations
\dontrun{
tune <- tuneM(midt, ncomp = 2, M.max = 100, inc = 10, N = 10)
tune}
}

\keyword{multivariate}
\keyword{dplot}
