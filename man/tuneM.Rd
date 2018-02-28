\name{tuneM}

\encoding{latin1}
\alias{tuneM}
\alias{print.tuneM}

\title{
tuneM
}

\description{
 A concise (1-5 lines) description of what the function does.
}

\usage{
tuneM(object, ncomp = 2, M.max = 30, inc = 5, N = 10, tol = 1e-06,
      show.plot = TRUE)

\method{print}{tuneM}(x, \ldots)
}

\arguments{
  \item{object}{an object of class \code{MIDTList} as returned by
        \code{\link{newMIDTList}}.}

  \item{x}{object of class inheriting from \code{tuneM}}

  \item{ncomp}{number of components to include in MFA.}

  \item{M.max}{integer corresponding to the maximum number of 
        imputed datasets. See Details.}

  \item{inc}{integer. \Rcode{M} number of imputations are considered 
        for \code{M = seq(inc, M.max, by = inc)}. See Details.}

  \item{N}{integer. Collections of size \code{N} are generated for 
        each number of imputations \code{M}. See Details.}

  \item{tol}{positive value, the tolerance used for assessing 
        stabilization.}

  \item{show.plot}{logical. Whether the plot will be displayed. 
        Default is \code{TRUE}.}

  \item{\ldots}{not used currently.}
}

\details{
 If necessary, more details than the description above
}

\value{
Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
objects to See Also as %\code{\link{help}},
}

%\examples{}

\keyword{ kwd1 }
\keyword{ kwd2 }
