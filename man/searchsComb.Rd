\name{searchsComb}

\encoding{latin1}
\alias{searchsComb}

\title{
Searchs the Element of a Combination
}

\description{
Searchs the element of a combination from all combinations
of the supplied vectors as created by
\code{\link{expand.grid}} without creating the combinations.
This function is called internally by \code{\link{MIMFA}}
and is not usually called directly by a user.
}

\usage{
searchsComb(args, idx)
}

\arguments{
  \item{args}{list containing the vectors sources for combinations.}

  \item{idx}{the index of the searched combination.}
}

\value{
The combination corresponding to \code{idx}.
}

\author{
Ignacio González
}

\examples{
comb <- expand.grid(height = seq(60, 80, 5), weight = seq(100, 300, 50),
                    sex = c("Male","Female"))
idx <- sample(nrow(comb), 1)
comb[idx, ]

searchsComb(list(height = seq(60, 80, 5), weight = seq(100, 300, 50),
            sex = c("Male","Female")), idx)
}

\keyword{utilities}
