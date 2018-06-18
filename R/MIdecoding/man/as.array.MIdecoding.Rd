\name{as.array.MIdecoding}
\alias{as.array.MIdecoding}

\title{Cast Bootstrap Estimates of Mutual Information into an Array}
\description{
  Extract an array of bootstrap estimates of Mutual Information (MI) organised
  by bootstrap and hyperparameter from the \code{data.frame} produced by
  \code{\link{MIdecoding}}.
}
\usage{
\method{as.array}{MIdecoding}(x, ...)
}

\arguments{
  \item{x}{an object of class \code{"MIdecoding"} as produced by
    \code{\link{MIdecoding}}.}
  \item{...}{additional arguments passed to or from previous methods.}
}

\note{
  This function is designed for use with the full output from
  \code{\link{MIdecoding}}, i.e., by specifying \code{opt.par.only=FALSE} in
  the arguments to \code{\link{MIdecoding}}.
}

\value{
  A numeric array of MI estimates (in bits) where the first dimension
  specifies bootstrap replicate and where remaining dimensions specify
  estimates for different classifier hyperparameters. Dimensions are named
  according to hyperparameter name, and elements along each dimension are
  named according to hyperparameter value.  Missing bootstraps are indicated
  by \code{NA} (e.g., as produced by non-zero \code{crossval} arguments in the
  call to \code{\link{MIdecoding}}).
}

\references{
  Granados, A. A., Pietsch, J. M. J., Cepeda-Humerez, S. A., Farquhar, I. L.,
  Tkacik, G., Swain, P. S. (2018) Distributed and dynamic intracellular
  organization of extracellular information. \emph{Proc Natl Acad Sci U S A}.
  \url{https://dx.doi.org/10.1073/pnas.1716659115}
}

\seealso{
  \code{\link{MIdecoding}} for generating bootstrap estimates of mutual
  information by decoding.
}
\examples{
## Load sample data
data(YeastStressTypeResponse)

## Obtain and summarise bootstrap estimates of MI by decoding
MIest <- MIdecoding(YeastStressTypeResponse,
                    params=list(ncomponents=1:10,
                                cost=c(0.1, 1, 10, 100)),
                    crossval=1, opt.par.only=FALSE)

## Cast to an array
MIest.array <- as.array(MIest)

## Mean MI obtained for each combination of hyperparameters
apply(MIest.array, 2:length(dim(MIest.array)), mean, na.rm=TRUE)
}

\keyword{htest}
\keyword{classif}
