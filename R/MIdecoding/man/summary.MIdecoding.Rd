\name{summary.MIdecoding}
\alias{summary.MIdecoding}
\alias{print.summary.MIdecoding}

\title{Summarise Bootstrap Estimates of Mutual Information}
\description{
  Extract and print summary statistics for the bootstrap estimates of Mutual
  Information (MI) produced by \code{\link{MIdecoding}} at the combination of
  hyperparameters resulting in the maximum MI.
}
\usage{
\method{summary}{MIdecoding}(object, ...)
\method{print}{summary.MIdecoding}(x, ...)
}

\arguments{
  \item{object}{an object of class \code{"MIdecoding"} as produced by
    \code{\link{MIdecoding}}.}
  \item{x}{an object of class \code{"summary.MIdecoding"} as produced by
    \code{summary.MIdecoding}.}
  \item{...}{additional arguments passed to or from previous methods.}
}

\details{
  The print method displays the mean and standard deviation of the bootstrap
  MI estimates, the hyperparameter values giving the highest MI, the median
  number of classification errors made by the chosen classifier, and the mean
  confusion matrix.

  If a named list was used in the call to \code{\link{MIdecoding}}, then the
  dimensions of the confusion matrix array will be appropriately named by
  class, and, in any case, each dimension is labelled according to whether it
  corresponds to the actual class (\code{"actual"}), the predicted class
  (\code{"predicted"}) or to the bootstraps (\code{"bootstraps"}).
}

\value{
  An object of class \code{summary.MIdecoding}, which is a \code{list} with
  elements:

  \item{MutInf}{the bootstrap MI estimates,}
  \item{confM}{the (normalised) confusion matrices corresponding to the MI
    estimates as a 3-dimensional array with bootstraps along the third
    dimension, and}
  \item{totalerrors}{the (absolute) number of errors made by the classifier
    for each bootstrap.}
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
                    params=list(ncomponents=1:10),
                    opt.par.only=FALSE)
summary(MIest)

## Calculate the standard deviation of elements in the confusion matrix
MIest.summary <- summary(MIest)
apply(MIest.summary$confM, 1:2, sd)
}

\keyword{htest}
\keyword{classif}
