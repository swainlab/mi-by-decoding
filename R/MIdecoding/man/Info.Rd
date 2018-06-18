\name{Info}
\alias{Info}

\title{Calculate Mutual Information from a Joint Probability Matrix}
\description{
  Calculate Mutual Information (MI) between, for example, the actual and
  predicted classes output by a machine learning classifier.
}
\usage{
Info(M)
}

\arguments{\item{M}{a unit-normalised square matrix}}

\note{
  This function expects a unit-normalised joint probability matrix, so
  unnormalised confusion matrices will not produce valid MI estimates.
}

\value{The mutual information (in bits) between joint discrete probability
  distributions.} 

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
## MI between uncorrelated two-state distributions
Info(matrix(rep(0.25, 4), 2)) # 0 bits

## MI between perfectly correlated two-state distributions
Info(matrix(c(0.5, 0, 0, 0.5), 2)) # 1 bit
}

\keyword{htest}
