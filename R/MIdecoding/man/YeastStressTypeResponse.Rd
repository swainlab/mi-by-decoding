\name{YeastStressTypeResponse}
\docType{data}
\alias{YeastStressTypeResponse}
\title{Yeast Msn2 Localisation in Response to Stress}
\description{
  Single-cell microscopy data from a microfluidics-based screening of the
  dynamics of nuclear translocation of Msn2, a transcription factor in
  \emph{S. cerevisiae} that is the target of evolutionarily conserved
  signalling pathways including the protein kinase A (PKA) and TOR kinase
  pathways. The data set includes single-cell time series describing the
  response for environmental transitions from rich media (2\% glucose) into
  three different stresses: carbon stress (low glucose), osmotic stress or
  oxidative stress. 
}
\usage{data(YeastStressTypeResponse)}
\format{
  A list of matrices named according to environmental condition (\code{"gluc"}
  for carbon stress, \code{"nacl"} for osmotic stress, \code{"oxid"} for
  oxidative stress, and \code{"rich"} for no stress), with the rows of each
  matrix corresponding to independent single-cell time series, and the columns
  corresponding to 21 observations taken at 2.5 minute intervals.
}
\note{
  This experimental data set forms a small subset of the complete data
  available at \url{https://dx.doi.org/10.7488/ds/2214}.
}
\source{\url{https://dx.doi.org/10.7488/ds/2214}}
\references{
  Granados, A. A., Pietsch, J. M. J., Cepeda-Humerez, S. A., Farquhar, I. L.,
  Tkacik, G., Swain, P. S. (2018) Distributed and dynamic intracellular
  organization of extracellular information. \emph{Proc Natl Acad Sci U S A}.
  \url{https://dx.doi.org/10.1073/pnas.1716659115}
}
\keyword{datasets}
