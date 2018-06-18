\name{MIdecoding}
\alias{MIdecoding}
\alias{print.MIdecoding}

\title{Estimate Mutual Information by Decoding}
\description{
  Generate bootstrap estimates for a lower bound on the Mutual Information
  (MI) between discrete states and stochastic time series by evaluating the
  performance of a classifier function (machine learning algorithm).
}
\usage{
MIdecoding(tcdata, Nbootstraps=25, classifier=MIdecodingClassifiers(),
           pca=classifier=="svmlinear", params=NULL,
           normalise=c("bootstrap", "none", "raw"),
           softnorm=TRUE, featurenorm=TRUE,
           crossval=0, n.opt.pars=5,
           Nbootstrap.batch=if(crossval>0) 12 else Nbootstraps,
           silent=FALSE, parallel=TRUE, uniform=TRUE, 
           opt.par.only=TRUE, ...)
\method{print}{MIdecoding}(x, ...)
}

\arguments{
  \item{tcdata}{
    a list of matrices, the stochastic time series for each state, where each
    row is a sample and columns are time points. All states must have the same
    number of time points (columns), but may differ in the number of samples.
  }
  \item{Nbootstraps}{the number of bootstrap estimates to make.}
  \item{classifier}{
    the classifier function to use. Defaults to \code{"svmlinear"}, a linear
    SVM classifier using the one-vs-one method to support multi-class
    classification (i.e., when number of classes is more than two). Other
    possible choices include \code{"svmrbf"} for a SVM classifier with the
    radial basis function kernel, \code{"rforest"} for the random forest
    classifier, and \code{"xgboost"} for the eXtreme Gradient Boosting
    classifier. For more information see the \sQuote{Details} section in
    \code{\link{MIdecodingClassifiers}}.
  }
  \item{pca}{
    whether or not the data should be transformed onto principal components
    before application of the classifier. This allows reduction of
    dimensionality for the SVM classifiers by specification of the
    \code{Ncomponents} parameter in \code{params}. Primarily intended for use
    with the \code{"svmlinear"} classifier. Components are determined by
    performing PCA across all training data.
  }
  \item{params}{
    a named list of hyperparameter values to try. If \code{NULL}, then a
    classifier-specific default list of hyperparameters is used. For more
    information see the \sQuote{Details} section of
    \code{\link{MIdecodingClassifiers}}.
  }
  \item{normalise}{
    whether and when data should be normalised. \code{"bootstrap"} (the
    default) re-normalises the data set for each bootstrap; \code{"none"}
    turns normalisation off; \code{"raw"} normalises the entire data set
    before bootstrapping.  Normalisation can improve classifier performance
    and is applied equally across all classes and across training and testing
    sets.
  }
  \item{softnorm}{
    whether normalisation should be hard (such that all data lies between -1
    and +1) or soft (subtracting the mean and dividing by two times the
    standard deviation).
  }
  \item{featurenorm}{
    whether normalisation should be performed feature-wise (i.e., per time
    point) or according to the aggregate of all features.
  }
  \item{crossval}{
    the number of rounds of four-fold cross-validation to run before picking
    optimal hyperparameters. If \code{crossval==0}, then cross-validation is
    not performed and the specified number of bootstraps are obtained for each
    combination of hyperparameters. Results from cross-validation bootstraps
    form part of the final output.
  }
  \item{n.opt.pars}{
    the number of optimal hyperparameter combinations for which the full
    complement of bootstraps should be obtained. Ignored if cross-validation
    is not performed.
  }
  \item{Nbootstrap.batch}{
    the number of bootstraps to calculate in each batch. After each batch, the
    hyperparameters that are optimal will be re-evaluated.
  }
  \item{silent}{whether to turn off output of progress.}
  \item{parallel}{whether to run bootstraps in parallel.}
  \item{uniform}{
    whether to train and test using equal numbers of samples in each class. If
    \code{FALSE} (not recommended), the number of samples per class depends on
    availability, and affects training weights and class priors in estimates
    of MI from test data.
  }
  \item{opt.par.only}{
    whether to subset final results to just the bootstraps of the
    highest-ranking combination of hyperparameters.
  }
  \item{...}{
    additional arguments to be passed on to the classifier. For more
    information see \sQuote{Details}.
  }
  \item{x}{an object of class \code{"MIdecoding"}}
}

\details{
  Given small sample sizes, it is difficult to estimate the Mutual Information
  (MI) between discrete states and stochastic time series due to the high
  dimensionality (large configuration space) of the time series. The MI by
  decoding algorithm attempts to overcome this limitation by performing the
  estimation in state space. A machine learning classifier (the
  \sQuote{decoder}) is trained on a subset of the data to transform time
  series into state space. By then treating confusion matrices derived from
  test data as providing joint probabilities over the state-space probability
  distributions, we can calculate MI. The data processing inequality then
  guarantees that this is a lower bound to the true MI between the states and
  time series. As such, the higher the performance of the trained classifier
  function, the tighter the bound becomes.

  In order to find the optimal classifier function, \code{MIdecoding} provides
  several alternative machine learning algorithms (described in
  \code{\link{MIdecodingClassifiers}}) and performs a heuristic search over
  pre-defined values of the various hyperparameters available to each of the
  algorithms. When values for multiple hyperparameters are specified in the
  \code{params} argument (as they are by default), a grid search over all
  possible combinations of hyperparameters is performed. If this space is
  large, or a high number of bootstraps is desired, then it is recommended to
  specify one or more rounds of cross validation (i.e., \code{crossval > 0}),
  in which case, bootstraps at all possible hyperparameter combinations are
  only obtained in these initial rounds. After that, additional bootstraps are
  obtained only for the \code{n.opt.pars} hyperparameter combinations with
  highest MI. Every \code{Nbootstrap.batch} bootstraps, the optimal
  hyperparameters are reevaluated, and bootstraps continue to be obtained
  until one combination of hyperparameters has at least \code{Nbootstrap}
  bootstrap replicates and has mean MI greater than any other combination of
  hyperparameters.

  If either the \code{"xgboost"} or \code{"rforest"} classifiers are chosen,
  then by default they additionally return the importance of each feature
  (e.g., each time point) in the classifier. If this information is not
  required, then some processing time can be saved by specifying
  \code{featrank=FALSE} as an additional argument (which gets passed on to the
  classifiers).
}

\value{
  An object of class \code{MIdecoding}, which is a \code{data.frame} with rows
  for each bootstrap, columns specifying the hyperparameter values used at
  the given bootstrap, and columns with:

  \item{MutInf}{the mutual information estimated from this bootstrap,}
  \item{confM}{the (normalised) confusion matrix for this bootstrap,}
  \item{totalerrors}{the total number of errors made by the classifier across
    all classes for this bootstrap, and}
  \item{featrank}{(for \code{"xgboost"} and \code{"rforest"} classifiers only)
    the importance of each feature.}
}

\references{
  Granados, A. A., Pietsch, J. M. J., Cepeda-Humerez, S. A., Farquhar, I. L.,
  Tkacik, G., Swain, P. S. (2018) Distributed and dynamic intracellular
  organization of extracellular information. \emph{Proc Natl Acad Sci U S A}.
  \url{https://dx.doi.org/10.1073/pnas.1716659115}
}

\seealso{
  \code{\link{Info}} for calculating mutual information from a confusion
  matrix, \code{\link{MIdecodingClassifiers}} for details on the classifiers
  that can be used with \code{MIdecoding} and their hyperparameters,
  \code{\link{summary.MIdecoding}} to obtain reformatted results for just the
  optimal set of hyperparameters, and \code{\link{as.array.MIdecoding}} for
  investigating classifier performance as a function of the hyperparameters.
}
\examples{
## Load sample data
data(YeastStressTypeResponse)

## Obtain bootstrap estimates of MI by decoding
MIest <- MIdecoding(YeastStressTypeResponse,
                    params=list(ncomponents=1:10))
MIest

## If the data set has high-order structure that can only be
## resolved with non-linear transformation, then it is more
## appropriate to choose a non-linear classifier (in which case
## pca=FALSE by default):
MIest.rbf <- MIdecoding(YeastStressTypeResponse,
                        classifier='svmrbf',
                        params=list(cost=100, gamma=10^(-4:-1)))
## Compare with linear SVM:
mean(MIest.rbf$MutInf)
mean(MIest$MutInf)

## Perform a wider search of parameter space, but use
## cross-validation to limit the number of bootstraps for
## poorly-performing parameter combinations, and include results
## for all parameters in output:
MIest.scan <- MIdecoding(YeastStressTypeResponse,
                         params=list(ncomponents=1:10,
                                     cost=c(0.1, 1, 10)),
                         crossval=1, opt.par.only=FALSE)
## Display the hyperparameters that were searched over
lapply(attr(MIest.scan, 'params'), unique)
}

\keyword{htest}
\keyword{classif}
\keyword{optimize}
