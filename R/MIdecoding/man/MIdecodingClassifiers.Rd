\name{MIdecodingClassifiers}
\alias{MIdecodingClassifiers}

\title{Decoding Classifiers Available in the MIdecoding Function}
\description{
  Outputs the machine learning classifiers available for estimation of Mutual
  Information (MI) by the \code{\link{MIdecoding}} function. Optionally
  outputs the default sets of hyperparameters used for each classifier.
}
\usage{MIdecodingClassifiers(withpars=FALSE)}

\arguments{
  \item{withpars}{
    whether to return for each available classifier the parameter values that
    define its default hyperparameter search space.
  }
}

\details{
  The best choice of classifier to use with the \code{\link{MIdecoding}}
  function is data-dependent, so it is worth performing preliminary tests with
  multiple classifiers to see which produces the highest estimated mutual
  information (MI), and thus a tighter lower bound on the true MI. 

  For many standard data sets, classification by a linear Support Vector
  Machine (SVM) on data transformed by Principal Component Analysis (PCA) will
  be sufficient, and can even outperform more advanced machine learning
  algorithms if the number of available samples is small. However, if a
  non-linear transformation would be required to successfully partition the
  data (e.g., if separation of the classes is in the frequency domain), then
  it will be necessary to use either SVM with the Radial Basis Function (RBF)
  kernel or one of the ensemble classifiers (Random Forest or eXtreme Gradent
  Boosting (XGBoost)).

  The classifiers currently available with \code{\link{MIdecoding}} and their
  associated hyperparameters are briefly introduced in the following.

  \subsection{Linear SVM}{
    The default classifier in the call to \code{\link{MIdecoding}} (specified
    as \code{classifier="svmlinear"}).

    SVM classification (with a linear kernel) using the one-vs-one method to
    support multi-class classification (i.e., when number of classes is more
    than two). Implementation is via the \code{\link[e1071]{svm}} function of
    the \pkg{e1071} package (which should already be installed since it is a
    requirement of this package). Further details of the algorithm can be
    found by consulting the documentation for the \code{\link[e1071]{svm}}
    function.

    Available hyperparameters:
    \describe{
      \item{\code{cost}}{cost of constraints violation.}
      \item{\code{ncomponents}}{
        (only if \code{pca=TRUE}) the number of principle components to use
        as input to the classifier.
      }
    }
  }

  \subsection{Nonlinear SVM}{
    Specified as \code{classifier="svmrbf"} in the call to
    \code{\link{MIdecoding}}.

    SVM classification with a RBF kernel using the one-vs-one method to support
    multi-class classification (i.e., when number of classes is more than two).
    Implementation is via the \code{\link[e1071]{svm}} function of the
    \pkg{e1071} package (which should already be installed since it is a
    requirement of this package). Further details of the algorithm can be found
    by consulting the documentation for the \code{\link[e1071]{svm}} function.

    Available hyperparameters:
    \describe{
      \item{\code{cost}}{cost of constraints violation.}
      \item{\code{gamma}}{\eqn{\gamma}{gamma} parameter of the RBF kernel.}
      \item{\code{ncomponents}}{
        (only if \code{pca=TRUE}) the number of principle components to use
        as input to the classifier.
      }
    }
  }

  \subsection{Random Forest}{
    Specified as \code{classifier="rforest"} in the call to
    \code{\link{MIdecoding}}.

    The Random Forest method of ensemble classification, with implementation via
    the \code{\link[randomForest]{randomForest}} function of the
    \pkg{randomForest} package (which needs to be installed in order to use this
    classifier). Further details of the algorithm can be found by consulting the
    documentation for the \code{\link[randomForest]{randomForest}} function.

    Available hyperparameters:
    \describe{
      \item{\code{ntree}}{number of trees grown.}
      \item{\code{mtry}}{number of predictors sampled for spliting at each
        node.}
    }
  }

  \subsection{XGBoost}{
    Specified as \code{classifier="xgboost"} in the call to
    \code{\link{MIdecoding}}.

    The eXtreme Gradient Boosting (XGBoost) method of ensemble classification,
    with implementation via the \code{\link[xgboost]{xgboost}} function of the
    \pkg{xgboost} package (which needs to be installed in order to use this
    classifier). Further details of the algorithm can be found by consulting the
    documentation for the \code{\link[xgboost]{xgboost}} function.

    Available hyperparameters:
    \describe{
      \item{\code{nrounds}}{maximum number of boosting iterations.}
      \item{\code{max.depth}}{maximum depth of a tree.}
    }
  }
}

\value{
  a character vector with the names of all available classifiers, or, if
  \code{withpars=TRUE}, a named list of lists of numeric vectors specifying,
  for each available classifier, the hyperparameter values scanned by default.
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
## List default hyperparameter values scanned by each classifier
MIdecodingClassifiers(TRUE)
}

\keyword{classif}
