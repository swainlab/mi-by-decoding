# Functions for estimating mutual information (MI) by decoding
#
# Copyright (C) 2018 Gasper Tkacik and Peter Swain
# 
# If you publish results that make use this software or the mutual 
# information by decoding algorithm, please cite:
# Granados, A.A., Pietsch, J.M.J., Cepeda-Humerez, S.A., Farquhar, I.L.,
# Tkacik, G., and Swain, P.S. (2018) Distributed and dynamic intracellular
# organization of extracellular information. Proc Natl Acad Sci U S A.
# https://dx.doi.org/10.1073/pnas.1716659115
# 
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Helper functions {{{

calcnormparams <- function(tcdata,soft=TRUE,perfeature=TRUE){
  params <- list()
  alldata <- Reduce(rbind,tcdata)
  if(soft){
    if(perfeature){
      params$offset <- apply(alldata,2,mean)
      params$scale <- 2*apply(alldata,2,sd)
    } else {
      params$offset <- rep(mean(alldata),ncol(alldata))
      params$scale <- rep(2*sd(alldata),ncol(alldata))
    }
  } else {
    if(perfeature){
      minval <- apply(alldata,2,min)
      maxval <- apply(alldata,2,max)
    } else {
      minval <- rep(min(alldata),ncol(alldata))
      maxval <- rep(max(alldata),ncol(alldata))
    }
    params$offset <- (maxval+minval)/2
    params$scale <- (maxval-minval)/2
  }
  # If there are any columns with zero variance, set the scale to 1 to avoid
  # division by zero:
  params$scale[params$scale==0] <- 1
  return(params)
}

applynorm <- function(tcdata,params){
  lapply(tcdata,function(x){
         nr <- nrow(x)
         return((x-matrix(rep(params$offset,nr),nr,byrow=TRUE))
          /matrix(rep(params$scale,nr),nr,byrow=TRUE)) })
}

exclude.cols <- c('confM','MutInf','totalerrors','featrank')

max.has.bootstraps <- function(m,parnames=setdiff(colnames(m),exclude.cols),Nbootstraps=100){
  if(is.null(m)) return(FALSE)
  parcols <- m[,parnames,drop=FALSE]
  meanvals <- tapply(m$MutInf,parcols,mean,na.rm=TRUE)
  tapply(m$MutInf,parcols,length)[[which.max(meanvals)]]>=Nbootstraps
}

get.max.pars <- function(m,parnames=setdiff(colnames(m),exclude.cols),n=1,Nbootstraps=NULL){
  parcols <- m[,parnames,drop=FALSE]
  meanvals <- tapply(m$MutInf,parcols,mean,na.rm=TRUE)
  pars <- tapply(split(parcols,seq(nrow(parcols))),parcols,
                 function(x) x[[1,drop=FALSE]],simplify=FALSE)
  parinds <- rev(order(meanvals))[1:min(n,length(pars))]
  if(!is.null(Nbootstraps)){
    parlens <- tapply(m$MutInf,parcols,length)
    parinds <- parinds[parlens[parinds]<Nbootstraps]
  }
  Reduce(rbind,pars[parinds,drop=FALSE])
}

get.max.nboots <- function(m,parnames=setdiff(colnames(m),exclude.cols)){
  parcols <- m[,parnames,drop=FALSE]
  meanvals <- tapply(m$MutInf,parcols,mean,na.rm=TRUE)
  tapply(m$MutInf,parcols,length)[[which.max(meanvals)]]
}

get.opt.rows <- function(parcols,varcols){
  minf <- tapply(varcols$MutInf,parcols,mean,na.rm=TRUE)
  maxind <- arrayInd(which.max(minf),dim(minf))[1,]
  parnames <- names(dimnames(minf))
  opt.pars <- sapply(structure(1:length(maxind),names=parnames),
                     function(x) dimnames(minf)[[x]][[maxind[[x]]]])
  Reduce(function(x,y) x & y,
         lapply(1:length(maxind), function(p)
                as.factor(parcols[[parnames[p]]])==opt.pars[[p]]))
}

# }}}

# Classifier wrappers {{{

trainandtestfns <-
  list(svmlinear=structure(function(dtrain,ltrain,dtest,params,class.weights,...)
       {
         # Parameters: ncomponents, cost
         if(is.null(params$ncomponents)) params$ncomponents <- ncol(dtrain)
         if(is.null(params$cost)) params$cost <- 1
         mapply(function(k,cost){
                if(k>ncol(dtrain)) return(rep(NA,nrow(dtest)))
                m <- svm(dtrain[,1:k],ltrain,scale=FALSE,kernel='linear',
                         class.weights=class.weights,cost=cost,...)
                predict(m,dtest[,1:k]) },
                params$ncomponents, params$cost, SIMPLIFY=FALSE)
       }, default.params=list(ncomponents=1:10,cost=2^seq(-5,15,2))),

       svmrbf=structure(function(dtrain,ltrain,dtest,params,class.weights,...)
       {
         # Parameters: ncomponents, cost, gamma
         if(is.null(params$ncomponents)) params$ncomponents <- ncol(dtrain)
         if(is.null(params$cost)) params$cost <- 1
         if(is.null(params$gamma)) params$gamma <- 1/ncol(dtrain)
         mapply(function(k,cost,g){
                if(k>ncol(dtrain)) return(rep(NA,nrow(dtest)))
                m <- svm(dtrain[,1:k],ltrain,scale=TRUE,kernel='radial',
                         gamma=g,cost=cost,class.weights=class.weights,...)
                predict(m,dtest[,1:k]) },
                params$ncomponents, params$cost, params$gamma, SIMPLIFY=FALSE)
       }, default.params=list(ncomponents=1:10,cost=2^seq(-5,15,2),
                              gamma=2^seq(-15,3,2))),

       svmlinearprob=structure(function(dtrain,ltrain,dtest,params,class.weights,...)
       {
         # Parameters: ncomponents, cost
         if(is.null(params$ncomponents)) params$ncomponents <- ncol(dtrain)
         if(is.null(params$cost)) params$cost <- 1
         Nclasses <- length(levels(ltrain))
         mapply(function(k,cost){
                if(k>ncol(dtrain))
                  return(matrix(rep(NA,nrow(dtest)*Nclasses),nrow=nrow(dtest)))
                m <- svm(dtrain[,1:k],ltrain,scale=TRUE,kernel='linear',
                         class.weights=class.weights,cost=cost,probability=TRUE,...)
                attr(predict(m,dtest[,1:k],probability=TRUE),'probabilities') },
                params$ncomponents, params$cost, SIMPLIFY=FALSE)
       }, default.params=list(ncomponents=1:10,cost=2^seq(-5,15,2))),

       rforest=structure(function(dtrain,ltrain,dtest,params,class.weights,
                                  featrank=TRUE,...)
       {
         # Parameters: ntree, mtry
         if(is.null(params$ntree)) params$ntree <- 500
         if(is.null(params$mtry)) params$mtry <- floor(sqrt(ncol(dtrain)))
         mapply(function(ntree,mtry){
                if(mtry>ncol(dtrain)) return(rep(NA,nrow(dtest)))
                m <- randomForest::randomForest(dtrain,ltrain,ntree=ntree,mtry=mtry,
                                  importance=featrank,classwt=class.weights,...)
                structure(predict(m,dtest,type='response'),
                          featrank=if(featrank)
                            randomForest::importance(m,type=1) else NULL) },
                params$ntree, params$mtry, SIMPLIFY=FALSE)
       }, default.params=list(ntree=seq(50,500,50),mtry=c(1,2,4,8,12,16,20))),

       rforestprob=structure(function(dtrain,ltrain,dtest,params,class.weights,
                                      featrank=TRUE,...)
       {
         # Parameters: ntree, mtry
         if(is.null(params$ntree)) params$ntree <- 500
         if(is.null(params$mtry)) params$mtry <- floor(sqrt(ncol(dtrain)))
         Nclasses <- length(levels(ltrain))
         mapply(function(ntree,mtry){
                if(mtry>ncol(dtrain)) 
                  return(matrix(rep(NA,nrow(dtest)*Nclasses),nrow=nrow(dtest)))
                m <- randomForest::randomForest(dtrain,ltrain,ntree=ntree,mtry=mtry,
                                  importance=featrank,classwt=class.weights,...)
                structure(predict(m,dtest,type='prob'),
                          featrank=if(featrank)
                            randomForest::importance(m,type=1) else NULL) },
                params$ntree, params$mtry, SIMPLIFY=FALSE)
       }, default.params=list(ntree=seq(50,500,50),mtry=c(1,2,4,8,12,16,20))),

       xgboost=structure(function(dtrain,ltrain,dtest,params,class.weights,
                                  featrank=TRUE,...)
       {
         # Parameters: nrounds, max.depth
         if(is.null(params$nrounds)) params$nrounds <- 200
         if(is.null(params$max.depth)) params$max.depth <- 6
         trainpars <- unique(params[,c('max.depth')]) # Find unique rows
         Nclasses <- length(levels(ltrain))
         if(Nclasses==2){
           ms <- lapply(trainpars, function(max.depth)
                        xgboost::xgboost(data=dtrain,label=as.numeric(ltrain)-1,
                                nrounds=max(params$nrounds),max.depth=max.depth,
                                eta=0.1,subsample=0.5,verbose=0,save_period=NULL,
                                objective='binary:logistic',...))
           p <- function(m,k) as.numeric(predict(m,dtest,ntreelimit=k)>0.5)+1
         } else {
           ms <- lapply(trainpars, function(max.depth)
                        xgboost::xgboost(data=dtrain,label=as.numeric(ltrain)-1,
                                nrounds=max(params$nrounds),max.depth=max.depth,
                                eta=0.1,subsample=0.5,verbose=0,save_period=NULL,
                                objective='multi:softmax',num.class=Nclasses,...))
           p <- function(m,k) predict(m,dtest,ntreelimit=k)+1
         }
         if(featrank){
           features <- as.character(1:ncol(dtrain))
           featranks <- 
             lapply(ms, function(m){ f <- xgboost::xgb.importance(features,model=m)
                    unname(structure(f$Gain,names=f$Feature)[features]) })
         }
         lvls <- levels(ltrain)
         mapply(function(nrounds,max.depth){
                ind <- which(trainpars==max.depth)
                structure(factor(lvls[p(ms[[ind]],nrounds)],levels=lvls),
                          featrank=if(featrank) featranks[[ind]] else NULL) },
                params$nrounds, params$max.depth, SIMPLIFY=FALSE)
       }, default.params=list(nrounds=seq(20,200,20),max.depth=3:10)),

       xgboostprob=structure(function(dtrain,ltrain,dtest,params,class.weights,
                                      featrank=TRUE,...)
       {
         # Parameters: nrounds, max.depth
         if(is.null(params$nrounds)) params$nrounds <- 200
         if(is.null(params$max.depth)) params$max.depth <- 6
         trainpars <- unique(params[,c('max.depth')]) # Find unique rows
         Nclasses <- length(levels(ltrain))
         if(Nclasses==2){
           ms <- lapply(trainpars, function(max.depth)
                        xgboost::xgboost(data=dtrain,label=as.numeric(ltrain)-1,
                                nrounds=max(params$nrounds),max.depth=max.depth,
                                eta=0.1,subsample=0.5,verbose=0,save_period=NULL,
                                objective='binary:logistic',...))
           p <- function(m,k){ pred1 <- predict(m,dtest,ntreelimit=k); cbind(pred1,1-pred1) }
         } else {
           ms <- lapply(trainpars, function(max.depth)
                        xgboost::xgboost(data=dtrain,label=as.numeric(ltrain)-1,
                                nrounds=max(params$nrounds),max.depth=max.depth,
                                eta=0.1,subsample=0.5,verbose=0,save_period=NULL,
                                objective='multi:softprob',num.class=Nclasses))
           p <- function(m,k) matrix(predict(m,dtest,ntreelimit=k),ncol=Nclasses,byrow=TRUE)
         }
         if(featrank){
           features <- as.character(1:ncol(dtrain))
           featranks <-
             lapply(ms, function(m){ f <- xgboost::xgb.importance(features,model=m)
                    unname(structure(f$Gain,names=f$Feature)[features]) })
         }
         mapply(function(nrounds,max.depth){
                ind <- which(trainpars==max.depth)
                structure(p(ms[[ind]],nrounds),
                          featrank=if(featrank) featranks[[ind]] else NULL) },
                params$nrounds, params$max.depth, SIMPLIFY=FALSE)
       }, default.params=list(nrounds=seq(20,200,20),max.depth=3:10)))

confMbyDecoding <- function(traindata,testdata,params,method=MIdecodingClassifiers(),...){
  method <- match.arg(method)
  if(!is.data.frame(params))
    stop('"params" must be a "data.frame" specifying parameter combinations to try')
  default.params <- lapply(trainandtestfns,attr,which='default.params')[[method]]
  if(!all(names(params) %in% names(default.params)))
    warning('Not all specified "params" are used by the "',method,'" classifier')

  # Ensure that the data sets are named
  if(is.null(names(traindata))) names(traindata) <- paste0('c',1:length(traindata))
  if(is.null(names(testdata))) names(testdata) <- names(traindata)
  if(!all(names(testdata)==names(traindata)))
    stop('Test and training data sets must be named identically')

  trainlabel <- factor(Reduce(c,mapply(function(x,nm) rep(nm,nrow(x)),
                                       traindata,names(traindata),SIMPLIFY=FALSE)))
  testlabel <- factor(Reduce(c,mapply(function(x,nm) rep(nm,nrow(x)),
                                      testdata,names(testdata),SIMPLIFY=FALSE)))

  testpreds <-
    trainandtestfns[[method]](Reduce(rbind,lapply(traindata,cbind)),trainlabel,
                              Reduce(rbind,lapply(testdata,cbind)),params,
                              sapply(traindata,nrow)/length(trainlabel),...)

  # Remove any NA parameters
  na.rows <- sapply(testpreds,function(x) any(is.na(x)))
  params <- params[!na.rows,,drop=FALSE]
  testpreds <- testpreds[!na.rows]

  # Calculate confusion matrix
  params$confM <-
    lapply(testpreds, function(testpred){
           if(is.null(dim(testpred)))
             confM <- t(sapply(levels(testlabel),function(rl){
                               sapply(levels(testpred),function(est){
                                      sum(est==testpred[testlabel==rl]) }) }))
           else
             confM <- t(sapply(levels(testlabel),function(rl)
                               apply(testpred[testlabel==rl,],2,mean)))
           confM/sum(confM) })

  # Calculate MI from confusion matrix
  params$MutInf <- sapply(params$confM, Info)

  # Additional parameters depending on the chosen classifier
  if(all(sapply(testpreds,function(x) is.null(dim(x)))))
    params$totalerrors <- sapply(testpreds,function(x) sum(x!=testlabel))
  if(all(sapply(testpreds,function(x) !is.null(attr(x,'featrank')))))
    params$featrank <- lapply(testpreds,attr,which='featrank')

  params
}

# }}}

# Exported functions {{{

MIdecodingClassifiers <- function(withpars=FALSE){
  if(withpars) lapply(trainandtestfns, function(x) attr(x,'default.params'))
  else names(trainandtestfns)
}

Info <- function(M){
  if(any(is.na(M)) || any(M<0)) return(NA)
  
  # For numerical stability when probabilities approach zero (this is
  # valid by L'Hopital's rule for the calculation of mutual information):
  if(any(M==0)) M <- M + min(M[M!=0])*10^-8

  Mlen <- ncol(M)
  L <- log2(M)
  P1 <- apply(M,2,sum)
  P2 <- apply(M,1,sum)

  return(sum(M*L) - sum(M*log2(t(matrix(rep(P1,Mlen),Mlen))*matrix(rep(P2,Mlen),Mlen))))
}

MIdecoding <- function(tcdata, Nbootstraps=25, classifier=MIdecodingClassifiers(),
                       pca=classifier=='svmlinear', params=NULL,
                       normalise=c('bootstrap','none','raw'),
                       softnorm=TRUE, featurenorm=TRUE,
                       crossval=0, n.opt.pars=5,
                       Nbootstrap.batch=if(crossval>0) 12 else Nbootstraps,
                       silent=FALSE, parallel=TRUE, uniform=TRUE,
                       opt.par.only=TRUE,...){
  classifier <- match.arg(classifier)
  normalise <- match.arg(normalise)

  # Ensure the required classifier libraries are loaded
  if(grepl('^xgboost',classifier)){
    if(!requireNamespace('xgboost', quietly=TRUE)){
      warning('the "xgboost" package must be installed to use XGBoost ',
              'classification; defaulting to "svmrbf"...')
      classifier <- 'svmrbf'
    }
  } else if(grepl('^rforest',classifier)){
    if(!requireNamespace('randomForest', quietly=TRUE)){
      warning('the "randomForest" package must be installed to use Random ',
              'Forest classification; defaulting to "svmrbf"...')
      classifier <- 'svmrbf'
    }
  }

  Ntimepoints <- sapply(tcdata,ncol)
  if(all(Ntimepoints==Ntimepoints[1])) Ntimepoints <- Ntimepoints[1]
  else stop('All classes must have the same number of time points')

  Nclasses <- length(tcdata)
  classnames <- names(tcdata)

  default.params <- attr(trainandtestfns[[classifier]],'default.params')
  if(is.null(params)) params <- default.params
  else if(is.numeric(params))
    params <- structure(list(params),names=names(default.params)[1])
  # Ignore specification of Ncomponents if PCA is off
  if(!pca) params$Ncomponents <- NULL
  params <- expand.grid(params)
  init.params <- params
  parnames <- names(params)

  # Use only 90% of the available cells for each bootstrap
  Ncells <- sapply(tcdata, function(x) floor(0.9*nrow(x)))
  if(uniform) Ncells <- rep(min(Ncells),Nclasses)
  
  Ntest <- round(Ncells/4)
  Ntrain <- Ncells - Ntest

  if(normalise=='raw')
    tcdata <- applynorm(tcdata,calcnormparams(tcdata,softnorm,featurenorm))

  lapplyfn <-
    if(grepl('^xgboost',classifier) || !parallel) lapply else mclapply

  crossval.bootstraps <- NULL
  crossval.agg <- NULL

  if(crossval>0){
    # Generate random samples for the bootstraps outside of the parallel loop
    bootstrap.samples <- lapply(1:crossval,function(b)
                                mapply(function(x,nc) sample.int(nrow(x),nc),
                                       tcdata,Ncells,SIMPLIFY=FALSE))
    if(!silent) cat('Cross-validating')
    crossval.bootstraps <-
      lapply(bootstrap.samples,function(b){
             if(!silent) cat('.')
             bdata <- mapply(function(x,ss) cbind(x[ss,]),tcdata,b,SIMPLIFY=FALSE)

             # Split the data for cross-validation
             Nfolds <- 4
             Ncells.fold <-
               lapply(sapply(bdata,nrow),function(n)
                      c(0,cumsum(rep(floor(n/Nfolds),Nfolds)
                                 + c(rep(1,n%%Nfolds),rep(0,Nfolds-n%%Nfolds)))))

             folds <- 
               lapplyfn(1:Nfolds,function(f){
                        # Split the data into training and test sets
                        traindata <- mapply(function(x,n){
                                            xmask <- rep(TRUE,nrow(x))
                                            xmask[(n[f]+1):(n[f+1])] <- FALSE
                                            cbind(x[xmask,]) },
                                            bdata,Ncells.fold,SIMPLIFY=FALSE)
                        testdata <- mapply(function(x,n){
                                           xmask <- rep(FALSE,nrow(x))
                                           xmask[(n[f]+1):(n[f+1])] <- TRUE
                                           cbind(x[xmask,]) },
                                           bdata,Ncells.fold,SIMPLIFY=FALSE)

                        if(pca){
                          pcs <- prcomp(Reduce(rbind,traindata))
                          traindata <- lapply(traindata,function(x) predict(pcs,x))
                          testdata <- lapply(testdata,function(x) predict(pcs,x))
                        }

                        if(normalise=='bootstrap'){
                          normparams <- calcnormparams(traindata,softnorm,featurenorm)
                          traindata <- applynorm(traindata,normparams)
                          testdata <- applynorm(testdata,normparams)
                        }

                        confMbyDecoding(traindata,testdata,params,classifier,...) }) })
    crossval.bootstraps <- Reduce(c,crossval.bootstraps)

    validnames <- Reduce(intersect,lapply(crossval.bootstraps,names))
    crossval.agg <- Reduce(rbind,lapply(crossval.bootstraps,function(x) x[,validnames]))
    params <- get.max.pars(crossval.agg,parnames,n=n.opt.pars,Nbootstraps=Nbootstraps)

    if(!silent) cat('\n')
  }

  if(!silent) cat('Bootstrapping')

  bootstraps <- crossval.agg
  max.nboots <- 0
  while(!max.has.bootstraps(bootstraps,parnames,Nbootstraps)){
    if(!silent) cat('.')
    # Generate random samples for the bootstraps outside of the parallel loop
    bootstrap.samples <- lapply(1:Nbootstrap.batch,function(b)
                                mapply(function(x,nc) sample.int(nrow(x),nc),
                                       tcdata,Ncells,SIMPLIFY=FALSE))

    extra.bootstraps <-
      lapplyfn(bootstrap.samples,function(b){
               bdata <- mapply(function(x,ss) cbind(x[ss,]),tcdata,b,SIMPLIFY=FALSE)

               # Split the data into training and test sets
               traindata <- mapply(function(x,n) cbind(x[1:n,]),bdata,Ntrain,SIMPLIFY=FALSE)
               testdata <- mapply(function(x,nte,ntr) cbind(x[(ntr+1):(ntr+nte),]),
                                  bdata,Ntest,Ntrain,SIMPLIFY=FALSE)

               if(pca){
                 pcs <- prcomp(Reduce(rbind,traindata))
                 traindata <- lapply(traindata,function(x) predict(pcs,x))
                 testdata <- lapply(testdata,function(x) predict(pcs,x))
               }

               if(normalise=='bootstrap'){
                 normparams <- calcnormparams(traindata,softnorm,featurenorm)
                 traindata <- applynorm(traindata,normparams)
                 testdata <- applynorm(testdata,normparams)
               }

               confMbyDecoding(traindata,testdata,params,classifier,...) })

    validnames <- Reduce(intersect,lapply(extra.bootstraps,names))
    if(!is.null(bootstraps)) validnames <- intersect(names(bootstraps),validnames)
    bootstraps <- 
      rbind(bootstraps,Reduce(rbind,lapply(extra.bootstraps,function(x) x[,validnames])))
    if(crossval>0) params <- get.max.pars(bootstraps,parnames,
                                          n=n.opt.pars,Nbootstraps=Nbootstraps)
    
    extra.max <- get.max.nboots(bootstraps,parnames)
    if(extra.max>max.nboots){
      max.nboots <- extra.max
      if(!silent) cat(max.nboots)
    }
  }
  if(!silent) cat('\n')

  Nparams <- nrow(init.params)
  varnames <- setdiff(names(bootstraps),parnames)

  if(opt.par.only){
    parcols <- bootstraps[,parnames,drop=FALSE]
    varcols <- bootstraps[,varnames,drop=FALSE]
    opt.rows <- get.opt.rows(parcols,varcols)
    bootstraps <- varcols[opt.rows,,drop=FALSE]
    init.params <- parcols[opt.rows,,drop=FALSE][1,,drop=FALSE]
    parnames <- character(0)
  }

  arguments <- 
    list(Nbootstraps=Nbootstraps, classifier=classifier, pca=pca,
         normalise=normalise, softnorm=softnorm, featurenorm=featurenorm,
         crossval=crossval, n.opt.pars=n.opt.pars, classnames=classnames,
         Nbootstrap.batch=Nbootstrap.batch, uniform=uniform)
  structure(bootstraps, class=c('MIdecoding','data.frame'),
            parcols=parnames, varcols=varnames, params=init.params,
            Nparams=Nparams, Nclasses=Nclasses, Ntimepoints=Ntimepoints,
            Ntest=Ntest, Ntrain=Ntrain, arguments=arguments)
}

print.MIdecoding <- function(x,...){
  if(is.null(attr(x,'Nclasses'))) return(print.data.frame(x))
  a <- attr(x,'arguments')
  cat('"MIdecoding" result for MI between ',attr(x,'Nclasses'),
      ' states and time series of length ',attr(x,'Ntimepoints'),'\n\n',
      'Classifier: "',a$classifier,'"',
      if(a$pca) ' with PCA' else character(0),'\n',
      'Hyperparameters: ',attr(x,'Nparams'),' combinations considered',
      ' (best has >= ',a$Nbootstraps,' bootstraps)\n\n',
      sep='')
  invisible(x)
}

summary.MIdecoding <- function(object, ...){
  a <- attr(object,'arguments')
  varcols <- object[,attr(object,'varcols'),drop=FALSE]
  optparams <- attr(object,'params')

  if(nrow(optparams)>1){
    if(length(attr(object,'parcols'))==0L)
      stop('"MIdecoding" result has been corrupted ("parcols" attribute missing)')
    parcols <- object[,attr(object,'parcols'),drop=FALSE]
    opt.rows <- get.opt.rows(parcols,varcols)
    optparams <- parcols[opt.rows,,drop=FALSE][1,,drop=FALSE]
    varcols <- varcols[opt.rows,,drop=FALSE]
  }

  varcols <- as.list(varcols)
  varcols$confM <- abind(varcols$confM, along=3)
  dimnames(varcols$confM) <-
    list(actual=a$classnames, predicted=a$classnames,
         bootstraps=1:dim(varcols$confM)[3])
  structure(varcols, class='summary.MIdecoding', optparams=optparams,
            Ntest=attr(object,'Ntest'), Ntrain=attr(object,'Ntrain'),
            Nclasses=attr(object,'Nclasses'), Ntimepoints=attr(object,'Ntimepoints'),
            Nparams=attr(object,'Nparams'), arguments=attr(object,'arguments'))
}

print.summary.MIdecoding <- function(x,...){
  a <- attr(x,'arguments')
  mi.stats <- format(c(mean(x$MutInf),sd(x$MutInf)), nsmall=1, digits=2)
  cat('\n"MIdecoding" summary for MI between ',attr(x,'Nclasses'),
      ' states and time series of length ',attr(x,'Ntimepoints'),'\n\n',
      'Lower bound on MI:\n  ', mi.stats[[1]], ' +/- ', mi.stats[[2]],
      ' (mean +/- S.D.; ',length(x$MutInf),' bootstraps)\n\n',
      'Classifier: "',a$classifier,'"',
      if(a$pca) ' with PCA' else character(0),'\n',
      'Best hyperparameters (out of ',attr(x,'Nparams'),
      ' combinations considered):\n', sep='')
  optparams <- attr(x,'optparams')
  rownames(optparams) <- ' '
  print(optparams)
  cat('\nMedian total classification errors:\n  ', median(x$totalerrors),
      ' (out of ', sum(attr(x,'Ntest')),' samples in test set)\n\n',
      'Mean confusion matrix:\n',sep='')
  confM <- apply(x$confM, 1:2, mean, na.rm=TRUE)
  names(dimnames(confM))[1] <- '  actual'
  print(confM)
  cat('\n')
  invisible(x)
}

as.array.MIdecoding <- function(x, ...){
  if(length(attr(x,'parcols'))==0L)
    return(structure(x$MutInf,dim=length(x$MutInf),
                     dimnames=list(bootstraps=1:length(x$MutInf))))

  parcols <- x[,attr(x,'parcols'),drop=FALSE]

  Nparams <- ncol(attr(x,'params'))
  MaxBootstraps <- max(as.vector(tapply(x$MutInf,parcols,length)))
  ragged <- tapply(x$MutInf,parcols,function(p) p)
  a <- aperm(abind(lapply(1:MaxBootstraps,function(n)
                          array(sapply(ragged,function(p) p[n]),
                                dim=dim(ragged),dimnames=dimnames(ragged))),
                   along=Nparams+1,use.dnns=TRUE),c(Nparams+1,1:Nparams))
  names(dimnames(a))[[1]] <- 'bootstraps'
  return(a)
}

#}}}

# vim:set foldmethod=marker:
