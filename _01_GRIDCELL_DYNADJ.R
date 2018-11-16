# ---------------------------------------------------------
# Dynamical Adjustment: Function repository for Prediction step in dynamical adjustment
# Based on the Elastic Net
# ---------------------------------------------------------

# Sebastian Sippel 
# 24.10.2018

require(xts)
require(zoo)
require(raster)
require(ncdf4)


## NEXT STEPS TO DO FOR METHOD DEVELOPMENT: 
# - (1) Option to cross-validate alpha
#         #test = lapply(X = seq(0, 1, 0.1), FUN=function(cur.alpha) cv.glmnet(x = X[lag0.train.idx,], y = as.numeric(Y.xts[lag0.train.idx]), nfolds = nfolds, foldid = crossclass, alpha = cur.alpha, lambda.min.ratio = 0.1))
#          #plot(sapply(1:10, FUN=function(i) min(test[[i]]$cvm)))
# - (2) Fix folds -> DONE!
# - (3) Implement lag training steps for daily/monthly
# - (4) At some point: Transferability from one model to another (but: later)



### TRAINING AND PREDICTION OF TIME SERIES (with same dataset):
#' @title Train a model for a climatic target variable Y based on circulation pattern X
#' @description 
#' @param Y.train time series of (univariate) target variable to train on (an XTS object)
#' @param X Circulation patterns predict
#' @param X.train Circulation patterns to train on (if NULL, X is used). If specified, MUST be of same length as train.years in Y.train
#' @param X.date.seq Sequence of Dates in X (which will be used for forward prediction)
#' @param train.years subset of years to train on for dynamical adjustment (if NULL: all years are used)
## #' @param na.rm Automatically remove rows (from training) if there is an NA in Y?
#' @details
#' Some details about the function
#' @return A deseasonalized (daily or monthly) time series in xts format.
#' @references None
#' @author Sebastian Sippel
#' @examples
train.dyn.adj.elastnet.annual <- function(Y.train, X, X.train = NULL, X.date.seq = NULL, 
                                          train.years = NULL, train.months = 1:12, add.mon = 2,
                                          alpha = 1, nfolds = 10, ret.cv.model=F, s = "lambda.min",
                                          lags.to.aggregate = list(1:2, 3:7, c(8:30)), n.pc = 10) {
  
  ## (0) SANITY CHECKS AND TRANSFORMATION ON INPUT DATA:
  ## ---------------------------------------------------
  {
    # CHECK IF Separate dataset of TRAIN Years are given (can EITHER GIVE train.years OR train on entire Y time):
    if (is.null(train.years)) train.years = as.numeric(unique(format(as.Date(time(Y.train)), "%Y")))
    
    # CHECK IF X is a raster and transform to matrix if necessary:
    if (!is.matrix(X)) {
       # area.vec = values(area(X)) / max(values(area(X)))  # -> EOF in pred. field could be weighted by area -> but rather not..
       if (is.null(X.date.seq)) X.date.seq = as.Date(substring(text = names(X), first = 2, last = 11), format = "%Y.%m.%d")
       X = t(values(X))
    }
    # CHECK if X is matrix and X.date.seq can be taken from matrix rownames:
    if (is.matrix(X) & is.null(X.date.seq)) X.date.seq = as.Date(substring(text = rownames(X), first = 2, last = 11), format = "%Y.%m.%d")
    
    # CHECK if X.train is given; and whether it is raster?
    if (!is.null(X.train) & !is.matrix(X.train)) X.train = t(values(X.train))
    # if (is.null(X.train)) X.train = X[which(X.date.seq %in% as.Date(time(Y.train))),]  ## DERIVE X.train from X
        
    # CHECK FOR NA's in ENTIRE COLUMNS in X (and remove if necessary):  
    # na.idx.col = apply(X = X, MARGIN=2, FUN=function(x) all(is.na(x)))
    # if( any(na.idx.col) ) X = X[,!na.idx.col]
    
    if (train.months[1] == "seasonal") train.months = c(1, 4, 7, 10)
  }
  
  
  ## (1) INCLUDE LAG PREDICTORS IF DESIRED:
  ## ------------------------------------------------
  {
  if (!is.null(lags.to.aggregate) & is.null(X.train)) {
    X = add.lag.predictors(X = X, X.train = X.train, 
                           svd.idx = which(as.numeric(format(time(Y.train), "%Y")) %in% train.years), 
                           lags.to.aggregate = lags.to.aggregate, n.pc = n.pc)
  } else if (!is.null(lags.to.aggregate) & !is.null(X.train)) {
    ret.list = add.lag.predictors(X = X, X.train = X.train, 
                           svd.idx = which(as.numeric(format(time(Y.train), "%Y")) %in% train.years), 
                           lags.to.aggregate = lags.to.aggregate, n.pc = n.pc)
    X = ret.list$X
    X.train = ret.list$X.train
  }
  }
  
  ## (2) TRAINING & PREDICTION STEPS:
  ## --------------------------------
  cur.glmnet = list()   # DEFINE FRAMEWORK FOR MODEL:
  Yhat = rep(NA, length(X.date.seq)) # xts(x = NULL, order.by = time(Y.train))
  
  for (i in train.months) {
    # print(i)
    
    ## TRAINING IDX:
    cur.mon = (i-add.mon):(i+add.mon)
    if (add.mon == 1) {
      if (i == 1) cur.mon = c(12,1,2)
      if (i == 12) cur.mon = c(11,12,1)
    } else if (add.mon == 2) {
      if (i == 1) cur.mon = c(11,12,1,2,3)
      if (i == 2) cur.mon = c(12,1,2,3,4)
      if (i == 11) cur.mon = c(9,10,11,12,1)
      if (i == 12) cur.mon = c(10,11,12,1,2)
    }
    Y.train.idx = which(as.numeric(format(as.Date(time(Y.train)), "%m")) %in% cur.mon & 
                             as.numeric(format(as.Date(time(Y.train)), "%Y")) %in% train.years)    # WHICH INDICES ARE IN PRESENT MONTHS AND IN train.years?
    if (any(is.na(Y.train[Y.train.idx]))) Y.train.idx=Y.train.idx[-which(is.na(Y.train[Y.train.idx]))]  # is any Y value NA?
    if (!is.null(lags.to.aggregate)) if (any((1:max(unlist(lags.to.aggregate))) %in% Y.train.idx)) Y.train.idx=Y.train.idx[-which(Y.train.idx %in% (1:max(unlist(lags.to.aggregate))))]
    

    ## (2.1) GLMNET TRAINING BASED ON CROSS-VALIDATION:
    crossclass = c(rep.row(x = 1:nfolds, n = ceiling(length(Y.train.idx)/10)))[1:length(Y.train.idx)]
    
    if (is.null(X.train)) {
      cur.glmnet[[i]] = cv.glmnet(x = X[Y.train.idx,], y = as.numeric(Y.train[Y.train.idx]), nfolds = nfolds, foldid = crossclass, alpha = alpha)
    } else {
      cur.glmnet[[i]] = cv.glmnet(x = X.train[Y.train.idx,], y = as.numeric(Y.train[Y.train.idx]), nfolds = nfolds, foldid = crossclass, alpha = alpha)
    }
    
    ## (2.2) PREDICTION STEP FOR RESPECTIVE MONTH:
    if ( all(train.months == c(1, 4, 7, 10)) ) {   # "SEASONAL" PREDICTION
      X.pred.idx = which(as.numeric(format(X.date.seq, "%m")) %in% cur.mon)    ## PREDICTION IDX: idx to fill with predicted value in X
    } else {
      X.pred.idx = which(as.numeric(format(X.date.seq, "%m")) == i)    ## PREDICTION IDX: idx to fill with predicted value in X
    }
    Yhat[X.pred.idx] = predict.cv.glmnet(cur.glmnet[[i]], newx = X[X.pred.idx,], s = s)  
  }
  
  
  ## (3) Return different quantities:
  ## --------------------------------
  if (ret.cv.model) return(cur.glmnet)
  
  Yhat.xts = xts(x = Yhat, order.by = X.date.seq)
  return(Yhat.xts)
}






#' @title Add lag predictors to circulation matrix X and X.train
add.lag.predictors <- function(X, X.train = NULL, svd.idx = NULL, lags.to.aggregate = list(1:2, 3:7, c(8:30)), n.pc = 10) {
  
  # svd.idx = which(as.numeric(format(time(Y.train), "%Y")) %in% train.years)
  
  # (0) RUN Singular value decomposition on the whole matrix:
  if (is.null(X.train)) {
    X.svd = svd( t(X[svd.idx,]) )
  } else {
    X.svd = svd( t(X.train[svd.idx,]) )
  }
  
  # (1) PROJECT anomaly pattern on entire X:
  X.pc = (t(X.svd$u) %*% t(X))[1:n.pc,] # first 10 PC's as lag predictors (!!)
  if (!is.null(X.train)) X.train.pc = (t(X.svd$u) %*% t(X.train))[1:n.pc,]
  # plot(X.svd$v[,1], X.pc[1,1:7305])
  
  # (2) EXTRACT averages of past lags:
  test = lapply(X = lags.to.aggregate, FUN=function(cur.lag) {
    X.lag = sapply(X = 1:n.pc, FUN=function(idx) c(rep(NA, cur.lag[1]), rollmean(x = X.pc[idx,], k = length(cur.lag), na.pad = T, align = "right")))[1:(dim(X)[1]),]
    return(X.lag)
  })
  X.lag = do.call(cbind, test)
  
  if (!is.null(X.train)) {
    test = lapply(X = lags.to.aggregate, FUN=function(cur.lag) {
      X.lag = sapply(X = 1:n.pc, FUN=function(idx) c(rep(NA, cur.lag[1]), rollmean(x = X.train.pc[idx,], k = length(cur.lag), na.pad = T, align = "right")))[1:(dim(X.train)[1]),]
      return(X.lag)
    })
    X.train.lag = do.call(cbind, test)
  }
  
  # (3) Return lag values for X and X.train
  if (is.null(X.train)) {
    return( cbind(X, X.lag) )
  } else {
    return(list(X = cbind(X, X.lag), X.train = cbind(X.train, X.train.lag)))
  }
}


