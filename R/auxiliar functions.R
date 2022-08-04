#' Local linear regression with stepwise forward and maximum number of features

locstep <- function(arg, x, y, h, dmax, varimp) #x must be a matrix
{
  d <- length(arg)
  n <- length(y)
  est <- rep(0,d)
  if (d > 1) {
    ker <- function(xx) {return((2 * pi)^(-d/2) * exp(-rowSums(xx^2)/2))}
    argu <- matrix(arg, dim(x)[1], d, byrow = TRUE)
    w <- ker(sqrt(varimp)*(x - argu)/h)/h^d
    weights <- w/sum(w)
    
    models <- leaps::regsubsets(y ~ ., data=data.frame(y, x), weights = weights, method="forward", nvmax=dmax)
    est0 <- coef(models, which.min(summary(models)$bic))[-1]
    est[match(names(est0), colnames(data.frame(x)))] = est0
    est <- c(est, coef(models, which.min(summary(models)$bic))[1] )
  } else {
    ker <- function(xx) {return(exp(-xx^2/2))}
    x <- matrix(x, length(x), 1)
    w <- ker(sqrt(varimp)*(x - arg)/h)/h^d
    weights <- w/sum(w)
    models <- lm(y ~ ., data=data.frame(y, x), weights = weights)
    est <- c(coef(models)[-1], coef(models)[1] )
  }
  return(est)
}

#' Local variable importance using Random Forest

varimpcal <- function(x, y){
  rf <- randomForest::randomForest(y ~ ., data=data.frame(y,x), ntree = 1000, keep.forest = FALSE, localImp=TRUE)
  localimp <- abs(t(rf$localImportance))/rowSums(abs(t(rf$localImportance)))
  globalimp <- matrix((abs(rf$importance[, 1])/sum(abs(rf$importance[, 1]))), dim(x)[1], ncol(x), byrow = TRUE)
  return(list(localimp = localimp, globalimp = globalimp))
}

#' Create cluster with coefficients and fit regressions for each cluster

regcluster <- function(k, xcoef, x, y, dmax = ncol(x)){ #xcoef without intercept
  coefs <- matrix(0, k, ncol(x) + 1)
  adjr2 <- matrix(NA, k, 1)
  d <- ncol(x)
  if (k%%1!=0 || k==0) {
    stop("k is not an integer positive value.")
  }
  else if (k > 1) {
    cluster <- unlist(as.list(kernlab::specc(xcoef, centers = k)))
    for (i in 1:k) {
      if (d == 1) {
        models <- lm(y ~ ., data=data.frame(y = y[cluster==i], x[cluster==i, ]))
        coefs[i, ] <- coef(models)[-1]
        coefs[i, d+1] <- coef(models)[1]
        adjr2[i] <- summary(models)$adj.r.squared
      } else if (d > 1) {
        models <- leaps::regsubsets(y ~ ., 
                                    data=data.frame(y = y[cluster==i], 
                                                    x[cluster==i, which(apply(x[cluster==i, ], 2, sd) != 0)]), 
                                    method="forward", nvmax=dmax)
        est0 <- coef(models, which.min(summary(models)$bic))[-1]
        coefs[i, match(names(est0), colnames(data.frame(x)))] = est0
        coefs[i, d+1] <- coef(models, which.min(summary(models)$bic))[1]
        adjr2[i]=summary(models)$adjr2[which.min(summary(models)$bic)]
      }
    }
    
    adjr2w <- weighted.mean(as.vector(adjr2), as.vector(table(cluster)))
    
  } else if (k==1) {
    cluster <- rep(1, nrow(x))
    if (d == 1) {
      models <- lm(y ~ ., data=data.frame(y = y, x))
      coefs <- coef(models)[-1] 
      coefs[d+1] <- coef(models)[1]
      adjr2w <- summary(models)$adj.r.squared
    } else if (d > 1) {
      models <- leaps::regsubsets(y ~ ., data=data.frame(y, x[, which(apply(x, 2, sd) != 0)]),
                                  method="forward", nvmax=dmax)
      est0 <- coef(models, which.min(summary(models)$bic))[-1]
      coefs[match(names(est0), colnames(data.frame(x)))] = est0
      coefs[d+1] <- coef(models, which.min(summary(models)$bic))[1]
      adjr2w <- summary(models)$adjr2[which.min(summary(models)$bic)]
    }
    
  } 
  return(list(coefs = coefs, adjr2 = adjr2w, cluster = cluster ))
}


