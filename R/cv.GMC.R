#' Cross validation for (group) GMC
#'
#' \code{cv.GMC} conduct cross validation for (group) GMC penalized regression over a grid of values for the regularization parameter lambda.
#' @param y The response vector.
#' @param X The design matrix.
#' @param alpha The convexity-preserving parameter taking on a value between 0 and 1. Default is 0.8.
#' @param group A vector of length ncol(X). Each component indicates which group the corresponding covariate variable belongs to.
#' @param group.multiplier The weight vector for the group GMC penalty.
#' @param nlambda The number of values for the regularization parameter. Default is 50. The function will automatically computes
#' a grid of 'nlambda' lambda values with the first one producing a zero solution.
#' @param lambdaSeq A user-specified sequence of decrasing lambda values. It is typically left unspecified and let the function automatically
#' computes a grid of 'nlambda' lambda values.
#' @param lambda.min The smallest lambda value, as a fraction of the maximum lambda.
#'                   Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param log.lambda Whether to compute the grid values of lambda on log scale (default) or linear scale. Default is TRUE.
#' @param nfolds The number of cross validation folds. Default is 5.
#' @param algorithm Which algorithm to use for solving the group GMC problem. PDHG (default) or FB.
#'                  We recommend using PDHG for large value and FB for small alpha.
#' @param maxiters The maximal number of iterations. Typically, PDHG only needs 100 but FB needs 1e4.
#' @param tol A tolerance for early stopping. Default is 1e-6.
#' @param seed Random seed for reproduccible results.
#' @param trace Whether to track the progress of the cross validation. Default is FALSE.
#' @return \code{cve} The cross-validation error for each value of lambda, averaged over nfolds.
#' @return \code{cvse} The standard error associated with each value of cve.
#' @return \code{lambdaSeq} The sequence of lambda values along which the cross-validation error was calculated.
#' @return \code{lambda_min} The value of lambda with the minimum cross-validation error.
#' @return \code{min} The index of lambda corresponding to lambda_min.
#' @return \code{lambda_1se} The value of lambda selected by the 1se rule.
#' @return \code{min_1se} The index of lambda corresponding to lambda_1se.
#' @author Xiaoqian Liu
#' @export
cv.GMC <- function(y, X, alpha=0.8, group=1:ncol(X), group.multiplier, nlambda=50, lambdaSeq,
                     lambda.min={if (nrow(X) > ncol(X)) 1e-3 else .05}, log.lambda=TRUE, nfolds=5,
                     algorithm='PDHG', maxiters=1e2, tol=1e-6, seed=1234, trace=FALSE) {


  # Get standardized X, y
  newy <- grpreg:::newY(y, family="gaussian")
  XG <- grpreg:::newXG(X, g=group, m=group.multiplier, attr(newy, 'm'), bilevel=FALSE)
  newX <- XG$X
  group1 <- XG$g      #new group vector
  K1 <- XG$m          #new group multiplier
  n <- length(y)


  # get the lambda sequence
  if(missing(lambdaSeq)){
    # compute the maximum lambda value
    gp <- unique(group1)
    J <- length(gp)
    lambda.max <- 0
    for (j in 1:J) {
      ind <- which(group1==gp[j])
      lambda.max <- max(lambda.max, norm(t(newX[, ind])%*%newy, "F")/(n*K1[j]))
    }
    # generate the lambda sequence
    if(log.lambda){
      lambdaSeq <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))
    }else{
      lambdaSeq <- seq(lambda.max, lambda.min*lambda.max, length=nlambda)
    }
  }else{
    nlambda <- length(lambdaSeq)
  }

  # Set up folds
  if (!missing(seed)) set.seed(seed)
  fold <- sample(1:n %% nfolds)
  fold[fold==0] <- nfolds



  # Do cross-validation
  E <- Y <- matrix(NA, nrow=length(y), ncol=nlambda)

  cv.args <- list()
  cv.args$alpha <- alpha
  cv.args$group <-  group1
  cv.args$group.multiplier <- K1
  cv.args$nlambda <-  nlambda
  cv.args$lambdaSeq <-  lambdaSeq
  cv.args$lambda.min <-  lambda.min
  cv.args$log.lambda <- log.lambda
  cv.args$algorithm <- algorithm
  cv.args$ShowTime <- FALSE
  cv.args$maxiters <- maxiters
  cv.args$tol <- tol

  for (i in 1:nfolds) {
    if (trace) cat("Starting CV fold #", i, sep="","\n")
    res <- cv_fold(i, newy, newX, fold, cv.args)
    Y[fold==i, ] <- res$yhat
    E[fold==i, ] <- res$loss
  }

  # Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[, ind, drop=FALSE]
  Y <- Y[, ind]
  lambda <- lambdaSeq[ind]

  # Return
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(n)
  min <- which.min(cve)


  #find the lambda.1se
  for (k in min:1) {
    if(cve[k]>cve[min]+cvse[min])
      break
  }
  lambda_1se <- lambda[k+1]
  min_1se <- k+1

  val <- list(cve=cve, cvse=cvse, lambdaSeq=lambda, min=min,
              min_1se=min_1se, lambda_min=lambda[min], lambda_1se=lambda_1se)

}



#' Cross-validation fold
#'
#' \code{cv_fold} fit the model for the i-th fold
#' @param i The i-th fold
#' @param y The response vector
#' @param X The design matrix
#' @param fold The fold indices
#' @param cv.args The arguments
#' @export
cv_fold <- function(i, y, X, fold, cv.args) {
  cv.args$X <- X[fold!=i, , drop=FALSE]
  cv.args$y <- y[fold!=i]
  fit.i <- do.call("GMC", cv.args)

  X2 <- X[fold==i, , drop=FALSE]
  y2 <- y[fold==i]

  L <- cv.args$nlambda
  loss <- yhat <- matrix(0, nrow = length(y2), ncol = L)

  for (l in 1:L) {
    yhat[, l] <-  X2%*%fit.i$beta[-1, l]
    loss[,l] <- (y2-yhat[, l])^2
  }

  list(loss=loss, yhat=yhat)
}
