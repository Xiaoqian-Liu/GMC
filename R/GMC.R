#' Linear regression with the (group) GMC penalty
#'
#' \code{GMC} fit regularization paths for linear regression with the (group) GMC penalty
#' over a grid of values for the regularization parameter lambda.
#'
#' @param y The response vector.
#' @param X The design matrix without an intercept. grGMC standardizes the data and includes an intercept by default.
#' @param alpha The convexity-preserving parameter alpha in the group GMC penalty taking on a value between 0 and 1. Default is 0.8.
#' @param group A vector of length ncol(X). Each component indicates which group the corresponding covariate variable belongs to.
#' @param group.multiplier The weight vector for the group GMC penalty.
#' @param nlambda The number of values for the regularization parameter. Default is 50. The function will automatically computes
#' a grid of 'nlambda' lambda values with the first one producing a zero solution.
#' @param lambdaSeq A user-specified sequence of decrasing lambda values. It is typically left unspecified and let the function automatically
#' computes a grid of 'nlambda' lambda values.
#' @param lambda.min The smallest lambda value, as a fraction of the maximum lambda.
#'                   Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param log.lambda Whether to compute the grid values of lambda on log scale (default) or linear scale. Default is TRUE.
#' @param algorithm Which algorithm to use for solving the group GMC problem. PDHG (default) or FB.
#'                  We recommend using PDHG for large value and FB for small alpha.
#' @param maxiters The maximal number of iterations. Typically, PDHG only needs 100 but FB needs 1e4.
#' @param tol A tolerance for early stopping. Default is 1e-6.
#' @param returnX Whether to return the standardized data.
#' @param ShowTime Whether to show the time cost for computing the solution path.
#' @return \code{beta} The fitted matrix of coefficients.
#' The number of rows equals to the number of coefficients plus one (intercept),
#' and the number of columns equals to the length of the lambda sequence.
#' @return \code{family} Only "gaussian" is allowed for now.
#' @return \code{penalty} We call our penalty as GMC.
#' @return \code{group} Same as the argument 'group'.
#' @return \code{lambdaSeq} The sequence of lambda values for the solution path.
#' @return \code{alpha} Same as the argument 'alpha'.
#' @return \code{iters} A vector containing the number of iterations
#' for the algorithm to compute each solution.
#' @return \code{bb} The fitted matrix of coefficients for the standardized data.
#'  This return is used as a warm start when having a sequence of lambda.
#' @return \code{group.multiplier} A named vector containing the multiplicative constant
#'  applied to each group's penalty.
#' @author Xiaoqian Liu
#' @export
#' @examples
#' set.seed(1234)
#' n <- 50
#' p <- 10
#' # true coefficients
#' beta <- rep(0, p)
#' beta[1:4] <- c(1,1,-1,-1)
#' # group information
#' group <- rep(1:5, each=p/5)
#' # set the group weights
#' gp <- unique(group)
#' J <- length(gp)
#' K <- rep(0,J)
#' for (j in 1:J) {
#'   K[j] <- sqrt(length(which(group==gp[j])))
#' }
#'
#' # design matrix
#' X <- matrix(rnorm(n*p), nrow = n)
#' # response
#' y <- as.vector(X%*%beta+rnorm(n, mean=0, sd = 0.5))
#'
#' # fit with grGMC function
#' fit_GMC <- GMC(y=y, X=X, alpha=0.9, group=group, group.multiplier=K, nlambda = 50, ShowTime=FALSE)

GMC <- function( y, X, alpha=0.8,  group=1:ncol(X),  group.multiplier, nlambda=50,
                  lambdaSeq, lambda.min={if (nrow(X) > ncol(X)) 1e-3 else .05}, log.lambda=TRUE,
                  algorithm = 'PDHG', maxiters = 1e2, tol = 1e-6, returnX=FALSE, ShowTime=TRUE){

  time <- proc.time()

  #construct new X and y as "grpreg"
  newy <- grpreg:::newY(y, family="gaussian")
  XG <- grpreg:::newXG(X, g=group, m=group.multiplier, attr(newy, 'm'), bilevel=FALSE)
  newX <- XG$X
  group1 <- XG$g      #new group vector
  K1 <- XG$m          #new group multiplier

  if (nrow(newX) != length(newy))
    stop("X and y do not have the same number of observations", call.=FALSE)



  # Fit the model
  n <- length(newy)
  p <- ncol(newX)



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


 # get Z for PDHG
  if(algorithm == 'PDHG'){
    Z <- alpha*t(newX)%*%newX/n
  }

  BetaMat <- matrix(0, nrow = p, ncol = nlambda)
  iters <- rep(0, nlambda)

  beta0 <- v0 <- BetaMat[, 1]
  for (i in 1:nlambda) {
    # set lambda
    lambda <- lambdaSeq[i]

    # run the algorithm
    if(algorithm == 'PDHG'){
      out_PDHG <- adaptive_pdhg(beta=beta0, v=v0, y=newy, X=newX, Z=Z, group=group1, K=K1, lambda = lambda,
                                 maxiters=maxiters, tol=tol)
      BetaMat[, i] <- beta0  <- as.vector(out_PDHG$x)
      v0 <- as.vector(out_PDHG$v)
      iters[i] <- out_PDHG$iters
    }else{
      out_FB <- FB(y=newy, X=newX, beta0=beta0, v0=v0, lambda=lambda, alpha=alpha, group=group1, K=K1, maxiters=maxiters, tol=tol)
      BetaMat[, i] <- beta0 <- out_FB$beta
      v0 <- out_FB$v
      iters[i] <- out_FB$iters
    }

  }

  b <- rbind(mean(y), BetaMat)   #add the intercept

  # Unstandardize
  b <- grpreg:::unorthogonalize(b, XG$X, XG$g)
  if (XG$reorder)  b[-1,] <- b[1+XG$ord.inv, ]
  beta <- grpreg:::unstandardize(b, XG)

  # Names
  varnames <- c("(Intercept)", XG$names)
  dimnames(beta) <- list(varnames, round(lambdaSeq, digits=4))

  val <- structure(list(beta = round(beta, 8),
                        family="gaussian",
                        penalty = "GMC",
                        group = factor(group),
                        lambdaSeq = lambdaSeq,
                        alpha=alpha,
                        iters=iters,
                        bb=BetaMat,
                        group.multiplier = XG$m),
                        class="GMC")

  if (returnX) {
    val$XG <- XG
    val$y <- newy
    if(algorithm == 'PDHG'){
    val$Z <- Z
    }
  }

  if(ShowTime){
    print(proc.time()-time)
  }

  return(val)

}
