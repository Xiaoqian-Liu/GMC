#' The forward-backward (FB) algorithm for solving the group GNC problem
#'
#' \code{FB} computes the solution of the group GMC problem using the forward-backward algorithm.
#' @param y The response vector.
#' @param X The design matrix
#' @param beta0 The initial value of beta.
#' @param v0 The initial value of v.
#' @param lambda The regularization parameter.
#' @param alpha The convexity-preserving parameter.
#' @param group A vector with the same length as beta0. Each component indicates which group
#'              the corresponding covariate variable belongs to.
#' @param K The weight vector for the group penalty.
#' @param maxiters The maximal number of iterations.
#' @param tol A tolerance for early stopping.
#' @useDynLib GMC
#' @export

FB <- function(y, X, beta0, v0, lambda, alpha, group, K, maxiters=1e3, tol=1e-4){

  n <- as.integer(nrow(X))
  p <- as.integer(ncol(X))

  X_data <- as.vector(X)

  XX <- mtmprod_C(X, X)
  XX_data <- as.vector(XX)

  Z <- alpha*XX
  Z_data <- as.vector(Z)
  if(alpha!=1){
      t <- 1.9*(1-alpha)/((1-2*alpha+2*alpha^2)*norm(XX, "2"))
  }else{
    t <- 1.9/norm(XX, "2")
  }


  storage.mode(t) = "double"
  storage.mode(X_data) = "double"
  storage.mode(XX_data) = "double"
  storage.mode(Z_data) = "double"
  storage.mode(y) = "double"
  storage.mode(beta0) = "double"
  storage.mode(v0) = "double"

  storage.mode(lambda) = "double"
  storage.mode(group) <- "integer"
  storage.mode(K) <- "double"
  storage.mode(maxiters) <- "integer"
  storage.mode(tol) <- "double"

  beta <- rep(0, p)
  v <- rep(0, p)
  iters <- 0
  storage.mode(beta) = "double"
  storage.mode(v) = "double"
  storage.mode(iters) <- "integer"

  sol=.C("FB_call", y=y, n=n, p=p, X_data=X_data, Z_data=Z_data, XX_data=XX_data,
         beta0=beta0, v0=v0, lambda=lambda, group=group, K=K, t=t, max_iters=maxiters,
         tol=tol, beta=beta, v=v, iters=iters, PACKAGE = 'GMC')

  return(list(beta=sol$beta, v=sol$v, iters=sol$iters))
}

