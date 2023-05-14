#' The proximal operator of the L2 norm
#'
#' \code{proxL2} computes the proximal operator of the L2 norm, tau*\| x \|_2
#' @param x the vector
#' @param tau the parameter
#' @useDynLib GMC
#' @export

proxL2 = function(x, tau){
  storage.mode(x) <- "double"
  storage.mode(tau) <- "double"

  x <- as.double(x)
  tau <- as.double(tau)
  p <- as.integer(length(x))
  x_prox <- as.double(rep(0, p))


  sol=.C("prox_L2", x=x, tau=tau, length=p,  x_prox=x_prox,  PACKAGE = 'GMC')

  return(sol$x_prox)
}



#' The matrix-vector multiplication using C code
#'
#' \code{mvprod_C} computes the computes the product of b=Ax
#' @param A the matrix multiplier
#' @param x the vector multiplier
#' @useDynLib GMC
#' @export

mvprod_C = function(A, x){

  a <- as.vector(A)
  x <- as.vector(x)
  m <- as.integer(nrow(A))
  n <- as.integer(ncol(A))

  storage.mode(a) <- "double"
  storage.mode(x) <- "double"

  b <- as.double(rep(0, m))


  sol=.C("test_mvprod", a=a, x=x, m=m, n=n, prod = b,  PACKAGE = 'GMC')

  return(sol$prod)
}




#' The matrix-vector multiplication using C code
#'
#' \code{mtvprod_C} computes the computes the product of b=A'x
#' @param A the matrix multiplier
#' @param x the vector multiplier
#' @useDynLib GMC
#' @export

mtvprod_C = function(A, x){

  a <- as.vector(A)
  x <- as.vector(x)
  m <- as.integer(nrow(A))
  n <- as.integer(ncol(A))

  storage.mode(a) <- "double"
  storage.mode(x) <- "double"

  b <- as.double(rep(0, n))


  sol=.C("test_mtvprod", a=a, x=x, m=m, n=n, prod = b, PACKAGE = 'GMC')
 #        PACKAGE = "GMC")

  return(sol$prod)
}




#' The matrix-matrix multiplication using C code
#'
#' \code{mmprod_C} computes the computes the product of two matrices C=AB
#' @param A one matrix multiplier
#' @param B tne other matrix multiplier
#' @useDynLib GMC
#' @export

mmprod_C = function(A, B){

  a <- as.vector(A)
  b <- as.vector(B)
  m <- as.integer(nrow(A))
  k <- as.integer(ncol(A))
  n <- as.integer(ncol(B))

  storage.mode(a) <- "double"
  storage.mode(b) <- "double"

  c <- as.double(rep(0, m*n))


  sol=.C("test_mmprod", a=a, b=b, c=c, m=m, k=k, n=n,  PACKAGE = 'GMC')

  C <- matrix(sol$c, nrow = m)
  return(C)
}


#' The matrix-matrix multiplication using C code
#'
#' \code{mtmprod_C} computes the computes the product of two matrices C=A'B
#' @param A one matrix multiplier
#' @param B tne other matrix multiplier
#' @useDynLib GMC
#' @export

mtmprod_C = function(A, B){

  a <- as.vector(A)
  b <- as.vector(B)
  k <- as.integer(nrow(A))
  m <- as.integer(ncol(A))
  n <- as.integer(ncol(B))

  storage.mode(a) <- "double"
  storage.mode(b) <- "double"

  c <- as.double(rep(0, m*n))


  sol=.C("test_mtmprod", a=a, b=b, c=c, m=m, k=k, n=n,  PACKAGE = 'GMC')

  C <- matrix(sol$c, nrow = m)
  return(C)
}
