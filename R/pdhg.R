#' Adaptive PDHG algorithm for group GMC
#'
#' @param beta The initial value for variable beta. If missing, set it as a zero vector.
#' @param v The initial value for variable v. If missing, set it as a zero vector.
#' @param y The response vector.
#' @param X The design matrix in group GMC.
#' @param Z The Z matrix in group GMC (Z = lambda*B'B/n).
#' @param group A vector with the same length as beta. Each component indicates which group
#'              he corresponding covariate variable belongs to.
#' @param K The weight vector for the group penalty or called group multiplier.
#' @param lambda The regularization parameter in group GMC.
#' @param maxiters The max number of iterations.
#' @param tol A tolerance for early stopping.
#' @param adaptive If 'TRUE', use adaptive method.
#'                 If 'FALSE', turn off the adaptivity and users must set tau and
#'                  sigma as the setpsize for the classical PDHG.
#' @param backtrack If 'TRUE', use backtracking method.
#' @param L The reciprical spectral radius of Z'Z. If missing, the algorithm will give an approcimate value.
#' @param tau  The intial stepsize for the primal variables.
#'             If missing, the algorithm will initialize it as the squared root of L.
#' @param sigma The intial stepsize for the dual variables.
#'              If missing, the algorithm will initialize it as L/tau.
#' @param a  Intial adaptive update strength for stepsizes.
#' @param eta  How fast does the adaptivity level decay.
#' @param Delta Update stepsizes when primal/dual ratio exceeds Delta.
#' @param gamma Used to determine when need to backtrack to maintain positivity conditions.
#' @param b Adaptivity parameter used for backtracking update.
#' @param stopRule Three choices in total. Default is ratioResidual.
#'            "ratioResidual": primal/maxPrimal<tol & dual/maxDual < tol.
#'            "residual": primal < tol & dual < tol.
#'            "iteration" iter > maxiters.
#' @param ShowTime Whether to show the time cost for computing the solution.
#' @useDynLib GMC
#' @export

adaptive_pdhg <- function(beta, v, y, X, Z, group, K, lambda, maxiters = 5e2, tol = 1e-3,
              adaptive = TRUE, backtrack = TRUE, L, tau, sigma,
              a = 0.5, eta = 0.95, Delta = 2, gamma = 0.75, b=0.95,
              stopRule = "ratioResidual", ShowTime=FALSE){

  # Approximate the reciprical spectral radius of A'A
  if(missing(L)){
    set.seed(1234)
    xr = rnorm(length(beta))
    Zx = mvprod_C(Z, xr) #Z%*%xr
    ZtZx = mtvprod_C(Z, Zx) #t(Z)%*%Zx
    L = 2*norm(Zx, "2")/norm(ZtZx, "2")
  }
  if(missing(tau)){
    tau = sqrt(L)
  }
  if(missing(sigma)){
    sigma = L/tau
  }

  if(missing(beta)){
    if(missing(v)){
      beta <- v <- rep(0, ncol(X))
    }else{
      beta <- v
    }
  }

  storage.mode(beta) <- "double"
  storage.mode(v) <- "double"

  y <- as.vector(y)
  storage.mode(y) <- "double"

  A_data <- as.vector(X)
  Z_data <- as.vector(Z)
  storage.mode(A_data) <- "double"
  storage.mode(Z_data) <- "double"

  n <- as.integer(nrow(X))
  p <-  as.integer(ncol(X))

  storage.mode(lambda) <- "double"
  group <- as.integer(group)
  storage.mode(K) <- "double"



  maxiters <- as.integer(maxiters)
  tol <- as.double(tol)
  adaptive <- as.integer(adaptive)
  backtrack <- as.integer(backtrack)

  tau <- as.double(tau)
  sigma <- as.double(sigma)
  L <- as.double(L)
  a <- as.double(a)
  eta <- as.double(eta)
  Delta <- as.double(Delta)
  gamma <- as.double(gamma)
  b <- as.double(b)

  if(stopRule == "ratioResidual"){
    Stoprule <- 1
  }
  if(stopRule == "residual"){
    Stoprule <- 2
  }
  if(stopRule == "iteration"){
    Stoprule <- 3
  }
  Stoprule <- as.integer(Stoprule)


  estimate_of_x <- double(p)
  estimate_of_v <- double(p)
  Num_of_iters <- as.integer(0)
  Presidual <- double(maxiters)
  Dresidual <- double(maxiters)
  tau_seq <- double(maxiters)
  sigma_seq <- double(maxiters)
  updates <- as.integer(0)


  time <- proc.time()
  sol=.C("test_PDHG", x = beta, v = v, y = y, A_data = A_data, Z_data = Z_data, n = n, p = p,
         group = group, K = K, lambda = lambda, maxiters = maxiters, tol = tol, adaptive = adaptive,
         backtrack = backtrack, tau = tau, sigma = sigma, L =L, a = a, eta = eta, Delta = Delta,
         gamma = gamma, b = b, stopRule = Stoprule, estimate_of_x = estimate_of_x,
         estimate_of_v = estimate_of_v, Num_of_iters = Num_of_iters, Presidual = Presidual,
         Dresidual = Dresidual, tau_seq = tau_seq, sigma_seq = sigma_seq, updates = updates,
         PACKAGE = 'GMC')

  if(ShowTime){
    print(proc.time()-time)
  }

  iters <- sol$Num_of_iters

  outs <- list()
  outs$x <- sol$estimate_of_x
  outs$v <- sol$estimate_of_v
  outs$tau <- sol$tau_seq[1:iters]
  outs$sigma <- sol$sigma_seq[1:iters]
  outs$p <- sol$Presidual[1:iters]
  outs$d <- sol$Dresidual[1:iters]
  outs$iters <- sol$Num_of_iters
  outs$updates <- sol$updates

  return(outs)

}



