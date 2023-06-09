% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FB.R
\name{FB}
\alias{FB}
\title{The forward-backward (FB) algorithm for solving the group GNC problem}
\usage{
FB(y, X, beta0, v0, lambda, alpha, group, K, maxiters = 1000, tol = 1e-04)
}
\arguments{
\item{y}{The response vector.}

\item{X}{The design matrix}

\item{beta0}{The initial value of beta.}

\item{v0}{The initial value of v.}

\item{lambda}{The regularization parameter.}

\item{alpha}{The convexity-preserving parameter.}

\item{group}{A vector with the same length as beta0. Each component indicates which group
the corresponding covariate variable belongs to.}

\item{K}{The weight vector for the group penalty.}

\item{maxiters}{The maximal number of iterations.}

\item{tol}{A tolerance for early stopping.}
}
\description{
\code{FB} computes the solution of the group GMC problem using the forward-backward algorithm.
}
