#' Ridge Regression
#'
#' \code{ridge_regression} returns the ridge regression coefficient estimates
#' for a sequence of regularization parameter values.
#'
#' @param y response variables
#' @param X design matrix
#' @param lambda vector of tuning parameters
#' @export
ridge_regression <- function(y, X, lambda) {
  n = nrow(X) ; p = ncol(X)
  if(length(y)!=n) stop("Response vector and design matrix are not comfortable")
  if(!all(lambda>=0)) stop("Tuning parameters should be nonnegative")
  X_svd <- svd(X) ; U <- X_svd$u ; d <- X_svd$d ; V <- X_svd$v
  betahat <- matrix(0 , ncol = length(lambda) , nrow = p)
  for(i in 1:length(lambda)){
    t <- (d/(d^2+lambda[i]))*drop(crossprod(U,matrix(y,ncol = 1)))
    betahat[,i] <-  crossprod(t(V) , t)
  }
  return(betahat)
}


#' Leave One Out
#'
#' \code{leave_one_out} returns the leave-one-out
#' for a sequence of regularization parameter values.
#'
#' @param y response variables
#' @param X design matrix
#' @param lambda vector of tuning parameters
#' @export
leave_one_out <- function(y, X, lambda) {
  hat_y <- crossprod(t(X) , ridge_regression(y , X, lambda))
  s <- svd(X); U <- s$u; d <- s$d
  Loo <- lambda
  for(i in 1:length(lambda)){
    h <- diag(U%*%diag( d^2/(d^2 + lambda[i]) )%*%t(U))
    Loo[i] <- mean(((y - hat_y[,i])/(1-h))^2)

    #H <- drop( rowSums((U^2)[,1:min(n,p)]%*%diag(d^2/(lambda[i] + d^2))))
    #Loo[i] <- mean( ((y-drop(hat_y[,i]))/(1-H))^2 )
  }
  return(Loo)
}





