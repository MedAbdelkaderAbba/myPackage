#' Gradient Step
#'
#' @param gradf handle to function that returns gradient of objective function
#' @param x current parameter estimate
#' @param t step-size
#' @export
gradient_step <- function(gradf, x, t) {
  return(x-t*gradf(x))
}




#' Gradient Descent (Fixed Step-Size)
#'
#' @param fx handle to function that returns objective function values
#' @param gradf handle to function that returns gradient of objective function
#' @param x0 initial parameter estimate
#' @param t step-size
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent_fixed <- function(fx, gradf,  x0, t, max_iter=1e2, tol=1e-5) {
  hat_x <- x0
  n = max_iter
  obj_values <- rep(0 , n)
  grad_values <- rep(0 , n)
  rel_change <- rep(0,n)
  frel_change <- rep(0,n)

  temp = tol +  1
  i = 1
  b = TRUE
  while(b){
    hatx <- gradient_step(gradf  ,x0 ,t)
    obj_values[i] <- fx(hatx)
    grad_values[i] <- sqrt(crossprod(gradf(hatx)))
    rel_change[i] <- sqrt(crossprod(hatx-x0)/crossprod(x0))
    frel_change[i] <- abs(fx(hatx) - fx(x0))/abs(fx(x0))
    temp <- min(frel_change[i], frel_change[i]*abs(fx(x0)))
    i<- i+1
    x0 <- hatx
    b <- (temp>=tol)&(i<=n)
    #browser()
  }

  m <- cbind(obj_values ,grad_values , rel_change , frel_change)
  colnames(m) <- c("fx_values","gradf","x_rel","fx_rel")
  return(list(x_star = hatx , m = m[1:(i-1),]))

}
#
# Your function should return
#
# - The final iterate value
# - The objective function values
# - The 2-norm of the gradient values
# - The relative change in the function values
# - The relative change in the iterate values



#' Backtracking
#'
#' @param fx handle to function that returns objective function values
#' @param x current parameter estimate
#' @param t current step-size
#' @param df the value of the gradient of objective function evaluated at the current x
#' @param alpha the backtracking parameter
#' @param beta the decrementing multiplier
#' @export
backtrack <- function(fx, x, t, df, alpha=0.5, beta=0.9) {
  cond <- ( fx(x - t*df) >= fx(x) - alpha*t*crossprod(df) )
  #if(cond==NaN)
  #browser()
  while(drop(cond)){
    t <- beta*t
    cond <- ( fx(x - t*df) >= fx(x) - alpha*t*crossprod(df) )
   #browser()
  }
  return(t)
}

#Your function should return the selected step-size.


#' Gradient Descent (Backtracking Step-Size)
#'
#' @param fx handle to function that returns objective function values
#' @param gradf handle to function that returns gradient of objective function
#' @param x0 initial parameter estimate
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent_backtrack <- function(fx, gradf, x0, max_iter=1e2, tol=1e-3) {
  hat_x <- x0
  n = max_iter
  obj_values <- rep(0 , n)
  grad_values <- rep(0 , n)
  rel_change <- rep(0,n)
  frel_change <- rep(0,n)
  stepsize = 1
  temp = tol +  1
  i = 1
  b = TRUE
  while(b){
    stepsize <- backtrack(fx , x0 , stepsize , df = gradf(x0) , alpha = 0.5 , beta = 0.9)
    hatx <- x0 - stepsize*gradf(x0)
    obj_values[i] <- fx(hatx)
    grad_values[i] <- sqrt(crossprod(gradf(hatx)))
    rel_change[i] <- sqrt(crossprod(hatx-x0)/crossprod(x0))
    frel_change[i] <- abs(fx(hatx) - fx(x0))/abs(fx(x0))
    temp <- min(frel_change[i], frel_change[i]*abs(fx(x0)))
    i<- i+1
    x0 <- hatx
    b <- (temp>=tol)&(i<=n)
  }

  m <- cbind(obj_values ,grad_values , rel_change , frel_change)
  colnames(m) <- c("fx_values","gradf","x_rel","fx_rel")
  return(list(x_star = hatx , m = m[1:(i-1),]))



}
# Your function should return
#
# - The final iterate value
# - The objective function values
# - The 2-norm of the gradient values
# - The relative change in the function values
# - The relative change in the iterate values



#' Gradient Descent
#'
#' @param fx handle to function that returns objective function values
#' @param gradf handle to function that returns gradient of objective function
#' @param x0 initial parameter estimate
#' @param t step-size
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent <- function(fx, gradf, x0, t=NULL, max_iter=1e2, tol=1e-3) {
  ifelse(is.null(t) , return(gradient_descent_backtrack(fx , gradf , x0 , max_iter , tol)) ,
         return(gradient_descent_fixed(fx , gradf , x0 , t  ,max_iter, tol)))
}
#
# Your function should return
#
# - The final iterate value
# - The objective function values
# - The 2-norm of the gradient values
# - The relative change in the function values
# - The relative change in the iterate values


#' Compute kth order differencing matrix
#'
#' @param k order of the differencing matrix
#' @param n Number of time points
#' @export
myGetDkn <- function(k, n) {
  if(k==1){
    return( Matrix::bandSparse(n-1,n , k = c(0,1) ,diagonals = cbind(rep(-1,n-1) ,rep( 1 , n-1 )) ) )
  }else{
    return(  Matrix::crossprod(Matrix::t(myGetDkn(1 , n-k+1)) , myGetDkn(k-1 , n)  ) )
  }
}

#' Objective Function for HP-filtering
#'
#' @param y response
#' @param theta regression coefficient vector
#' @param Dkn sparse differencing matrix
#' @param lambda regularization parameter
#' @export
fx_hp <- function(y, theta, Dkn, lambda=0) {
  return(0.5*drop( crossprod(y-theta) + lambda*crossprod(crossprod(t(Dkn),theta)) ))
}

#' Gradient for HP-filtering
#'
#' @param y response
#' @param theta regression coefficient vector
#' @param Dkn sparse differencing matrix
#' @param lambda regularization parameter
#' @export
gradf_hp <- function(y, theta, Dkn, lambda=0) {
  return( - y + theta + lambda* crossprod(Dkn , crossprod(t(Dkn) , theta) ))
}



