#####
##### Derivatives File
#####


#' Directional Derivative Exact
#'
#' \code{dd_exact} computes the exact directional derivative of the multivariate
#' function f, at the point a, in the direction b.
#'
#' @param gradf handle to function that returns the gradient of the function f
#' @param a point at which the directional derivative is evaluated
#' @param b the direction vector
#' @export
dd_exact <- function(gradf, a, b) {
return(sum(diag(crossprod(gradf(a) , b))))
}








#' Directional Derivative Approximate
#'
#' \code{dg_approx} computes an approximate directional derivative of the multivariate
#' function f, at the point a, in the direction b.
#'
#' @param f handle to function that returns the function f
#' @param a point at which the directional derivative is evaluated
#' @param b the direction vector
#' @param h small displacement
#' @export
dg_approx <- function(f, a, b, h = 1e-13) {
return( (f(a+ h*b) - f(a))/h )
}
