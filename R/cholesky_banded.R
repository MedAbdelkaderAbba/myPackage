###
### Step 1
###
#' Check Matrix is k-banded
#'
#' \code{check_k_banded} returns a Boolean variable indicating whether or
#' not a matrix is k-banded.
#'
#' @param A The input matrix
#' @param k bandwidth parameter
#' @export
check_k_banded <- function(A, k) {
  b <- TRUE
  if( k<0 || !(k - floor(k)==0) ) stop("Argument k is not a nonnegative integer")
  n <- nrow(A)
  p <- ncol(A)
  m <- min(n, p)
  B <- matrix(T, n, p )
  for(i in 1:n) B[i, seq(max(i-k, 1) , min(i+k, p)) ] <- FALSE
  if(crossprod(A[B])!=0) b <-FALSE
  return(b)
}

###
### Step 2
###


#' Banded Cholesky
#'
#' \code{chol_banded} computes the Cholesky decomposition of a banded
#' positive definite matrix.
#'
#' @param A The input matrix
#' @param k bandwidth parameter
#' @param checks_off Boolean variable to turn on / off checks
#' @export
#'
chol_banded <- function(A, k, checks_off=FALSE) {

  if(!checks_off){
    if(!check_k_banded(A , k)) stop("The input matrix is not k banded")
    if(!isSymmetric(A)) stop("Input matrix need be symmetric")
    A <- Matrix::Matrix(Matrix::triu(A), sparse = T) # In case we pass the test for the 1st time We only need the upper triangular part

  }

  if(A[1,1]<=0) stop("Input matrix is not positive definite")
  A[1,1] <- sqrt(A[1,1])
  #if(length(A[-1,-1]) > 1){
    #print(A[1,1])
    A[1,2:min(k+1,ncol(A))] <- A[1,2:min(k+1,ncol(A))]/A[1,1]
    A[-1 , -1] <- A[-1 , -1] - Matrix::triu(Matrix::Matrix(Matrix::tcrossprod(A[1,2:ncol(A)])))
    if(length(A[-1,-1]) > 1){
      A[-1 , -1] <- chol_banded(A[-1 ,-1] , k , checks_off = TRUE)
    }else{
      A[-1,-1] <- sqrt(A[-1,-1])
    }


  return(A)

}




