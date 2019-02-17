###
### Part 2 - Question 1
###

#' Sweep k
#'
#' \code{sweep_k} applies the sweep operator to a symmetric matrix
#' on the kth diagonal entry if it is possible.
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index on which to sweep
#' @export
#' We will compute the operations on the upper half only.
sweep_k <- function(A, k) {
if(diag(A)[k] == 0) { stop("The k-th diagonal entry is equal to 0") }
if(!isSymmetric(A)) {stop("The input matrix is not symmetric." )}
n <- ncol(A)
s <- diag(A)[k]
B <- A
for(j in 1:n){
  for(i in 1:j){

    ifelse(j!=k & i!=k , B[i , j] <- B[i , j] - (A[ i ,  k ] * A[ j ,k ] / A[k , k]) ,
                                     B[i , j] <- A[i , j] /A[k , k] )

  }

}

B[k,k] = -B[k,k]/s
B[lower.tri(B)] <- t(B)[lower.tri(B)]

return(B)
}

###
### Part 2 - Step 2
###
#' Inverse Sweep k
#'
#' \code{isweep_k} applies the inverse sweep operator to a symmetric matrix
#' on the kth diagonal entry if it is possible.
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index on which to sweep
#' @export
isweep_k <- function(A, k) {
  if(diag(A)[k] == 0) { stop("The k-th diagonal entry is equal to 0") }
  if(!isSymmetric(A)) {stop("The input matrix is not symmetric." )}
  n <- ncol(A)
  s<- diag(A)[k]
  B <- A
  for(j in 1:n){
    for(i in 1:j){
      ifelse(j!=k & i!=k , B[i , j] <- B[i , j] - A[i ,k  ] *
               A[  j , k  ] / A[k , k] ,
             B[i , j] <- - A[i , j] /A[k , k] )
    }

  }

  B[k,k] = B[k,k]/s
  B[lower.tri(B)] <- t(B)[lower.tri(B)]

  return(B)


}



###
### Part 2 - step 3
###

#' Sweep
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index entry set on which to sweep
#' @export
sweep <- function(A, k=NULL) {
  if(is.null(k)) {
    for (i in 1 : ncol(A)  ) { A <- sweep_k(A , i ) }
  }else{
      for(i in 1:length(k)) {A <- sweep_k(A , k[i]) }
  }
  return(A)
}

###
### Part 2 - step 4
###

#' Sweep
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index entry set on which to sweep
#' @export
isweep <- function(A, k=NULL) {
  if(is.null(k)) {
    for (i in 1 : ncol(A)  ) { A <- isweep_k(A , i ) }
  }else{
    for(i in 1:length(k)) {A <- isweep_k(A , k[i]) }
  }
  return(A)
}


