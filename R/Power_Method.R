###
### Power Method file
###

#' Power Method for dense matrices
#'
#' \code{power_method_dense} applies the power method to estimate
#' the eigenvector of a dense matrix associated with its eigenvector of largest magnitude.
#'
#' @param A The input matrix
#' @param x initial guess
#' @param max_iter maximum number of iterations
#' @param tol relative error stopping criterion
#' @export
power_method_dense <- function(A, x, max_iter, tol) {
stopifnot(ncol(A)==length(x))
x_new <- drop(A%*%x)
x_new <- x_new/norm(x , type = '2')
i = 1
repeat{
  x<- x_new
  x_new <- drop(A%*%x)
  x_new <- x_new/norm(x_new , type = '2')
  i<- i+1
  if( norm(x_new + x, type = '2') < tol*norm(x , type = '2') | norm(x_new - x, type = '2') <
      tol*norm(x , type = '2') | i>= max_iter){ break }

}
print(i)
return(x_new)
}







#####
##### Step 2
#####

#' Power Method for sparse matrices
#'
#' \code{power_method_sparse} applies the power method to estimate
#' the eigenvector of a sparse matrix associated with its eigenvector of largest magnitude.
#'
#' @param A The input matrix
#' @param x initial guess
#' @param max_iter maximum number of iterations
#' @param tol relative error stopping criterion
#' @export
power_method_sparse <- function(A, x, max_iter, tol) {
  stopifnot(ncol(A)==length(x) )
  if(class(A)!= "dgCMatrix"){ stop("The input matrix is not sparse") }
  x_new <- drop(A%*%x)
  x_new <- x_new/norm(x , type = '2')
  i = 1
  repeat{
    x<- x_new
    x_new <- drop(A%*%x)
    x_new <- x_new/norm(x_new , type = '2')
    i<- i+1
    if( norm(x_new + x, type = '2') < tol*norm(x , type = '2') | norm(x_new - x, type = '2') <
        tol*norm(x , type = '2') | i>= max_iter){ break }

  }
  print(i)
  return(x_new)
}


#####
##### Step 3
#####

#' Power Method for low rank matrices
#'
#' \code{power_method_low_rank} applies the power method to estimate
#' the eigenvector of a low rank matrix associated with its eigenvector of largest magnitude.
#'
#' @param U The left input factor matrix
#' @param V The right input factor matrix
#' @param x initial guess
#' @param max_iter maximum number of iterations
#' @param tol relative error stopping criterion
#' @export
power_method_low_rank <- function(U, V, x, max_iter, tol) {

  if(ncol(U)!= ncol(V)){ stop("The factor matrices do not have the number of columns") }
  x_new <- drop(U%*%crossprod(V , x))
  x_new <- x_new/norm(x , type = '2')
  i = 1
  repeat{
    x<- x_new
    x_new <- drop(U%*%crossprod(V,x))
    x_new <- x_new/norm(x_new , type = '2')
    i<- i+1
    if( norm(x_new + x, type = '2') < tol*norm(x , type = '2') | norm(x_new - x, type = '2') <
        tol*norm(x , type = '2') | i>= max_iter){ break }

  }
  print(i)
  return(x_new)
}



#####
##### Step 4
#####

#' Power Method for sparse + low rank matrices
#'
#' \code{power_method_sparse_plus_low_rank} applies the power method to estimate
#' the eigenvector of a sparse + low rank matrix associated with its eigenvector of largest magnitude.
#'
#' @param S sparse input matrix term
#' @param U The left input factor matrix term
#' @param V The right input factor matrix term
#' @param x initial guess
#' @param max_iter maximum number of iterations
#' @param tol relative error stopping criterion
#' @export
power_method_sparse_plus_low_rank <- function(S, U, V, x, max_iter, tol) {
  stopifnot(ncol(S)==length(x) )
  if(class(S)!= "dgCMatrix"){ stop("The input matrix is not sparse") }
  if(ncol(U)!= ncol(V)){ stop("The factor matrices do not have the number of columns") }
  x_new <- drop(S%*%x + U%*%crossprod(V , x))
  x_new <- x_new/norm(x , type = '2')
  i = 1
  repeat{
    x<- x_new
    x_new <- drop(S%*%x + U%*%crossprod(V , x))
    x_new <- x_new/norm(x_new , type = '2')
    i<- i+1
    if( norm(x_new + x, type = '2') < tol*norm(x , type = '2') | norm(x_new - x, type = '2') <
        tol*norm(x , type = '2') | i>= max_iter){ break }

  }
  print(i)
  return(x_new)
}
