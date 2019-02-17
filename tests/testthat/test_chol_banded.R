####
#### Running few checks
####

k <- 3; m <- 10
A <- matrix(rnorm(m**2), m, m)
#A <- A + t(A)
A <- round(10*A)/10
for (i in 1:m) {
  for (j in 1:m) {
    if (abs(i - j) > k) A[i,j] <- 0
  }
}
A[lower.tri(A)] <- 0
A <- crossprod(A)
l <- sample(1 , 1:nrow(A))
P <- A - 2*A[l,l]*tcrossprod(diag(m)[,l])

test_that("check_k_banded returns the appropriate error message", {
  expect_error(check_k_banded(A , -1) , "Argument k is not a nonnegative integer")
  expect_error(check_k_banded(A , 3.5) , "Argument k is not a nonnegative integer")
})

test_that("check_k_banded correctly works", {
  expect_equal(check_k_banded(A , 1) , FALSE)
  expect_equal(check_k_banded(A , 3) , TRUE)
})

test_that("chol_banded returns appropriate error message" , {
  expect_error(chol_banded(as.matrix(Matrix::tril(A)) , 3 ) , "Input matrix need be symmetric" )
  expect_error(chol_banded(A , 2 ) , "The input matrix is not k banded" )
  expect_error(chol_banded(P, 3 ) ,
               "Input matrix is not positive definite" )
})

test_that("chol_banded works correctly" , {
expect_equal(crossprod(as.matrix(chol_banded(A , 3))) , A)
})
