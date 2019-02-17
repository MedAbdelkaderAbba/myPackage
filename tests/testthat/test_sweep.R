###
###
###
n <- 10
u <- matrix(rnorm(n), ncol=1)
A <- tcrossprod(u)
diag(A) <- diag(A) + 1
k<-1

#test_sweep<- function(A , k){
  test_that("sweep_k and isweep_k cancel each other",{
    expect_equal(A , isweep_k(sweep_k(A ,k) , k) )
  })
  test_that("sweep and isweep cancel each other",{
    expect_equal(A , isweep(sweep(A ,k) , k) )
  })
  test_that("Sweep operator functions correctly",{
    expect_equal(max(abs(crossprod(A , sweep(A)) + diag(ncol(A)))) , 0 )
  })
#}

