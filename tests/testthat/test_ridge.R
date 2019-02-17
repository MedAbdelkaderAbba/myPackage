u <- rnorm(10)
X <- matrix(rnorm(54) , ncol = 6 , nrow = 9)
l <- c(rpois(5,2) , -.5 , rgamma(5 , 2,2))
y <- rnorm(nrow(X))
test <- test_that("ridge_regression returns the appropriate error message",{
  expect_error( ridge_regression(u , X , lambda = rep(2,5)) , "Response vector and design matrix are not comfortable" )
  expect_error(ridge_regression(y , X , lambda = l) ,"Tuning parameters should be nonnegative" )
})

test <- test_that("ridge_regression returns the correct values of the coefficients",{

  expect_equal(crossprod( crossprod(X) + 2*diag(ncol(X)) , ridge_regression( y , X , lambda = 2) ) , crossprod(X , y))
})
