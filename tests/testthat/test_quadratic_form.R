
nn <- 10
kVec <- runif(nn)
kAllZeroVec <- rep(0, nn)

# quadratic matrix of with appropriate dimension to kVec
kMat <- matrix(rnorm(nn^2), ncol = nn)
# this matrix is positive definite
kPosDefMat <- t(kMat) %*% kMat

eigen.results <- eigen(kPosDefMat)

test_that("Test_quadratic_form", {
  expect_equal(quadratic_form(kMat, kVec), 
               c(t(kVec) %*% kMat %*% kVec))
  # 0 times 0 is zero
  expect_equal(0, quadratic_form(matrix(0), 0))
  # first argument must be a matrix
  expect_error(quadratic_form(kVec, kVec))
  # for positive def matrix the quadratic form is positive
  expect_true(quadratic_form(kPosDefMat, kVec) > 0)
  # for eigenvalue it must hold that x' A x = lambda
  expect_equal(quadratic_form(kPosDefMat, eigen.results$vector[,1]), eigen.results$value[1])
})

