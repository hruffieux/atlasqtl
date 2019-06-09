source("main.R")

context("Checking whether example has converged.")

test_that("The algorithm has converged", {
  expect_equal(vb$converged, TRUE)
})
