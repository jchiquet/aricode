library(testthat)
library(lifecycle)
library(aricode)

test_that("sortPairs is deprecated", {
  c1 <- rep(1:3, each = 3)
  c2 <- rep(4:6, each = 3)
  lifecycle::expect_deprecated(sortPairs(c1, c2))
})

test_that("clustComp is deprecated", {
  c1 <- rep(1:3, each = 3)
  c2 <- rep(4:6, each = 3)
  lifecycle::expect_deprecated(clustComp(c1, c2))
})

test_that("MARIraw is deprecated", {
  c1 <- rep(1:3, each = 3)
  c2 <- rep(4:6, each = 3)
  lifecycle::expect_deprecated(MARIraw(c1, c2))
})
