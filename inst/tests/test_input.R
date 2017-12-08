library(testthat)
library(aricode)

test_that("Testing supposedly handled input types", {

  c1 <- rep(1:3,each=3) * 1.5
  c2 <- rep(4:6,each=3) * 2.7

  cat("\n-numeric type-")
  c1.numeric <- c1
  c2.numeric <- c2
  ref.object <- list(spMat = NULL, levels = list(c1 = unique(c1.numeric), c2 = unique(c2.numeric)), nij = rep(3,3), ni. =  rep(3,3),  n.j = rep(3,3))
  expect_that(sortPairs(c1.numeric,c2.numeric), equals(ref.object))

  cat("\n-integer type-")
  c1.integer <- as.integer(c1)
  c2.integer <- as.integer(c2)
  ref.object <- list(spMat = NULL, levels = list(c1 = unique(c1.integer), c2 = unique(c2.integer)), nij = rep(3,3), ni. =  rep(3,3),  n.j = rep(3,3))
  expect_that(sortPairs(c1.integer,c2.integer), equals(ref.object))

  cat("\n-character type-")
  c1.char <- as.character(c1)
  c2.char <- as.character(c2)
  ref.object <- list(spMat = NULL, levels = list(c1 = unique(c1.char), c2 = unique(c2.char)), nij = rep(3,3), ni. =  rep(3,3),  n.j = rep(3,3))
  expect_that(sortPairs(c1.char,c2.char), equals(ref.object))

  cat("\n-factor type-")
  c1.factor <- as.factor(c1.char)
  c2.factor <- as.factor(c2.char)
  ref.object <- list(spMat = NULL, levels = list(c1 = levels(c1.factor), c2 = levels(c2.factor)), nij = rep(3,3), ni. =  rep(3,3),  n.j = rep(3,3))
  expect_that(sortPairs(c1.factor,c2.factor), equals(ref.object))

  cat("\n-different types-")
  ref.object <- list(spMat = NULL, levels = list(c1 = unique(c1.char), c2 = unique(c2.factor)), nij = rep(3,3), ni. =  rep(3,3),  n.j = rep(3,3))
  expect_that(sortPairs(c1.char,c2.factor), equals(ref.object))
})

test_that("Testing error for not handled input types", {

  c1 <- rep(1:3,each=3) * 1.5
  c2 <- rep(4:6,each=3) * 2.7

  cat("\n-NA-")
  c2NA <- c2; c2NA[1] <- NA
  expect_error(sortPairs(c1,c2NA))

  cat("\n-different sizes-")
  expect_error(sortPairs(c1,c2[-1]))

  cat("\n-list type-")
  expect_error(sortPairs(as.list(c1),as.list(c2)))

  cat("\n-matrix type-")
  expect_error(sortPairs(as.matrix(c1),as.matrix(c2)))

  cat("\n-data.frame type-")
  expect_error(sortPairs(as.data.frame(c1),as.data.frame(c2)))

})
