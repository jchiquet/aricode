library(testthat)
library(aricode)

## define the ARI as in the mclust package
adjustedRandIndex <- function (x, y)
{
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y))
    stop("arguments must be vectors of the same length")
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1)))
    return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}

test_that("Testing coherence of the adjusted Rand Index", {

  cat("\n-large random vectors- ")
  n <- 1e5
  c1 <- as.numeric(sample(1:(n/100), n, replace=TRUE))
  c2 <- as.numeric(sample(1:(n/100), n, replace=TRUE))
  expect_equal(ARI(c1,c2), adjustedRandIndex(c1,c2))

  cat("\n-real data- ")
  data(iris)
  cl <- cutree(hclust(dist(iris[,-5])), 4)
  expect_equal(ARI(cl,iris$Species), adjustedRandIndex(cl,iris$Species))

  cat("\n-completely equal vectors with no groups-")
  c1 <- 1:100
  c2 <- 1:100
  expect_equal(ARI(c1,c2), adjustedRandIndex(c1,c2))

  cat("\n-completely equal vectors with one groups-")
  c1 <- rep(1,100)
  c2 <- rep(2,100)
  expect_equal(ARI(c1,c2), adjustedRandIndex(c1,c2))

  cat("\n-completely different vectors with one groups-")
  c1 <- c(rep(0,99),1)
  c2 <- rep(1,100)
  expect_equal(ARI(c1,c2), adjustedRandIndex(c1,c2))
})
