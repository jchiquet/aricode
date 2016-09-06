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

n <- 10000
c1 <- as.numeric(sample(1:(n/100), n, replace=TRUE))
c2 <- as.numeric(sample(1:(n/100), n, replace=TRUE))

test_that("Testing cohenrence with mclust function adjustedRandIndex", {
  expect_equal(ARI(c1,c2), adjustedRandIndex(c1,c2))
})

