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

EMI <- function(counts) {
    s_emi <- 0
    n <- sum(counts)
    s1 <- margin.table(counts, 1)
    s2 <- margin.table(counts, 2)
    for(i in 1:nrow(counts)){
        for (j in 1:ncol(counts)){
            ai <- s1[i]
            bj <- s2[j]
            min_nij <- max(1, ai+bj-n)
            max_nij <- min(ai,bj)
            if (min_nij >= max_nij) next
            n.ij <- seq(min_nij, max_nij)   #sequence of consecutive numbers
            t1<- (n.ij / n) * log((n.ij * n) / (ai*bj))
            t2 <- exp(lfactorial(ai) + lfactorial(bj) + lfactorial(n - ai) + lfactorial(n - bj) - lfactorial(n) - lfactorial(n.ij) - lfactorial(ai - n.ij) - lfactorial(bj - n.ij) - lfactorial(n - ai - bj + n.ij))
            emi <- sum(t1*t2)
            if (is.nan(emi)) browser()
            s_emi <- s_emi + emi
        }
    }
    s_emi
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

test_that("Testing coherence of the expected mutual information", {

  cat("\n-real data- ")
  data(iris)
  cl <- cutree(hclust(dist(iris[,-5])), 4)
  counts <- table(cl,iris$Species)
  ni. <- rowSums(counts)
  n.j <- colSums(counts)
  expect_equal(aricode:::expected_MI(ni., n.j), EMI(counts))

  cat("\n-completely equal vectors with one groups-")
  c1 <- rep(1,100)
  c2 <- rep(2,100)
  counts <- table(c1, c2)
  ni. <- rowSums(counts)
  n.j <- colSums(counts)
  expect_equal(aricode:::expected_MI(ni., n.j), EMI(counts))

  cat("\n-completely different vectors with one groups-")
  c1 <- c(rep(0,99),1)
  c2 <- rep(1,100)
  counts <- table(c1, c2)
  ni. <- rowSums(counts)
  n.j <- colSums(counts)
  expect_equal(aricode:::expected_MI(ni., n.j), EMI(counts))
})
