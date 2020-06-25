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

Chi2_ref <- function(c1, c2) {
  as.numeric(chisq.test(c1, c2, correct=F)$stat[1])
}

## Martina's code
ARI_estimated <- function(c1, c2) {
  # c1, c2, two classifications of the same observations
  n_kl <- table(c1,c2)

  n   <- sum(n_kl)
  n_k <- rowSums(n_kl)
  n_l <- colSums(n_kl)

  T1 <- sum(choose(n_kl,2))
  T2 <- 2*n + sum(outer(n_k, n_l) * n_kl) - sum(n_kl^2) - sum(n_k^2) -  sum(n_l^2)
  T3 <- 1/(6*choose(n,4))*(sum(outer(n_k^2, n_l^2)) - 4*T1 - 4*T2 - n*2*(sum(choose(n_k,2)) + sum(choose(n_l,2))) - n^2)/4

  RI_estim <- (1/choose(n, 2))*T1
  ARI_estim <- RI_estim-T3
 ARI_estim
}

test_that("Testing coherence of the adjusted Rand Index", {

  ## "\n-large random vectors-
  n <- 1e5
  c1 <- as.numeric(sample(1:(n/100), n, replace=TRUE))
  c2 <- as.numeric(sample(1:(n/100), n, replace=TRUE))
  expect_equal(ARI(c1,c2), adjustedRandIndex(c1,c2))

  ## "\n-real data-
  data(iris)
  cl <- cutree(hclust(dist(iris[,-5])), 4)
  expect_equal(ARI(cl,iris$Species), adjustedRandIndex(cl,iris$Species))

  ## "\n-completely equal vectors with no groups-")
  c1 <- 1:100
  c2 <- 1:100
  expect_equal(ARI(c1,c2), adjustedRandIndex(c1,c2))

  ## "\n-completely equal vectors with one groups-")
  c1 <- rep(1,100)
  c2 <- rep(2,100)
  expect_equal(ARI(c1,c2), adjustedRandIndex(c1,c2))

  ## "\n-completely different vectors with one groups-")
  c1 <- c(rep(0,99),1)
  c2 <- rep(1,100)
  expect_equal(ARI(c1,c2), adjustedRandIndex(c1,c2))
})

test_that("Testing coherence of the expected mutual information", {

  ## "\n-real data-
  data(iris)
  cl <- cutree(hclust(dist(iris[,-5])), 4)
  counts <- table(cl,iris$Species)
  ni. <- rowSums(counts)
  n.j <- colSums(counts)
  expect_equal(aricode:::expected_MI(ni., n.j), EMI(counts))

  ## "\n-completely equal vectors with one groups-")
  c1 <- rep(1,100)
  c2 <- rep(2,100)
  counts <- table(c1, c2)
  ni. <- rowSums(counts)
  n.j <- colSums(counts)
  expect_equal(aricode:::expected_MI(ni., n.j), EMI(counts))

  ## "\n-completely different vectors with one groups-")
  c1 <- c(rep(0,99),1)
  c2 <- rep(1,100)
  counts <- table(c1, c2)
  ni. <- rowSums(counts)
  n.j <- colSums(counts)
  expect_equal(aricode:::expected_MI(ni., n.j), EMI(counts))
})


test_that("Testing coherence of the Chi2 statistics information", {

  ## "\n-large random vectors-
  n <- 1e5
  c1 <- as.numeric(sample(1:(n/100), n, replace=TRUE))
  c2 <- as.numeric(sample(1:(n/100), n, replace=TRUE))
  expect_equal(ARI(c1,c2), adjustedRandIndex(c1,c2))

  ## "\n-moderate random-
  n <- rpois(lambda=100, n=1) +   3
  k1 <- rpois(lambda=5, n=1)+2; k2 <- rpois(lambda=5, n=1)+2
  c1 <- sample(1:k1, n, replace=T)
  c2 <- sample(1:k2, n, replace=T)
  expect_equal(aricode::Chi2(c1,c2), Chi2_ref(c1,c2))

  ## "\n-real data-
  data(iris)
  cl <- cutree(hclust(dist(iris[,-5])), 4)
  expect_equal(aricode::Chi2(cl,iris$Species), Chi2_ref(cl,iris$Species))

  ## "\n-completely equal vectors with no groups-")
  c1 <- 1:100
  c2 <- 1:100
  expect_equal(aricode::Chi2(c1,c2), Chi2_ref(c1,c2))
})

test_that("Testing coherence of the MARI", {

  ## "\n-large random vectors-
  n <- 1e5
  c1 <- as.numeric(sample(1:(n/100), n, replace=TRUE))
  c2 <- as.numeric(sample(1:(n/100), n, replace=TRUE))
  expect_equal(aricode::MARIraw(c1,c2), ARI_estimated(c1,c2))

  ## "\n-moderate random-
  n <- rpois(lambda=100, n=1) +   3
  k1 <- rpois(lambda=5, n=1)+2; k2 <- rpois(lambda=5, n=1)+2
  c1 <- sample(1:k1, n, replace=T)
  c2 <- sample(1:k2, n, replace=T)
  expect_equal(aricode::MARIraw(c1,c2), ARI_estimated(c1,c2))

  ## "\n-real data-
  data(iris)
  cl <- cutree(hclust(dist(iris[,-5])), 4)
  expect_equal(aricode::MARIraw(c1,c2), ARI_estimated(c1,c2))

  ## "\n-completely equal vectors with no groups-")
  c1 <- 1:100
  c2 <- 1:100
  expect_equal(aricode::MARIraw(c1,c2), ARI_estimated(c1,c2))

  ## "\n-completely equal vectors with one groups-")
  c1 <- rep(1,100)
  c2 <- rep(2,100)
  expect_equal(aricode::MARIraw(c1,c2), ARI_estimated(c1,c2))

  ## "\n-completely different vectors with one groups-")
  c1 <- c(rep(0,99),1)
  c2 <- rep(1,100)
  expect_equal(aricode::MARIraw(c1,c2), ARI_estimated(c1,c2))

})

test_that("Testing coherence of the MARI with ARI when equivalent", {
  n <- rpois(lambda=100, n=1) +   3
  k1 <- rpois(lambda=5, n=1)+2; k2 <- rpois(lambda=5, n=1)+2
  class1 <- sample(1:k1, n, replace=T)
  class2 <- sample(1:k2, n, replace=T)
  ## expect_equal or almost equal in test_that
  expect_equal(ARI(class1, class2), MARI(class1, class2), tolerance = 1e-3)
})
