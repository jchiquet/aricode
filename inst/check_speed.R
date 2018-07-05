library(microbenchmark)
library(aricode)
library(mclust)

n <- 10^5
c1 <- sample.int(100, size = n, replace=TRUE)
c2 <- sample.int(100, size = n, replace=TRUE)

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

res <- microbenchmark(aricode = ARI(c1, c2),
                      R = adjustedRandIndex(c1, c2),
                      mclust = mclust::adjustedRandIndex(c1, c2), times=100L)
ggplot2::autoplot(res)

