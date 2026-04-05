library(microbenchmark)
library(aricode)
library(mclust)

n <- 10^7
c1 <- sample.int(100, size = n, replace = TRUE)
c2 <- sample.int(100, size = n, replace = TRUE)
res1 <- aricode:::getRank(c1)
res2 <- aricode:::getRank(c2)
mylevels <- list(c1 = res1$index, c2 = res2$index)

## define the ARI as in the mclust package
adjusted_rand_index <- function(x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y)) {
    stop("arguments must be vectors of the same length")
  }
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1))) {
    return(1)
  }
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c) / (a + b + c + d)) /
    ((a + b + a + c) / 2 - (a + b) * (a + c) / (a + b + c + d))
  return(ARI)
}

res <- microbenchmark(
  aricode = aricode::ARI(c1, c2),
  R = adjusted_rand_index(c1, c2),
  mclust = mclust::adjustedRandIndex(c1, c2), times = 20L
)
ggplot2::autoplot(res)


res2 <- microbenchmark::microbenchmark(
  aricode:::std_SortPairs(c1, c2, length(mylevels$c1), length(mylevels$c2)),
  aricode:::cpp_SortPairs(c1, c2, length(mylevels$c1), length(mylevels$c2)),
  times = 20L
)
ggplot2::autoplot(res2)

c1 <- sample(1:(n / 200), n, replace = TRUE)
c2 <- c1
i_change <- sample(1:n, n / 50, replace = FALSE)
c2[i_change] <- c2[rev(i_change)]
profvis::profvis(sortPairs(c1, c2))

c1 <- as.numeric(c1)
c2 <- as.numeric(c2)
profvis::profvis(sortPairs(c1, c2))

c1 <- as.character(c1)
c2 <- as.character(c2)
profvis::profvis(sortPairs(c1, c2))
