library(microbenchmark)
library(aricode)

n <- 10^6
c1 <- sample.int(100, size = n, replace=TRUE)
c2 <- sample.int(100, size = n, replace=TRUE)
res_old <- sortPairs(c1, c2)
res_new <- sortPairs_new(c1, c2)

sum(abs(res_new$nij - res_old$nij))
sum(abs(res_new$ni. - res_old$ni.[-1]))
sum(abs(res_new$n.j - res_old$nj.[-1]))
res1 <- microbenchmark(order(c1, c2), sortPairs(c1, c2), times=100L)
ggplot2::autoplot(res1)

n <- 10^6
c1 <- sample.int(1000, size = n, replace = TRUE)
c2 <- sample.int(1000, size = n, replace = TRUE)

res2 <- microbenchmark(order(c1, c2), sortPairs(c1, c2), times=100L)
ggplot2::autoplot(res2)
