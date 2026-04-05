#' Sort Pairs
#'
#' A function to sort pairs of integers or factors and identify the pairs between two classifications
#'
#' @param c1 a vector of length n with value between 0 and N1 < n containing the first classification.
#' Must be a vector of integers or characters, a numeric, a factor, but not a list. Avoid character.
#' @param c2 a vector of length n with value between 0 and N2 < n containing the second classification.
#' Must be a vector of integers or characters, a numeric or a factor, but not a list. Avoid character.
#' @param spMat logical: send back the contingency table as sparsely encoded
#' (cost more than the algorithm itself). Default is FALSE
#' @details
#' Pair sorting, which is at the heart of computing all clustering comparison measures, has been carefully
#' optimized. Hence, even basic R operations (checking for the presence of NAs, type conversion, or
#' constructing a sparse contingency matrix as an output) have non negligible cost compare to the pair
#' sorting itself. For optimal performance, please provide the vectors as integers or factors without any NAs.
#'
#' @import Matrix
#' @export
sortPairs <- function(c1, c2, spMat = FALSE) {
  stopifnot("c1 and c2 must have the same length." = length(c1) == length(c2))

  stopifnot(
    "c1 and c2 must be vectors or factors but not lists." =
      (is.vector(c1) || is.factor(c1)) && !is.list(c1) &&
        (is.vector(c2) || is.factor(c2)) && !is.list(c2)
  )

  if (is.character(c1) || is.character(c2)) {
    if (length(c1) > 1e4 || length(c2) > 1e4) {
      warning("Converting a long vector of characters to integers is much slower than the pair sorting itself. Consider using integers or factors, or converting them beforehand for better performance")
    }
    c1 <- as.integer(c1)
    c2 <- as.integer(c2)
  }

  if (anyNA(c1) || anyNA(c2)) stop("NA are not supported.")

  if (is.factor(c1) && is.factor(c2)) {
    mylevels <- list(c1 = levels(c1), c2 = levels(c2))
    c1 <- as.integer(c1) - 1L
    c2 <- as.integer(c2) - 1L
  } else {
    ## getRank is O(n) if max(c1)-min(c1) and max(c2)-min(c2) is of order length(c1)=length(c2)
    ## NOTE: getRank does not assume c1 and c2 are between 0 and n
    c1 <- as.integer(c1)
    c2 <- as.integer(c2)
    res1 <- getRank(c1)
    res2 <- getRank(c2)
    mylevels <- list(c1 = res1$index, c2 = res2$index)
    c1 <- res1$translated # here ranks are in [0, n)
    c2 <- res2$translated # here ranks are in [0, n)
  }

  out <- std_SortPairs(c1, c2, length(mylevels$c1), length(mylevels$c2))

  if (spMat) {
    spOut <- sparseMatrix(
      i = out$pair_c1,
      j = out$pair_c2,
      x = out$pair_nb,
      dims = sapply(mylevels, length),
      dimnames = mylevels, index1 = FALSE
    )
  } else {
    spOut <- NULL
  }

  res <- list(
    spMat = spOut,
    levels = mylevels,
    nij = out$pair_nb,
    ni. = out$c1_nb,
    n.j = out$c2_nb,
    pair_c1 = out$pair_c1,
    pair_c2 = out$pair_c2
  )
  res
}

#' Adjusted Rand Index
#'
#' A function to compute the adjusted rand index between two classifications
#'
#' @inheritParams sortPairs
#'
#' @return a scalar with the adjusted Rand index.
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' ARI(cl, iris$Species)
#' @export
ARI <- function(c1, c2) {
  ## get pairs using C
  ## ensure that values of c1 and c2 are between 0 and n1
  res <- sortPairs(c1, c2)

  ## get ARI using pairs
  N <- length(c1)

  stot <- sum(choose(res$nij, 2), na.rm = TRUE)
  srow <- sum(choose(res$ni., 2), na.rm = TRUE)
  scol <- sum(choose(res$n.j, 2), na.rm = TRUE)

  expectedIndex <- (srow * scol) / (choose(N, 2))
  maximumIndex <- (srow + scol) / 2

  if (expectedIndex == maximumIndex && stot != 0) {
    res <- 1
  } else {
    res <- (stot - expectedIndex) / (maximumIndex - expectedIndex)
  }
  res
}

#' Rand Index
#'
#' A function to compute the rand index between two classifications
#'
#' @inheritParams sortPairs
#' @return a scalar with the rand index.
#' @seealso \code{\link{ARI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' RI(cl, iris$Species)
#' @export
RI <- function(c1, c2) {
  ## get pairs using C
  ## ensure that values of c1 and c2 are between 0 and n1
  res <- sortPairs(c1, c2)

  ## get ARI using pairs
  N <- length(c1)

  ## return the rand-index
  res <- 1 + (sum(res$nij^2) - (sum(res$ni.^2) + sum(res$n.j^2)) / 2) / choose(N, 2)
  res
}

#' Modified Adjusted Rand Index
#'
#' A function to compute a modified adjusted rand index between two classifications as proposed by Sundqvist et al. in prep, based on a multinomial model.
#'
#' @inheritParams sortPairs
#' @return a scalar with the modified ARI.
#' @seealso \code{\link{ARI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' MARI(cl, iris$Species)
#' @export
MARI <- function(c1, c2) {
  ## get pairs using C
  ## ensure that values of c1 and c2 are between 0 and n1
  res <- sortPairs(c1, c2)
  N <- length(c1)

  stot <- sum(choose(res$nij, 2), na.rm = TRUE)
  srow <- sum(choose(res$ni., 2), na.rm = TRUE)
  scol <- sum(choose(res$n.j, 2), na.rm = TRUE)

  ## using Lemma 3.3
  ## triplets
  T1 <- 2 * N
  T2 <- sum(res$nij * res$ni.[res$pair_c1 + 1] * res$n.j[res$pair_c2 + 1], na.rm = TRUE)
  T3 <- -sum(res$nij^2, na.rm = TRUE) - sum(res$ni.^2, na.rm = TRUE) - sum(res$n.j^2, na.rm = TRUE)

  ## quadruplets (and division by 6 choose(N, 4)
  expectedIndex <- (srow * scol - stot - (T1 + T2 + T3)) / (6 * choose(N, 4))

  ## return the rand-index
  expectedIndex <- expectedIndex * choose(N, 2) ## RESCALE SO THAT THE CODE IS EQUIVALENT TO THE ARI
  maximumIndex <- (srow + scol) / 2
  if (expectedIndex == maximumIndex & stot != 0) {
    res <- 1
  } else {
    res <- (stot - expectedIndex) / (maximumIndex - expectedIndex)
  }
  res
  ## return the adjusted (and divided) rand-index
  res
}

#' raw Modified Adjusted Rand Index
#'
#' A function to compute a modified adjusted rand index between two classifications as proposed by Sundqvist et al. in prep, based on a multinomial model. Raw means, that the index is not divided by the (maximum - expected) value.
#'
#' @inheritParams sortPairs
#' @return a scalar with the modified ARI without the division by the (maximum - expected)
#' @seealso \code{\link{ARI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' MARIraw(cl, iris$Species)
#' @export
MARIraw <- function(c1, c2) {
  ## get pairs using C
  ## ensure that values of c1 and c2 are between 0 and n1
  res <- sortPairs(c1, c2)
  N <- length(c1)

  stot <- sum(choose(res$nij, 2), na.rm = TRUE)
  srow <- sum(choose(res$ni., 2), na.rm = TRUE)
  scol <- sum(choose(res$n.j, 2), na.rm = TRUE)

  ## using Lemma 3.3
  ## triplets
  T1 <- 2 * N
  T2 <- sum(res$nij * res$ni.[res$pair_c1 + 1] * res$n.j[res$pair_c2 + 1], na.rm = TRUE)
  T3 <- -sum(res$nij^2, na.rm = TRUE) - sum(res$ni.^2, na.rm = TRUE) - sum(res$n.j^2, na.rm = TRUE)

  ## quadruplets (and division by 6 choose(N, 4)
  expectedIndex <- (srow * scol - stot - (T1 + T2 + T3)) / (6 * choose(N, 4))

  ## return the rand-index
  res <- (stot / choose(N, 2)) - expectedIndex
  res
}

#' Chi-square statistics
#'
#' A function to compute the Chi-2 statistics
#'
#' @inheritParams sortPairs
#' @return a scalar with the chi-square statistics.
#' @seealso \code{\link{ARI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' Chi2(cl, iris$Species)
#' @export
Chi2 <- function(c1, c2) {
  ## get pairs using C
  ## ensure that values of c1 and c2 are between 0 and n1
  res <- sortPairs(c1, c2)
  N <- length(c1)

  res <- N * sum(res$nij^2 / (res$ni.[res$pair_c1 + 1] * res$n.j[res$pair_c2 + 1]))
  res <- res - N
  res
}


#' Frobenius norm
#'
#' A function to compute the Frobenius norm between two classification as defined in Lajugie et al. 2014 and Arlot et al 2019
#'
#' @inheritParams sortPairs
#' @return a scalar with the chi-square statistics.
#' @seealso \code{\link{ARI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @references
#'   - Rémi Lajugie, Francis Bach, and Sylvain Arlot. "Large-margin metric learning for constrained partitioning problems." International Conference on Machine Learning. PMLR, 2014.
#'   - Sylvain Arlot , Alain Celisse, and Zaid Harchaoui. "A kernel multiple change-point algorithm via model selection." Journal of machine learning research 20.162 (2019): 1-56.
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' Frobenius(cl, iris$Species)
#' @export
Frobenius <- function(c1, c2) {
  ## get pairs using C
  ## ensure that values of c1 and c2 are between 0 and n1
  res <- sortPairs(c1, c2)

  out <- length(res$ni.) + length(res$n.j) - 2 * sum(res$nij^2 / (res$ni.[res$pair_c1 + 1] * res$n.j[res$pair_c2 + 1]))
  out
}

#' Entropy
#'
#' A function to compute the empirical entropy for two vectors of classification and the joint entropy
#'
#' @inheritParams sortPairs
#' @return a list with the two conditional entropies, the joint entropy and output of sortPairs.
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' entropy(cl, iris$Species)
#' @export
entropy <- function(c1, c2) {
  res <- sortPairs(c1, c2)

  N <- length(c1)

  H.UV <- -sum(res$nij * log(res$nij)) / N + log(N)
  H.U <- -sum(res$ni. * log(res$ni.)) / N + log(N)
  H.V <- -sum(res$n.j * log(res$n.j)) / N + log(N)

  res <- list(UV = H.UV, U = H.U, V = H.V, sortPairs = res)
  res
}

#' Measures of similarity between two classification
#'
#' A function various measures of similarity between two classifications
#'
#' @inheritParams sortPairs
#' @return a list with all the measures available
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{ARI}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' clustComp(cl, iris$Species)
#' @export
clustComp <- function(c1, c2) {
  H <- entropy(c1, c2)

  MI <- -H$UV + H$U + H$V
  VI <- H$UV - MI
  NVI <- 1 - MI / H$UV
  ID <- max(H$U, H$V) - MI
  NID <- 1 - MI / max(H$U, H$V)
  NMI <- MI / max(H$U, H$V)
  EMI <- expected_MI(as.integer(H$sortPairs$ni.), as.integer(H$sortPairs$n.j))

  res <- list(
    RI = RI(c1, c2),
    ARI = ARI(c1, c2),
    MI = -H$UV + H$U + H$V,
    AMI = (-H$UV + H$U + H$V - EMI) / (max(H$U, H$V) - EMI),
    VI = H$UV - MI,
    NVI = 1 - MI / H$UV,
    ID = max(H$U, H$V) - MI,
    NID = 1 - MI / max(H$U, H$V),
    NMI = MI / max(H$U, H$V),
    Chi2 = Chi2(c1, c2),
    MARI = MARI(c1, c2),
    MARIraw = MARIraw(c1, c2),
    Frobenius = Frobenius(c1, c2)
  )
  res
  res
}

#' Adjusted Mutual Information
#'
#' A function to compute the adjusted mutual information between two classifications
#'
#' @inheritParams sortPairs
#' @return a scalar with the adjusted rand index.
#' @seealso \code{\link{ARI}}, \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' AMI(cl, iris$Species)
#' @export
AMI <- function(c1, c2) {
  H <- entropy(c1, c2)
  MI <- -H$UV + H$U + H$V
  EMI <- expected_MI(as.integer(H$sortPairs$ni.), as.integer(H$sortPairs$n.j))

  res <- (MI - EMI) / (max(H$U, H$V) - EMI)
  res
}

#' Normalized mutual information (NMI)
#'
#' A function to compute the NMI between two classifications
#'
#' @inheritParams sortPairs
#' @param variant a string in ("max", "min", "sqrt", "sum", "joint"): different variants of NMI. Default use "max".
#' @return a scalar with the normalized mutual information .
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{ARI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' NMI(cl, iris$Species)
#' @export
NMI <- function(c1, c2, variant = c("max", "min", "sqrt", "sum", "joint")) {
  variant <- match.arg(variant)

  H <- entropy(c1, c2)
  MI <- -H$UV + H$U + H$V

  D <- switch(variant,
    "max" = max(H$U, H$V),
    "sqrt" = sqrt(H$U * H$V),
    "min" = min(H$U, H$V),
    "sum" = .5 * (H$U + H$V),
    "joint" = H$UV
  )
  res <- MI / D
  res
}

#' Normalized information distance (NID)
#'
#' A function to compute the NID between two classifications
#'
#' @inheritParams sortPairs
#' @return a scalar with the normalized information distance .
#' @seealso \code{\link{RI}}, \code{\link{NMI}}, \code{\link{NVI}}, \code{\link{ARI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' NID(cl, iris$Species)
#' @export
NID <- function(c1, c2) {
  H <- entropy(c1, c2)
  MI <- -H$UV + H$U + H$V
  res <- 1 - MI / max(H$U, H$V)
  res
}

#' Normalized variation of information (NVI)
#'
#' A function to compute the NVI between two classifications
#'
#' @inheritParams sortPairs
#' @return a scalar with the normalized variation of information.
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NMI}}, \code{\link{ARI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' NVI(cl, iris$Species)
#' @export
NVI <- function(c1, c2) {
  H <- entropy(c1, c2)
  MI <- -H$UV + H$U + H$V
  res <- 1 - MI / H$UV
  res
}
