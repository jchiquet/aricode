#' Sort Pairs
#'
#' A function to sort pairs of integers or factors and identify the pairs between two classifications
#'
#' @param c1 A vector of length $n$ with values between 0 and $N_1 < n$ representing the
#' first classification. Supported types: integer, numeric, or factor. Avoid character
#' vectors for better performance. Must not be a list.
#' @param c2 A vector of length $n$ with values between 0 and $N_2 < n$ representing the
#' second classification. Supported types: integer, numeric, or factor. Avoid character
#' vectors for better performance. Must not be a list.
#' @param spMat Logical. If \code{TRUE}, returns the contingency table as a
#' sparse matrix. Note: sparse encoding may be more computationally expensive
#' than the algorithm itself. Default is \code{FALSE}.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item \strong{spMat}: A sparsely encoded contingency matrix (only if \code{spMat = TRUE}).
#'   \item \strong{levels}: A list containing the retained levels for each classification.
#'   \item \strong{nij}: A vector of positive pair counts.
#'   \item \strong{ni., n.j}: Vectors of class counts for \code{c1} and \code{c2}, respectively.
#'   \item \strong{pair_c1, pair_c2}: Integer vectors specifying the classes in \code{c1} and \code{c2}
#'   corresponding to the counts in \code{nij}. These provide the row and column indices for
#'   the contingency matrix.
#' }
#' @details
#' Pair sorting, which is at the heart of computing all clustering comparison measures, has been carefully
#' optimized. Hence, even basic R operations (checking for the presence of NAs, type conversion, or
#' constructing a sparse contingency matrix as an output) have non-negligible cost compared to the pair
#' sorting itself. For optimal performance, please provide the vectors as integers or factors without any NAs.
#'
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' out <- sort_pairs(cl, iris$Species)
#'
#' @import Matrix
#' @export
sort_pairs <- function(c1, c2, spMat = FALSE) {
  stopifnot("c1 and c2 must have the same length." = length(c1) == length(c2))
  stopifnot("c1 and c2 must have positive length." = length(c1) > 0)

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
    c1 <- as.integer(c1)
    c2 <- as.integer(c2)
    ## get_rank is O(n) if max(c1)-min(c1) and max(c2)-min(c2) is of order n=length(c1)=length(c2)
    res1 <- get_rank(c1)
    res2 <- get_rank(c2)
    mylevels <- list(c1 = res1$index, c2 = res2$index)
    c1 <- res1$translated # here ranks are in [0, n)
    c2 <- res2$translated
  }

  out <- std_sort_pairs(c1, c2, length(mylevels$c1), length(mylevels$c2))

  if (spMat) {
    spOut <- Matrix::sparseMatrix(
      i = out$pair_c1,
      j = out$pair_c2,
      x = out$count_pair,
      dims = sapply(mylevels, length),
      dimnames = mylevels, index1 = FALSE
    )
  } else {
    spOut <- NULL
  }

  res <- list(
    spMat = spOut,
    levels = mylevels,
    nij = out$count_pair,
    ni. = out$count_c1,
    n.j = out$count_c2,
    pair_c1 = out$pair_c1,
    pair_c2 = out$pair_c2
  )
  res
}


#' Sort Pairs
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' A function to sort pairs of integers or factors and identify the pairs between two classifications
#'
#' Just a change in the function name: please use [sort_pairs()].
#'
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' out <- sortPairs(cl, iris$Species)
#' # ->
#' out <- sort_pairs(cl, iris$Species)
#' @keywords internal
#' @importFrom lifecycle badge deprecate_warn
#' @export
sortPairs <- function(c1, c2, spMat = FALSE) {
  lifecycle::deprecate_warn("1.1.0", "sortPairs()", "sort_pairs()")
  sort_pairs(c1, c2, spMat)
}

#' Adjusted Rand Index
#'
#' A function to compute the adjusted rand index between two classifications
#'
#' @inheritParams sort_pairs
#'
#' @param sorted_pairs optional output of function sort_pairs (if already computed).
#' If `NULL` (the default), will be called internally
#'
#' @return a scalar with the adjusted Rand index.
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' ARI(cl, iris$Species)
#' @export
ARI <- function(c1, c2, sorted_pairs = NULL) {
  if (is.null(sorted_pairs)) {
    sorted_pairs <- sort_pairs(c1, c2)
  }

  N <- length(c1)
  stot <- sum(choose(sorted_pairs$nij, 2), na.rm = TRUE)
  srow <- sum(choose(sorted_pairs$ni., 2), na.rm = TRUE)
  scol <- sum(choose(sorted_pairs$n.j, 2), na.rm = TRUE)

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
#' A function to compute the Rand index between two classifications
#'
#' @inheritParams ARI
#'
#' @return a scalar with the Rand index.
#' @seealso \code{\link{ARI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' RI(cl, iris$Species)
#' @export
RI <- function(c1, c2, sorted_pairs = NULL) {
  if (is.null(sorted_pairs)) {
    sorted_pairs <- sort_pairs(c1, c2)
  }

  N <- length(c1)
  res <- 1 + (sum(sorted_pairs$nij^2) - (sum(sorted_pairs$ni.^2) + sum(sorted_pairs$n.j^2)) / 2) / choose(N, 2)
  res
}

#' Modified Adjusted Rand Index
#'
#' A function to compute a modified adjusted rand index between two classifications as proposed by Sundqvist et al. (2023), based on a multinomial model.
#'
#' @inheritParams ARI
#' @param raw Boolean: should the raw version of the MARI be computed? Default to `FALSE`.
#' @return a scalar with the modified ARI.
#' @seealso \code{\link{ARI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @references Sundqvist, Martina, Julien Chiquet, and Guillem Rigaill. "Adjusting the
#'  adjusted Rand Index: A multinomial story." Computational Statistics 38.1
#'  (2023): 327-347.
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' MARI(cl, iris$Species)
#' @export
MARI <- function(c1, c2, sorted_pairs = NULL, raw = FALSE) {
  if (is.null(sorted_pairs)) {
    sorted_pairs <- sort_pairs(c1, c2)
  }

  N <- length(c1)
  stot <- sum(choose(sorted_pairs$nij, 2), na.rm = TRUE)
  srow <- sum(choose(sorted_pairs$ni., 2), na.rm = TRUE)
  scol <- sum(choose(sorted_pairs$n.j, 2), na.rm = TRUE)

  ## using Lemma 3.3
  ## triplets
  T1 <- 2 * N
  T2 <- sum(as.double(sorted_pairs$nij) * as.double(sorted_pairs$ni.[sorted_pairs$pair_c1 + 1]) * as.double(sorted_pairs$n.j[sorted_pairs$pair_c2 + 1]), na.rm = TRUE)
  T3 <- -sum(sorted_pairs$nij^2, na.rm = TRUE) - sum(sorted_pairs$ni.^2, na.rm = TRUE) - sum(sorted_pairs$n.j^2, na.rm = TRUE)

  ## quadruplets (and division by 6 choose(N, 4)
  expectedIndex <- (srow * scol - stot - (T1 + T2 + T3)) / (6 * choose(N, 4))

  ## return the rand-index
  if (raw) {
    res <- (stot / choose(N, 2)) - expectedIndex
  } else {
    expectedIndex <- expectedIndex * choose(N, 2) ## RESCALE SO THAT THE CODE IS EQUIVALENT TO THE ARI
    maximumIndex <- (srow + scol) / 2
    if (expectedIndex == maximumIndex & stot != 0) {
      res <- 1
    } else {
      res <- (stot - expectedIndex) / (maximumIndex - expectedIndex)
    }
  }
  res
}

#' raw Modified Adjusted Rand Index
#'
#' A function to compute a modified adjusted rand index between two classifications as proposed by Sundqvist (2023), based on a multinomial model.
#' Raw means that the index is not divided by the (maximum - expected) value.
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' The function MARI now owns an argument function `raw` if one wishes to compute the raw version of MARI.
#'
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' out <- MARIraw(cl, iris$Species)
#' # ->
#' out <- MARI(cl, iris$Species, raw = TRUE)
#' @keywords internal
#' @importFrom lifecycle badge deprecate_warn
#' @export
MARIraw <- function(c1, c2, sorted_pairs = NULL) {
  lifecycle::deprecate_warn("1.1.0", "clustComp()", "compare_clustering()")
  MARI(c1, c2, sorted_pairs, raw = TRUE)
}


#' Chi-square statistics
#'
#' A function to compute the Chi-2 statistic
#'
#' @inheritParams ARI
#' @return a scalar with the Chi-square statistic.
#' @seealso \code{\link{ARI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' Chi2(cl, iris$Species)
#' @export
Chi2 <- function(c1, c2, sorted_pairs = NULL) {
  if (is.null(sorted_pairs)) {
    sorted_pairs <- sort_pairs(c1, c2)
  }
  N <- length(c1)
  res <- N * sum(sorted_pairs$nij^2 / (as.double(sorted_pairs$ni.[sorted_pairs$pair_c1 + 1]) * as.double(sorted_pairs$n.j[sorted_pairs$pair_c2 + 1])))
  res <- res - N
  res
}


#' Frobenius norm
#'
#' A function to compute the Frobenius norm between two classifications as defined in Lajugie et al. 2014 and Arlot et al 2019
#'
#' @inheritParams ARI
#' @return a scalar with the Frobenius norm.
#' @seealso \code{\link{ARI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @references
#'   - Rémi Lajugie, Francis Bach, and Sylvain Arlot. "Large-margin metric learning for constrained partitioning problems." International Conference on Machine Learning. PMLR, 2014.
#'   - Sylvain Arlot , Alain Celisse, and Zaid Harchaoui. "A kernel multiple change-point algorithm via model selection." Journal of machine learning research 20.162 (2019): 1-56.
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' Frobenius(cl, iris$Species)
#' @export
Frobenius <- function(c1, c2, sorted_pairs = NULL) {
  if (is.null(sorted_pairs)) {
    sorted_pairs <- sort_pairs(c1, c2)
  }
  out <- length(sorted_pairs$ni.) + length(sorted_pairs$n.j) - 2 * sum(sorted_pairs$nij^2 / (as.double(sorted_pairs$ni.[sorted_pairs$pair_c1 + 1]) * as.double(sorted_pairs$n.j[sorted_pairs$pair_c2 + 1])))
  out
}

#' Entropy
#'
#' A function to compute the empirical entropy for two vectors of classification and the joint entropy
#'
#' @inheritParams ARI
#' @return a list with the two conditional entropies, the joint entropy and output of sort_pairs.
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' entropy(cl, iris$Species)
#' @export
entropy <- function(c1, c2, sorted_pairs = NULL) {
  if (is.null(sorted_pairs)) {
    sorted_pairs <- sort_pairs(c1, c2)
  }

  N <- length(c1)

  H.UV <- -sum(sorted_pairs$nij * log(sorted_pairs$nij)) / N + log(N)
  H.U <- -sum(sorted_pairs$ni. * log(sorted_pairs$ni.)) / N + log(N)
  H.V <- -sum(sorted_pairs$n.j * log(sorted_pairs$n.j)) / N + log(N)

  res <- list(UV = H.UV, U = H.U, V = H.V, sort_pairs = sorted_pairs)
  res
}

#' Normalized mutual information (NMI)
#'
#' A function to compute the NMI between two classifications
#'
#' @inheritParams ARI
#' @param variant a string in ("max", "min", "sqrt", "sum", "joint"): different variants of NMI. Default use "max".
#' @return a scalar with the normalized mutual information .
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{ARI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' NMI(cl, iris$Species)
#' @export
NMI <- function(c1, c2, variant = c("max", "min", "sqrt", "sum", "joint"), sorted_pairs = NULL) {
  variant <- match.arg(variant)

  if (is.null(sorted_pairs)) {
    sorted_pairs <- sort_pairs(c1, c2)
  }

  H <- entropy(c1, c2, sorted_pairs)
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
#' @inheritParams ARI
#' @return a scalar with the normalized information distance .
#' @seealso \code{\link{RI}}, \code{\link{NMI}}, \code{\link{NVI}}, \code{\link{ARI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' NID(cl, iris$Species)
#' @export
NID <- function(c1, c2, sorted_pairs = NULL) {
  if (is.null(sorted_pairs)) {
    sorted_pairs <- sort_pairs(c1, c2)
  }

  H <- entropy(c1, c2, sorted_pairs)
  MI <- -H$UV + H$U + H$V
  res <- 1 - MI / max(H$U, H$V)
  res
}

#' Normalized variation of information (NVI)
#'
#' A function to compute the NVI between two classifications
#'
#' @inheritParams ARI
#' @return a scalar with the normalized variation of information.
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NMI}}, \code{\link{ARI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' NVI(cl, iris$Species)
#' @export
NVI <- function(c1, c2, sorted_pairs = NULL) {
  if (is.null(sorted_pairs)) {
    sorted_pairs <- sort_pairs(c1, c2)
  }

  H <- entropy(c1, c2, sorted_pairs)
  MI <- -H$UV + H$U + H$V
  res <- 1 - MI / H$UV
  res
}

#' Adjusted Mutual Information
#'
#' A function to compute the adjusted mutual information between two classifications
#'
#' @inheritParams ARI
#' @return a scalar with the adjusted rand index.
#' @seealso \code{\link{ARI}}, \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' AMI(cl, iris$Species)
#' @export
AMI <- function(c1, c2, sorted_pairs = NULL) {
  if (is.null(sorted_pairs)) {
    sorted_pairs <- sort_pairs(c1, c2)
  }

  H <- entropy(c1, c2, sorted_pairs)
  MI <- -H$UV + H$U + H$V
  EMI <- expected_MI(as.integer(H$sort_pairs$ni.), as.integer(H$sort_pairs$n.j))

  res <- (MI - EMI) / (max(H$U, H$V) - EMI)
  res
}

#' Measures of similarity between two classification
#'
#' A function for computing all the measures of similarity implemented in this package at once.
#' Include (A)RI, (N)MI, (N)VI, (N)ID, Chi2, MARI, Frobenius
#'
#' @inheritParams ARI
#' @param AMI Boolean: should the AMI be computed (more costly than all other measures)? Default is `FALSE`.
#' @return a list with all the measures available
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{ARI}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' compare_clustering(cl, iris$Species)
#' @export
compare_clustering <- function(c1, c2, sorted_pairs = NULL, AMI = FALSE) {
  if (is.null(sorted_pairs)) {
    sorted_pairs <- sort_pairs(c1, c2)
  }

  H <- entropy(c1, c2, sorted_pairs)

  MI <- -H$UV + H$U + H$V
  VI <- H$UV - MI
  NVI <- 1 - MI / H$UV
  ID <- max(H$U, H$V) - MI
  NID <- 1 - MI / max(H$U, H$V)
  NMI <- MI / max(H$U, H$V)

  res <- list(
    RI = RI(c1, c2, sorted_pairs),
    ARI = ARI(c1, c2, sorted_pairs),
    MI = -H$UV + H$U + H$V,
    NMI = MI / max(H$U, H$V),
    VI = H$UV - MI,
    NVI = 1 - MI / H$UV,
    ID = max(H$U, H$V) - MI,
    NID = 1 - MI / max(H$U, H$V),
    Chi2 = Chi2(c1, c2, sorted_pairs),
    MARI = MARI(c1, c2, sorted_pairs),
    MARIraw = MARIraw(c1, c2, sorted_pairs),
    Frobenius = Frobenius(c1, c2, sorted_pairs)
  )
  if (AMI) {
    EMI <- expected_MI(as.integer(H$sort_pairs$ni.), as.integer(H$sort_pairs$n.j))
    res$AMI <- (-H$UV + H$U + H$V - EMI) / (max(H$U, H$V) - EMI)
  }
  res
}

#' Measures of similarity between two classification
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' A function for computing all the measures of similarity implemented in this package at once.
#' Include ARI, RI, MI VI, NVI, ID, NID, NMI, Chi2, MARI, Frobenius
#'
#' Just a change in the function name: please use [compare_clustering()].
#'
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[, -5])), 4)
#' out <- sort_pairs(cl, iris$Species)
#' # ->
#' out <- sort_pairs(cl, iris$Species)
#' @keywords internal
#' @importFrom lifecycle badge deprecate_warn
#' @export
clustComp <- function(c1, c2, sorted_pairs = NULL, AMI = FALSE) {
  lifecycle::deprecate_warn("1.1.0", "clustComp()", "compare_clustering()")
  compare_clustering(c1, c2, sorted_pairs, AMI)
}
