#' @title Sort Pairs
#' @description A function to sort pairs of integers or factors and identify the pairs
#'
#' @param c1 a vector of length n with value between 0 and N1 < n
#' @param c2 a vector of integer of length n with value between 0 and N2 < n
#' @param spMat logical: send back the contigency table a sparsely encoded (cost more than the algorithm itself). Default is FALSE
#' @import Matrix
#' @useDynLib aricode
#' @export
sortPairs <- function(c1, c2, spMat = FALSE){

  if (anyNA(c1) | anyNA(c2))
    stop("NA are not supported.")

  if (( (!is.vector(c1) & !is.factor(c1)) | is.list(c1)) | ((!is.vector(c2) & !is.factor(c2)) | is.list(c2)))
    stop("c1 and c2 must be vectors or factors but not lists.")

  if (length(c1) != length(c2))
    stop("the two vectors must have the same length.")

  n <- length(c1)

  ## if c1 and c2 are integer
  if (is.integer(c1) & is.integer(c2)) {
    mylevels <- list(c1 = unique(c1), c2 = unique(c2))
    c1 <- c1 - min(c1)
    c2 <- c2 - min(c2)
    ## if the range is not adapted to the C code
    if (!(max(c1) <= n-1 & max(c2) <= n-1)) {
      c1 <- as.integer(factor(c1, levels = mylevels$c1)) - 1L
      c2 <- as.integer(factor(c2, levels = mylevels$c2)) - 1L
    }
  ## if factor, force conversion to integer
  } else if (is.factor(c1) & is.factor(c2)) {
    mylevels <- list(c1 = levels(c1), c2 = levels(c2))
    c1 <- as.integer(c1) - 1L
    c2 <- as.integer(c2) - 1L
  } else {
    ## if neither a factor nor an integer or different of types force to factor then integer
    mylevels <- list(c1 = unique(c1), c2 = unique(c2))
    c1 <- as.integer(factor(c1, levels = mylevels$c1)) - 1L
    c2 <- as.integer(factor(c2, levels = mylevels$c2)) - 1L
  }

  result <- .C("c_SortPairs",
    c1         = as.integer(c1),
    c2         = as.integer(c2),
    new_c1     = integer(n),
    new_c2     = integer(n),
    pair_c1    = integer(n),
    pair_c2    = integer(n),
    pair_count = integer(n),
    count1     = integer(n),
    count2     = integer(n),
    n          = as.integer(n),
    nzero      = as.integer(1), PACKAGE="aricode")

  if (spMat) {
    spOut <- sparseMatrix(i=result$pair_c1[1:(result$nzero+1)],
                          j=result$pair_c2[1:(result$nzero+1)],
                          x=result$pair_count[1:(result$nzero+1)],
                          dims=sapply(mylevels,length),
                          dimnames = mylevels, index1=FALSE)
  } else {
    spOut <- NULL
  }

  return(list(spMat = spOut,
              levels = mylevels,
              nij = result$pair_count[1:(result$nzero+1)],
              ni. = result$count1[which(result$count1 > 0)],
              n.j = result$count2[which(result$count2 > 0)])
  )
}

#' @title Adjusted Rand Index
#' @description A function to compute the adjusted rand index between two classifications
#'
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @return a scalar with the adjusted rand index.
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[,-5])), 4)
#' ARI(cl,iris$Species)
#' @export
ARI <- function(c1, c2){

    ## get pairs using C
    ## ensure that values of c1 and c2 are between 0 and n1
    res <- sortPairs(c1, c2)

    ## get ARI using pairs
    N <- length(c1)

    stot <- sum(choose(res$nij, 2), na.rm=TRUE)
    srow <- sum(choose(res$ni., 2), na.rm=TRUE)
    scol <- sum(choose(res$n.j, 2), na.rm=TRUE)

    expectedIndex <-(srow*scol)/(choose(N,2))
    maximumIndex <- (srow+scol)/2

    if (expectedIndex == maximumIndex & stot != 0) {
      return(1)
    } else {
      return((stot-expectedIndex)/(maximumIndex-expectedIndex))
    }
}

#' @title Rand Index
#' @description A function to compute the rand index between two classifications
#'
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @return a scalar with the rand index.
#' @seealso \code{\link{ARI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[,-5])), 4)
#' RI(cl,iris$Species)
#' @export
RI <- function(c1, c2){
  ## get pairs using C
  ## ensure that values of c1 and c2 are between 0 and n1
  res <- sortPairs(c1, c2)

  ## get ARI using pairs
  N <- length(c1)

  ## return the rand-index
  return(1 + (sum(res$nij^2) - (sum(res$ni.^2) + sum(res$n.j^2))/2)/choose(N,2))
}

#' @title Entropy
#' @description A function to compute the empirical entropy for two vectors of classification and the joint entropy
#'
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @return a list with the two conditional entropies and the joint entropy.
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[,-5])), 4)
#' entropy(cl,iris$Species)
#' @export
entropy <- function(c1, c2){
  res <- sortPairs(c1, c2)

  N <- length(c1)

  H.UV <- - sum(res$nij * log(res$nij))/N + log(N)
  H.U  <- - sum(res$ni. * log(res$ni.))/N + log(N)
  H.V  <- - sum(res$n.j * log(res$n.j))/N + log(N)

  return(list(UV = H.UV, U = H.U, V = H.V))
}

#' @title measures of similarity between two classification
#' @description A function various measures of similarity between two classifications
#'
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @return a list with the RI, ARI, NMI, NVI and NID.
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{ARI}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[,-5])), 4)
#' clustComp(cl,iris$Species)
#' @export
clustComp <- function(c1, c2) {

  H   <- entropy(c1,c2)

  MI  <- - H$UV + H$U + H$V
  VI  <- H$UV - MI
  NVI <- 1 - MI/H$UV
  ID  <- max(H$U, H$V) - MI
  NID <- 1 - MI / max(H$U, H$V)
  NMI <- MI / max(H$U, H$V)

  return(list(RI  = RI(c1,c2)             ,
              ARI = ARI(c1,c2)            ,
              MI  = - H$UV + H$U + H$V    ,
              VI  = H$UV - MI             ,
              NVI = 1 - MI/H$UV           ,
              ID  = max(H$U, H$V) - MI    ,
              NID = 1 - MI / max(H$U, H$V),
              NMI = MI / max(H$U, H$V)
  ))
}

#' @title Normalized mutual information (NMI)
#' @description A function to compute the NMI between two classifications
#'
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @return a scalar with the normalized mutual information .
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{ARI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[,-5])), 4)
#' NMI(cl,iris$Species)
#' @export
NMI <- function(c1, c2, variant=c("max", "min", "sqrt", "sum", "joint")) {

  variant <- match.arg(variant)

  H   <- entropy(c1,c2)

  MI  <- - H$UV + H$U + H$V

  D <- switch(variant,
          "max" = max(H$U, H$V),
          "sqrt" = sqrt(H$U * H$V),
          "min" = min(H$U, H$V),
          "sum" = .5*(H$U + H$V),
         "joint" = H$UV)
  return(MI / D)
}

#' @title Normalized information distance (NID)
#' @description A function to compute the NID between two classifications
#'
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @return a scalar with the normalized information distance .
#' @seealso \code{\link{RI}}, \code{\link{NMI}}, \code{\link{NVI}}, \code{\link{ARI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[,-5])), 4)
#' NID(cl,iris$Species)
#' @export
NID <- function(c1, c2) {

  H   <- entropy(c1,c2)

  MI  <- - H$UV + H$U + H$V
  return(1 - MI / max(H$U, H$V))
}

#' @title Normalized variation of information (NVI)
#' @description A function to compute the NVI between two classifications
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @return a scalar with the normalized variation of information.
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NMI}}, \code{\link{ARI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[,-5])), 4)
#' NVI(cl,iris$Species)
#' @export
NVI <- function(c1, c2) {

  H   <- entropy(c1,c2)
  MI  <- - H$UV + H$U + H$V
  return(1 - MI/H$UV)
}
