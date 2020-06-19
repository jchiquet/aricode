#' @title Sort Pairs
#'
#' @description A function to sort pairs of integers or factors and identify the pairs
#' @param c1 a vector of length n with value between 0 and N1 < n
#' @param c2 a vector of length n with value between 0 and N2 < n
#' @param spMat logical: send back the contingency table as sparsely encoded (cost more than the algorithm itself). Default is FALSE
#' @import Matrix
#' @export
sortPairs <- function(c1, c2, spMat=FALSE){
  if (anyNA(c1) | anyNA(c2))
    stop("NA are not supported.")

  if (((!is.vector(c1) & !is.factor(c1)) | is.list(c1)) | ((!is.vector(c2) & !is.factor(c2)) | is.list(c2)))
    stop("c1 and c2 must be vectors or factors but not lists.")

  if (length(c1) != length(c2))
    stop("the two vectors must have the same length.")

  n <- length(c1)

  ## if c1 and c2 are integer
  if (is.integer(c1) & is.integer(c2)) {
    mylevels <- list(c1 = unique(c1), c2 = unique(c2))
    #c1 <- c1 - min(c1)
    #c2 <- c2 - min(c2)
    ## if the range is not adapted to the C code
    #if (!(max(c1) <= n-1 & max(c2) <= n-1)) { 
    ## HERE WE MAKE ENSURE THAT ALL INDEX FROM 1 to K and 1 to L are present which is 
    ## usefull for sparseMatrix and the calculation of some criteria linking: n_k. n_.l to n_kl
    ## TODO add a skip parameter if we can ensure that this is not needed?
      c1 <- as.integer(factor(c1, levels = mylevels$c1)) - 1L
      c2 <- as.integer(factor(c2, levels = mylevels$c2)) - 1L
    #} 
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


  i_order <- order(c1, c2, method="radix") - 1L
  out <- countPairs(c1, c2, i_order)

  if (spMat) {
    spOut <- sparseMatrix(i=out$pair_c1,
                          j=out$pair_c2,
                          x=out$pair_nb,
                          dims=sapply(mylevels,length),
                          dimnames = mylevels, index1=FALSE)
  } else {
    spOut <- NULL
  }

  res <- list(spMat = spOut,
              levels = mylevels,
              nij = out$pair_nb,
              ni. = out$c1_nb,
              n.j = out$c2_nb,
              pair_c1 = out$pair_c1,
              pair_c2 = out$pair_c2
  )
  res
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
      res <- 1
    } else {
      res <- (stot-expectedIndex)/(maximumIndex-expectedIndex)
    }
    res
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
  res <- 1 + (sum(res$nij^2) - (sum(res$ni.^2) + sum(res$n.j^2))/2)/choose(N,2)
  res
}

#' @title Modified Adjusted Rand Index
#' @description A function to compute a modified adjusted rand index between two classifications as proposed by Sundqvist et al. in prep, based on a multinomial model. 
#'
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @return a scalar with the modified ARI.
#' @seealso \code{\link{ARI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[,-5])), 4)
#' MARI(cl,iris$Species)
#' @export
MARI <- function(c1, c2){
  ## get pairs using C
  ## ensure that values of c1 and c2 are between 0 and n1
  res <- sortPairs(c1, c2)
  N <- length(c1)
  ##

  stot <- sum(choose(res$nij, 2), na.rm=TRUE)
  srow <- sum(choose(res$ni., 2), na.rm=TRUE)
  scol <- sum(choose(res$n.j, 2), na.rm=TRUE)

  ## using Lemma 3.3
  ## triplets
  T1 <- 2*N
  T2 <- sum(res$nij * res$ni.[res$pair_c1+1] * res$n.j[res$pair_c2+1], na.rm=TRUE)
  T3 <- -sum(res$nij^2, na.rm=TRUE) - sum(res$ni.^2, na.rm=TRUE) - sum(res$n.j^2, na.rm=TRUE)
  
  ## quadruplets (and division by 6 choose(N, 4)
  expectedIndex <- (srow*scol - stot - (T1+T2+T3)) / (6 *choose(N, 4))

  ## return the rand-index
  expectedIndex <- expectedIndex * choose(N, 2) ## RESCALE SO THAT THE CODE IS EQUIVALENT TO THE ARI
  maximumIndex <- (srow+scol)/2 
  if (expectedIndex == maximumIndex & stot != 0) {
    res <- 1
  } else {
    res <- (stot-expectedIndex)/(maximumIndex-expectedIndex)
  }
    res
  ## return the adjusted (and divided) rand-index
  res
}

#' @title raw Modified Adjusted Rand Index
#' @description A function to compute a modified adjusted rand index between two classifications as proposed by Sundqvist et al. in prep, based on a multinomial model. Raw means, that the index is not divided by the (maximum - expected) value.
#'
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @return a scalar with the modified ARI.
#' @seealso \code{\link{ARI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[,-5])), 4)
#' MARI.raw(cl,iris$Species)
#' @export
MARI.raw <- function(c1, c2){
  ## get pairs using C
  ## ensure that values of c1 and c2 are between 0 and n1
  res <- sortPairs(c1, c2)
  N <- length(c1)
  ##

  stot <- sum(choose(res$nij, 2), na.rm=TRUE)
  srow <- sum(choose(res$ni., 2), na.rm=TRUE)
  scol <- sum(choose(res$n.j, 2), na.rm=TRUE)

  ## using Lemma 3.3
  ## triplets
  T1 <- 2*N
  T2 <- sum(res$nij * res$ni.[res$pair_c1+1] * res$n.j[res$pair_c2+1], na.rm=TRUE)
  T3 <- -sum(res$nij^2, na.rm=TRUE) - sum(res$ni.^2, na.rm=TRUE) - sum(res$n.j^2, na.rm=TRUE)
  
  ## quadruplets (and division by 6 choose(N, 4)
  expectedIndex <- (srow*scol - stot - (T1+T2+T3)) / (6 *choose(N, 4))

  ## return the rand-index
  res <- (stot / choose(N, 2)) - expectedIndex
  res
}

#' @title Chi-square statistics
#' @description A function to compute the Chi-2 statistics
#'
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @return a scalar with the chi-square statistics.
#' @seealso \code{\link{ARI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[,-5])), 4)
#' CHI2(cl,iris$Species)
#' @export
CHI2 <- function(c1, c2){
  ## get pairs using C
  ## ensure that values of c1 and c2 are between 0 and n1
  res <- sortPairs(c1, c2)
  N <- length(c1)
  
  res <- N* sum(res$nij^2 / (res$ni.[res$pair_c1+1] * res$n.j[res$pair_c2+1]) )
  res <- res - N
  res
}


#' @title Entropy
#' @description A function to compute the empirical entropy for two vectors of classification and the joint entropy
#'
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @return a list with the two conditional entropies, the joint entropy and output of sortPairs.
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

  res <- list(UV = H.UV, U = H.U, V = H.V, sortPairs = res)
  res
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
  
  EMI <- expected_MI(as.integer(H$ni.), as.integer(H$n.j))


  res <- list(RI  = RI(c1,c2)             ,
              ARI = ARI(c1,c2)            ,
              MI  = - H$UV + H$U + H$V    ,
              AMI = (- H$UV + H$U + H$V - EMI) / (max(H$U,H$V) - EMI),
              VI  = H$UV - MI             ,
              NVI = 1 - MI/H$UV           ,
              ID  = max(H$U, H$V) - MI    ,
              NID = 1 - MI / max(H$U, H$V),
              NMI = MI / max(H$U, H$V),
	      ## new
	      CHI2 = CHI2(c1,c2),
  	      MARI = MARI(c1,c2),
	      MARI.raw = MARI.raw(c1,c2)
  )
  res
}

#' @title Adjusted Mutual Information
#' @description A function to compute the adjusted mutual information between two classifications
#'
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @return a scalar with the adjusted rand index.
#' @seealso \code{\link{ARI}}, \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[,-5])), 4)
#' AMI(cl,iris$Species)
#' @export
AMI <- function(c1, c2){

  H   <- entropy(c1,c2)
  MI  <- - H$UV + H$U + H$V
  EMI <- expected_MI(as.integer(H$ni.), as.integer(H$n.j))

  res <- (MI - EMI) / (max(H$U,H$V) - EMI)
  res
}

#' @title Normalized mutual information (NMI)
#' @description A function to compute the NMI between two classifications
#'
#' @param c1 a vector containing the labels of the first classification. Must be a vector of characters, integers, numerics, or a factor, but not a list.
#' @param c2 a vector containing the labels of the second classification.
#' @param variant a string in ("max", "min", "sqrt", "sum", "joint"): different variants of NMI. Default use "max".
#' @return a scalar with the normalized mutual information .
#' @seealso \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{ARI}}, \code{\link{clustComp}}
#' @examples
#' data(iris)
#' cl <- cutree(hclust(dist(iris[,-5])), 4)
#' NMI(cl,iris$Species)
#' @export
NMI <- function(c1, c2, variant = c("max", "min", "sqrt", "sum", "joint")) {

  variant <- match.arg(variant)

  H   <- entropy(c1,c2)

  MI  <- - H$UV + H$U + H$V

  D <- switch(variant,
          "max" = max(H$U, H$V),
          "sqrt" = sqrt(H$U * H$V),
          "min" = min(H$U, H$V),
          "sum" = .5*(H$U + H$V),
         "joint" = H$UV)
  res <- MI / D
  res
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
  res <- 1 - MI / max(H$U, H$V)
  res
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
  res <- 1 - MI/H$UV
  res
}




