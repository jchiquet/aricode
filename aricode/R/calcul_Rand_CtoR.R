#' @title Sort Pairs
#' @description A function to sort pairs of integers or factors
#' @import Matrix
#' @useDynLib aricode
#' @export
sortPairs <- function(
  ### Fonction to sort each individual and then identifiy paris
  c1,
  ### a vector of integer of length n with value between 0 and N1 < n
  c2
  ### a vector of integer of length n with value between 0 and N2 < n
  ### if N1 or N2 > n should lead to segfault
  ){

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
        ## if neither a factor nor an integer, force to factor then integer
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
                 nzero      = as.integer(1), package="aricode")

    return(list(spMat = sparseMatrix(i=result$pair_c1[1:result$nzero],
                          j=result$pair_c2[1:result$nzero],
                          x=result$pair_count[1:result$nzero],
                          dims=sapply(mylevels,length),
                          dimnames = mylevels, index1=FALSE),
                pair_count = result$pair_count,
                cla1_count = result$count1,
                cla2_count = result$count2
    ))
}

#' @title Adjusted Rand Index
#' @description A function to compute the adjusted rand index
#'
#' @export
ARI <- function(c1, c2){

    ##TODO: handle bordeline and trivial cases

    ## get pairs using C
    ## ensure that values of c1 and c2 are between 0 and n1
    result <- sortPairs(c1, c2)

    ## get ARI using pairs
    N <- length(c1)

    nij <- result$pair_count[which(result$pair_count > 1)]
    ni. <- result$cla1_count[which(result$cla1_count > 1)]
    n.j <- result$cla2_count[which(result$cla2_count > 1)]

    stot <- sum(choose(nij, 2), na.rm=TRUE)
    srow <- sum(choose(ni., 2), na.rm=TRUE)
    scol <- sum(choose(n.j, 2), na.rm=TRUE)

    expectedIndex <-(srow*scol)/(choose(N,2))
    maximumIndex <- (srow+scol)/2

    if (expectedIndex == maximumIndex & stot != 0) {
      return(1)
    } else {
      return((stot-expectedIndex)/(maximumIndex-expectedIndex))
    }
}

#' @title Rand Index
#' @description A function to compute the rand index
#'
#' @export
RI <- function(c1, c2){
  ## get pairs using C
  ## ensure that values of c1 and c2 are between 0 and n1
  result <- sortPairs(c1, c2)

  ## get ARI using pairs
  N <- length(c1)

  nij <- result$pair_count[which(result$pair_count > 1)]
  ni. <- result$cla1_count[which(result$cla1_count > 1)]
  n.j <- result$cla2_count[which(result$cla2_count > 1)]

  ## return the rand-index
  return(1 + (sum(nij^2) - (sum(ni.^2) + sum(n.j^2))/2)/choose(N,2))
}

#' @title entropy
#' @description A function to compute class 1, 2 and pair entropy
#'
#' @export
entropy <- function(c1, c2){
  result <- sortPairs(c1, c2)

  N <- length(c1)

  nij <- result$pair_count[which(result$pair_count > 1)]
  ni. <- result$cla1_count[which(result$cla1_count > 1)]
  n.j <- result$cla2_count[which(result$cla2_count > 1)]

  H.UV <- - sum(nij * log(nij))/N + log(N)
  H.U  <- - sum(ni. * log(ni.))/N + log(N)
  H.V  <- - sum(n.j * log(n.j))/N + log(N)

  return(list(UV = H.UV, U = H.U, V = H.V))
}

#' @title measures of similarity between classification
#' @description A function various measures of similarity between two classifications
#'
#' @export
clustCompM <- function(c1, c2) {

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
#' @description A function various measures of similarity between two classifications
#'
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
#' @description A function various measures of similarity between two classifications
#'
#' @export
NID <- function(c1, c2) {

  H   <- entropy(c1,c2)

  MI  <- - H$UV + H$U + H$V
  return(1 - MI / max(H$U, H$V))
}

#' @title Normalized variation of information (NVI)
#' @description A function various measures of similarity between two classifications
#'
#' @export
NVI <- function(c1, c2) {

  H   <- entropy(c1,c2)

  MI  <- - H$UV + H$U + H$V
  return(1 - MI/H$UV)
}
