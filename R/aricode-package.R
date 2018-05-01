#' aricode
#'
#' A package for efficient computations of standard clustering comparison measures. Available measures are described in the paper of Vinh et al, JMLR, 2009 (see reference below).
#'
#' Traditional implementations (e.g., function \code{adjustedRandIndex} of package \code{mclust}) are in Omega(n + u v) where n is the size of the vectors the classifications of which are to be compared, u and v are the respective number of classes in each vectors. Here, the implementation is in Theta(n), plus the gain of speed due to the C code.
#'
#' @section Functions in aricode:
#' The functions included in aricode are:
#' \itemize{
#'  \item{ARI:}{ computes the adjusted rand index}
#'  \item{RI:}{ computes the rand index}
#'  \item{NVI:}{ computes the normalized variation information}
#'  \item{NID:}{ computes the normalized information distance}
#'  \item{NMI:}{ computes the normalized mutual information}
#'  \item{entropy:}{ computes the conditional and joint entropies}
#'  \item{clustComp:}{ computes all clustering comparison measures at once}
#' }
#' @author Julien Chiquet \email{julien.chiquet@@gmail.com}
#' @author Guillem Rigaill \email{guillem.rigaill@@evry.inra.fr}
#' @references Nguyen Xuan Vinh, Julien Epps, and James Bailey. "Information theoretic measures for clusterings comparison: Variants, properties, normalization and correction for chance." Journal of Machine Learning Research 11.Oct (2010): 2837-2854. as described in Vinh et al (2009)
#' @seealso \code{\link{ARI}}, \code{\link{RI}}, \code{\link{NID}}, \code{\link{NVI}}, \code{\link{NMI}}, \code{\link{entropy}}, \code{\link{clustComp}}
#' @docType package
#' @name aricode
#' @useDynLib aricode
#' @importFrom Rcpp sourceCpp
NULL
