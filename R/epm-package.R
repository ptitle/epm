##' @title EcoPhyloMapper (epm)
##' @description An R package that facilitates the aggregation of species'
##' geographic ranges from vector or raster spatial data, and that enables
##' the calculation of various morphological and phylogenetic metacommunity metrics
##' across geography.
##'
##' A detailed wiki for the R package can be found on the epm github page:
##' \url{https://github.com/ptitle/epm/wiki#table-of-contents}
##' 
##' @author Pascal O. Title, Donald L. Swiderski, Miriam L. Zelditch
##' 
##' @name epm
##' @docType package
##' 
##' 
##' @useDynLib epm
##' @importFrom Rcpp evalCpp
##' @importFrom Rcpp sourceCpp
##'
##' @import methods
##' 
##' 
##' @importFrom graphics axis grconvertX grconvertY par identify rect
##' @importFrom stats cophenetic cov dist prcomp runif sd setNames
##' @importFrom utils setTxtProgressBar txtProgressBar
##' @importFrom grDevices gray
##' @importFrom stats aggregate
##' 
##' 
##' 
##' 
NULL