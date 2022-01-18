##' @title Partial Disparity
##'
##' @description Calculate species-specific partial disparity, relative to some
##'   group mean.
##'
##' @param dat matrix of multivariate morphological data
##' @param groupMean if \code{NULL}, calculated from \code{dat}, otherwise can
##'   be provided as a vector of mean values
##'
##' @details Calculates partial disparities, as in Foote 1993. By default, the
##'   group mean is calculated from the full input data.
##'
##'
##' @return numeric vector
##'
##' @author Pascal Title
##'
##' @examples
##' tamiasTraits[1:5, 1:5]
##' getSpPartialDisparities(tamiasTraits)
##'
##'
##' @export



# if groupMean is NULL, then calculated from dat, otherwise can be a vector of length ncol(dat)
getSpPartialDisparities <- function(dat, groupMean = NULL) {
	
	# calculate partial disparity as the squared distance of each species from the group mean, 
	# weighted by 1 / N - 1, where N = overall group size
	
	if (is.null(groupMean)) {
		groupMean <- colMeans(dat)
	} else {
		if (length(groupMean) != ncol(dat)) {
			stop('groupMean must have same length as ncol dat.')
		}
	}
	
	res <- utils::tail(as.matrix(dist(rbind(dat, groupMean)))[, 1:nrow(dat)], 1) ^ 2 * (1 / (nrow(dat) - 1))
	res <- setNames(as.vector(res), colnames(res))
	return(res)
}
