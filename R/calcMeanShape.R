##' @title Calculate mean shape per grid cell
##'
##' @description For an epmGrid object that contains geometric
##' morphometric shape coordinates, calculate the per-grid-cell
##' mean shape. 
##'
##' @param x object of class \code{epmGrid}
##'
##' @details This function will ignore cells that are empty.
##'
##' @return a list with 2 elements:
##' 	(1) matrix where nrow = number of grid cells and ncol = 
##' the number of data columns. Each row is a vector of mean shape coordinates.
##' 	(2) a matrix of xy coordinates corresponding to those grid cells.
##'
##' @author Pascal Title
##' 
##' @examples
##' \donttest{
##' tamiasEPM
##' tamiasEPM <- addTraits(tamiasEPM, tamiasTraits)
##'	meanshape <- calcMeanShape(tamiasEPM)
##' 
##' head(meanshape[[1]])
##' head(meanshape[[2]])
##' }
##' @export



calcMeanShape <- function(x) {
	
	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}

	if (!inherits(x[['data']], c('matrix', 'data.frame'))) {
		stop('This function expects multivariate data.')
	}
	
	# expand to full cell species list
	cellComms <- expandSpeciesCellList(x)
	
	cellComms <- intersectList(cellComms, rownames(x[['data']]))
	ind <- which(sapply(cellComms, anyNA) == FALSE)
	cellComms <- cellComms[ind]
	
	meanShapeMat <- pbapply::pblapply(cellComms, function(y) colMeans(x[['data']][y, ]))
	meanShapeMat <- do.call(rbind, meanShapeMat)
	
	# fetch xy coordinates as first and second columns
	if (inherits(x[[1]], 'SpatRaster')) {
	    coordmat <- terra::xyFromCell(x[[1]], ind)
	} else if (inherits(x[[1]], 'sf')) {
	    coordmat <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(x[[1]][ind,])))
	}
	colnames(coordmat) <- c('x', 'y')
		
	return(list(meanShapeMat, coordmat))
}