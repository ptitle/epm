##' @title Get extent of list
##'
##' @description Given a list of SpatialPolygons, return an extent
##' object that encompasses all items. 
##'
##' @param shapes a list of SpatialPolygons or simple features
##'
##' @return An object of class \code{bbox}. 
##'
##' @author Pascal Title
##' 
##' @examples
##' getExtentOfList(tamiasPolyList)
##' 
##' @export


getExtentOfList <- function(shapes) {
	
	if (inherits(shapes, c('sf', 'sfc')) & !inherits(shapes, 'list')) {
		shapes <- list(shapes)
	}

	if (inherits(shapes[[1]], c('SpatialPolygons', 'SpatialPolygonsDataFrame'))) {	
		shapes <- lapply(shapes, function(x) sf::st_as_sf(x))
	} 
	
	if (inherits(shapes[[1]], c('sf', 'sfc'))) {
		x <- lapply(shapes, function(x) sf::st_bbox(x))
	} else {
		stop('shapes object not recognized.')
	}

	minLong <- min(sapply(x, function(x) x$xmin, simplify = TRUE))
	maxLong <- max(sapply(x, function(x) x$xmax, simplify = TRUE))
	minLat <- min(sapply(x, function(x) x$ymin, simplify = TRUE))
	maxLat <- max(sapply(x, function(x) x$ymax, simplify = TRUE))

	res <- x[[1]]
	res[[1]]<- minLong
	res[[2]] <- minLat
	res[[3]] <- maxLong
	res[[4]] <- maxLat
	
	return(res)
}	
