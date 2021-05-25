##' @title Get extent of list
##'
##' @description Given a list of SpatialPolygons, return an extent
##' object that encompasses all items. 
##'
##' @param shapes a list of SpatialPolygons or simple features
##' @param format if \code{raster}, then returns a raster extent object,
##' 	otherwise, returns a sf formatted bbox object.
##'
##' @return An object of class \code{extent} or \code{bbox}. 
##'
##' @author Pascal Title
##' 
##' @examples
##' getExtentOfList(tamiasPolyList)
##' 
##' @export


getExtentOfList <- function(shapes, format = 'raster') {
	
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

	if (format == 'raster') {
		res <- raster::extent(shapes[[1]])
		res@xmin <- minLong
		res@xmax <- maxLong
		res@ymin <- minLat
		res@ymax <- maxLat
	} else {
		res <- x[[1]]
		res[[1]]<- minLong
		res[[2]] <- minLat
		res[[3]] <- maxLong
		res[[4]] <- maxLat
	}
	
	return(res)
}	
