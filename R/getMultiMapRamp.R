##' @title Extract min and max for multiple epmGrids
##'
##' @description Extracts the range of values across a list of input
##' 	objects for use in plotting
##'
##' @param x list of objects of class \code{epmGrid}, \code{RasterLayer}
##' 	\code{SpatRaster} or \code{sf} objects.
##'
##' @details If the user would like to plot multiple epmGrid objects
##' with a standardized color ramp, then the returned values from this function
##' can be supplied to \code{\link{plot.epmGrid}}. Also works with RasterLayer
##' and sf objects. For sf object, only one attribute can be specified.  
##'
##' @return a numeric vector of length 2: overall min and max value. 
##'
##' @author Pascal Title
##' 
##' @examples
##' library(raster)
##' tamiasEPM
##'
##' # create a dummy raster for demonstration purposes.
##' ras <- raster()
##' values(ras) <- runif(ncell(ras), min = 0, max = 40)
##' 
##' getMultiMapRamp(list(tamiasEPM, ras))
##'  
##' @export


getMultiMapRamp <- function(x) {
	
	if (!inherits(x, 'list')) {
		x <- list(x)
	}
	
	if (!all(unique(sapply(x, class)) %in% c('RasterLayer', 'SpatRaster', 'epmGrid', 'sf', 'SpatialPointsDataFrame', 'SpatialPolygonsDataFrame'))) {
		stop('Input cannot be handled.')
	}
	
	rangeList <- vector('list', length(x))
	for (i in 1:length(x)) {
		
		if (inherits(x[[i]], 'epmGrid')) {
			tmp <- x[[i]][[1]]
			datCol <- attributes(x[[i]])$metric
			if (length(datCol) > 1) {
				stop(paste('List item', i, 'contains more than one attribute.'))
			}
			if (inherits(tmp, 'sf')) {
				rangeList[[i]] <- range(tmp[[datCol]], na.rm = TRUE)
			} else {
				rangeList[[i]] <- range(terra::minmax(tmp[[datCol]]))
			}
		}
		
		if (inherits(x[[i]], 'sf')) {
			tmp <- x[[i]]
			datCol <- setdiff(colnames(tmp), attributes(tmp)$sf_column)
			if (length(datCol) > 1) {
				stop(paste('List item', i, 'contains more than one attribute.'))
			}
			rangeList[[i]] <- range(tmp[[datCol]], na.rm = TRUE)
		}
		
		if (inherits(x[[i]], 'RasterLayer')) {
			x[[i]] <- as(x[[i]], 'SpatRaster')
		}
		
		if (inherits(x[[i]], 'SpatRaster')) {
			rangeList[[i]] <- range(terra::minmax(x[[i]]))
		}
		
		if (inherits(x[[i]], c('SpatialPointsDataFrame', 'SpatialPolygonsDataFrame'))) {
			if (ncol(tmp@data) > 1) {
				stop(paste('List item', i, 'contains more than one attribute.'))
			}
			rangeList[[i]] <- range(x[[i]]@data[,1], na.rm = TRUE)
		}
	}
		
	# find global min and max
	return(c(sapply(rangeList, min), sapply(rangeList, max)))
		
}
