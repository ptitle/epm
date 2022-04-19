##' @title Convert raster to sf grid
##'
##' @description Convert a raster to sf polygons object, matching the attributes
##'  of the target object. 
##'
##' @param x rasterLayer or rasterStack or SpatRaster
##' @param target epmGrid or sf object
##' @param fun function for summarizing raster cells to polygons
##' @param crop if TRUE, the raster will be cropped to the bounding box of the target
##'
##' @details 
##' By default, raster cells that overlap with target grid cell polygons
##' 	will be averaged. If target is a raster grid, then \code{terra::resample}
##'		is used. 
##' 
##'
##' @return sf polygons object
##'
##' @author Pascal Title
##'
##' @examples
##' library(terra)
##' 
##' # We have a terra grid object (for example, climate data read in as a raster)
##' # Here, we are just generating some random data for demo
##' env <- rast(vect(tamiasEPM[[1]]), resolution = 100000)
##' env[] <- sample(1:100, ncell(env), replace = TRUE)
##' plot(env)
##'
##' # Now, if we are interested in doing analyses of environmental data in relation to 
##' # the epmGrid data we have, we want to convert the env data to the same grid structure
##' # where the cells align and where raster grid values are resampled and averaged.
##'
##' newgrid <- rasterToGrid(env, target = tamiasEPM, fun = 'mean')
##' plot(newgrid)
##' 
##' 
##' 
##' 
##' 
##' @export



rasterToGrid <- function(x, target, fun = 'mean', crop = TRUE) {
	
	if (!inherits(x, c('RasterLayer', 'rasterStack', 'SpatRaster'))) {
		stop('x must be either a rasterLayer, rasterStack or SpatRaster.')
	}
	if (!inherits(target, c('epmGrid', 'sf'))) {
		stop('target must be a epmGrid object, or a sf polygon object')
	}
	if (inherits(x, c('RasterLayer', 'rasterStack'))) {
		x <- as(x, 'SpatRaster')
	}
	
	if (inherits(target, 'epmGrid')) {
		target <- target[[1]]
	}

	if (inherits(target, 'sf')) {
	
		## terra::cells(x, as(target, 'SpatVector'))
		
		if (crop) {
			e <- sf::st_bbox(target)
			e <- e[c('xmin', 'xmax', 'ymin', 'ymax')]
			e <- terra::ext(e)
			x <- terra::crop(x, e)
		}
			
		# convert raster to sf grid
		# if (names(x) == '') names(x) <- 'dat'
		# d <- terra::as.data.frame(terra::as.polygons(x, dissolve = FALSE), geom = "hex")
		# d$geometry <- structure(as.list(d$geometry), class = "WKB")
		# x2 <- sf::st_as_sf(d, crs = x@ptr$get_crs("wkt"))
		
		# # aggregate x by target's polygons, summarizing with supplied function	
		# res <- aggregate(x2, by = target, fun)
		
		df <- terra::extract(x, terra::vect(target), list = FALSE)
		newvals <- aggregate(df, by = list(df$ID), FUN = fun)
		res <- sf::st_sf(vals = newvals[,3], geometry = sf::st_geometry(target))
		
	} else if (inherits(target, 'SpatRaster')) {
		
		res <- terra::resample(x, target)
		
	} else {
		stop('target format not recognized.')
	}
	
	return(res)
}





