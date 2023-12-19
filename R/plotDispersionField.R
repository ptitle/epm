##' @title Plot dispersion fields
##'
##' @description For a set of specified coordinates, plot a richness map for 
##' the species that are found at those coordinates. 
##'
##' @param x object of class \code{epmGrid}
##' @param coords coordinates as either a spatial points object (sp or sf), 
##' 	a matrix/dataframe with two columns or a numeric vector of c(long, lat).
##' @param plotCoords boolean; should the coordinates be plotted as well?
##' @param legend boolean; should legend be included?
##' @param col either a vector of color names that will be interpolated, 
##' or a color ramp function that takes an integer 
##' (see for example \code{\link{colorRampPalette}}).
##'	@param basemap if \code{'none'}, then only the grid is plotted. 
##'		If \code{'worldmap'}, then vector map is plotted.
##' 	If \code{'interactive'}, then the \code{mapview} package is used.
##' @param lwd grid cell border width
##' @param borderCol color for grid cell borders
##' @param alpha opacity of all colors and borders, ranging from 0 
##' (fully transparent) to 1 (fully opaque)
##' @param includeFrame boolean; include frame around plot?
##' @param use_tmap boolean; if FALSE, plotting will be done via sf instead 
##' of tmap package
##' @param add logical, add to existing plot?
##'
##' @return Nothing is returned.
##'
##' @details Assemblage dispersion fields represent an overlapping of 
##' geographic ranges for the taxa that occur in the focal grid cells. 
##' 
##' @author Pascal Title
##'
##' @references
##' Graves, G. R., & Rahbek, C. (2005). Source pool geometry and the assembly 
##' of continental avifaunas. Proceedings of the National Academy of Sciences, 
##' 102(22), 7871–7876.
##' 
##' @examples
##' # plotDispersionField(tamiasEPM, c(-1944951, 69588.74))
##' plotDispersionField(tamiasEPM, c(-1944951, 69588.74), use_tmap = FALSE)
##' 
##' @export


plotDispersionField <- function(x, coords, plotCoords = TRUE, legend = TRUE, col, lwd = 0.5, basemap = 'worldmap', borderCol = 'black', alpha = 1, includeFrame = FALSE, use_tmap = TRUE, add = FALSE) {
	
	# x = tamiasEPM; coords = c(-1413764.0, 573610.8); legend = TRUE; basemap = 'worldmap'; lwd = 0.15; borderCol = 'black'; alpha = 1; includeFrame = FALSE; use_tmap = TRUE; add = FALSE; plotCoords = TRUE

	if (!inherits(x, 'epmGrid')) {
		stop('Object must be of class epmGrid')
	}
	
	if (!basemap %in% c('worldmap', 'interactive', 'none')) {
		stop("basemap argument must be 'none', 'worldmap' or 'interactive'.")
	}
	
	if (use_tmap & !requireNamespace('tmap', quietly = TRUE)) {
		message('\ttmap package not installed -- defaulting to sf/terra plot')
		use_tmap <- FALSE
	}
	
	if (add) {
		use_tmap <- FALSE
		basemap <- 'none'
		legend <- FALSE
	}
	
	siteTaxa <- extractFromEpmGrid(x, spatial = coords)
	
	# unroll the compact species cell list
	gridcells <- expandSpeciesCellList(x)
	
	# intersect it with our taxa of interest
	gridcells <- intersectList(gridcells, siteTaxa)
	richvec <- lengths(gridcells)
	richvec[sapply(gridcells, anyNA)] <- 0
	
	if (inherits(x[[1]], 'sf')) {
		grid <- x[[1]][richvec > 0, ]
		grid$dispersionCount <- richvec[richvec > 0]
	} else {
		grid <- terra::rast(x[[1]], nlyrs = 1)
		grid[richvec > 0] <- richvec[richvec > 0]
	}
	
	ncol <- 1000
	isInt <- FALSE
	if (max(richvec) <= 10) {
		ncol <- max(richvec)
		isInt <- TRUE
	}

	if (missing(col)) {
		# default color ramp: viridis turbo, but without the most extreme values
		colramp <- function(n) viridisLite::turbo(n = n, begin = 0.1, end = 0.9)
	} else {
		if (inherits(col, 'function')) {
			colramp <- col
		} else if (inherits(col, 'character') & length(col) > 1) {
				colramp <- grDevices::colorRampPalette(col)
		} else {
			stop('color not understood.')
		}
	}

	colors <- colramp(ncol)
	colors <- grDevices::adjustcolor(colors, alpha.f = alpha)
	borderColor <- grDevices::adjustcolor('black', alpha.f = alpha)
	
	if (!inherits(coords, 'sf')) {
		if (is.vector(coords)) {
			coords <- matrix(coords, nrow = 1, ncol = 2)
		}
		if (is.matrix(coords)) {
			coords <- as.data.frame(coords)
		}
		coords <- sf::st_as_sf(coords, coords = 1:2, crs = sf::st_crs(x[[1]]))
	}


	if (use_tmap) {
		
		if (inherits(x[[1]], 'sf')) {

		    if (basemap != 'interactive' & getOption("tmap.mode") == 'view') {
		        tmap::tmap_mode('plot')
		    }
		
			if (basemap == 'worldmap') {

				# wrld <- sf::st_transform(worldmap, crs = sf::st_crs(x[[1]]))
				map <- tmap::tm_shape(worldmap, is.master = FALSE, projection = sf::st_crs(x[[1]]), bbox = x[[1]]) + tmap::tm_borders(lwd = 0.5)
			}
			
			if (basemap == 'none') {
				map <- NULL
			}

			map <- map + tmap::tm_shape(grid) + tmap::tm_fill('dispersionCount', palette = colramp(ncol), legend.show = legend, title = '', style = 'pretty', as.count = TRUE, midpoint = NA, alpha = alpha) + tmap::tm_borders(col = borderCol, lwd = lwd, alpha = alpha) + tmap::tm_layout(frame = includeFrame, legend.outside = TRUE)
			
			# add the focal coordinate(s)
			if (plotCoords) {
				map <- map + tmap::tm_shape(coords) + tmap::tm_dots(col = 'white', size = 0.05, shape = 3)
			}
			
		} else {
			# tmap with terra grid

			# get bounding box of grid cells with data. We will do this so that we can leave out of the plot those grid cells that are empty and that may take up quite a bit of geographic space. Only an issue for rasters.
			datBB <- sf::st_bbox(sf::st_as_sf(as.data.frame(terra::xyFromCell(grid, which(!is.na(terra::values(grid))))), coords = 1:2, crs = sf::st_crs(x[[1]])))
					    
		    if (basemap != 'interactive' & getOption("tmap.mode") == 'view') {
		        tmap::tmap_mode('plot')
		    }
					
			if (basemap == 'worldmap') {

				# wrld <- sf::st_transform(worldmap, crs = sf::st_crs(x[[1]]))
				map <- tmap::tm_shape(worldmap, is.master = FALSE, projection = sf::st_crs(x[[1]]), bbox = datBB) + tmap::tm_borders(lwd = 0.5)
			}
			
			if (basemap == 'none') {
				map <- NULL
			}
			
			if (basemap == 'interactive') {
				map <- NULL
				tmap::tmap_mode('view')
			}

			map <- map + tmap::tm_shape(grid, bbox = datBB) + tmap::tm_raster(palette = colramp(ncol), legend.show = legend, title = '', style = 'pretty', as.count = TRUE, midpoint = NA, alpha = alpha) + tmap::tm_layout(frame = includeFrame, legend.outside = TRUE)
			
			# add the focal coordinate(s)
			if (plotCoords) {
				map <- map + tmap::tm_shape(coords) + tmap::tm_dots(col = 'white', size = 0.05, shape = 3)
			}
		}
		
		map
		
	} else {
		
		# not tmap
		if (inherits(x[[1]], 'sf')) {
		
			cols <- colramp(max(richvec))[grid$dispersionCount]
			
			plot(sf::st_geometry(grid), axes = includeFrame, main = NULL, key.pos = NULL, lwd = lwd, reset = FALSE, col = cols, border = borderColor, add = add)
		
		} else {
			
			# use terra plotting
			
			datBB2 <- terra::ext(grid)
			datBB2[1] <- datBB['xmin']
			datBB2[2] <- datBB['xmax']
			datBB2[3] <- datBB['ymin']
			datBB2[4] <- datBB['ymax']

			terra::plot(terra::crop(grid, datBB2), col = colramp(ncol), axes = includeFrame, legend = FALSE, plg = list(shrink = 0.7, title = ''), alpha = alpha, add = add)
			
		}
		
		# add the focal coordinate(s)
		if (plotCoords) {
			plot(coords, add = TRUE, col = 'white', cex = 0.5, pch = 3)
		}
		
		if (legend) {
			if (isInt) {
				nTicks <- max(richvec)
			} else {
				nTicks <- 3
			}
			addLegend(grid['dispersionCount'], location = 'topright', ramp = colramp, isInteger = isInt, ncolors = ncol, nTicks = nTicks)
		}
		
		if (basemap == 'worldmap') {
			# add map for context
			wrld <- sf::st_transform(worldmap, crs = sf::st_crs(x[[1]]))
			grXY <- graphics::par("usr")
			graphics::clip(grXY[1], grXY[2], grXY[3], grXY[4]) # this ensures that world map is constrained to plot region
			graphics::plot(wrld, add = TRUE, lwd = lwd)
		}
	}
}




	




