##' @title Interactively choose extent
##'
##' @description Given a list of polygons or point occurrences, sets up
##' 	an interactive plot	to allow the user to draw the desired extent. 
##' This can be used to	define the extent in \code{\link{createEPMgrid}}.
##'
##' @param polyList a list of Simple Feature polygons or points.
##'
##' @param cellType either \code{hexagon} or \code{square}.
##' 
##' @param bb c(xmin, xmax, ymin, ymax) to limit the extent for 
##' the interactive plot.
##'
##' @details This function returns both a sf polygon and the same polygon
##' 		as a WKT string. Either can be supplied to \code{\link{createEPMgrid}}
##' 		as the extent. A recommended strategy is to use this function to find
##' 		your extent, and to copy/paste the WKT string into your R script so that 
##' 		you can retain it for future use, and maintain reproducibility. 
##' 		See example.
##'
##'		What is chosen for \code{cellType} has no effect on what you might choose in 
##'		\code{\link{createEPMgrid}}. Square cells will probably be fastest. If hexagons
##'		are selected, grid cell points are plotted instead of polygons to speed up plotting.
##'
##'		You may see the message \code{Failed to compute min/max, no valid pixels found in 
##' 		sampling. (GDAL error 1)} . This just means that a species did not register in any 
##' 		grid cells. This can be ignored.
##'
##' 	The basemap is from \url{https://www.naturalearthdata.com/}. 
##'
##' @return A list with a polygon, and its WKT string 
##'
##' @author Pascal Title
##' 
##' @examples
##' if (interactive()) {
##' 	ex <- interactiveExtent(tamiasPolyList)
##'
##' 	# You can use this as the extent in createEPMgrid
##' 	grid <- createEPMgrid(tamiasPolyList, resolution = 50000, extent = ex$wkt)
##'
##' 	# One way to make your code reproducible would be to copy/paste the wkt 
##' 	# in your code for future use:
##' 	ex <- interactiveExtent(tamiasPolyList)
##' 	ex$wkt
##' 	customExtent <- "POLYGON ((-2238201 3532133, -2675450 1722657, -2470677 -317634, 
##' 	-1863632 -1854074, -521614.8 -2170280, -349356.8 799040.9, -2238201 3532133))"
##'
##' 	grid <- createEPMgrid(tamiasPolyList, resolution = 50000, extent = customExtent)
##'
##'}
##' @export




interactiveExtent <- function(polyList, cellType = 'square', bb = NULL) {
	
	if (!inherits(polyList, 'Spatial') & inherits(polyList[[1]], 'Spatial')) {
		# if class SpatialPolygons, convert to sf
		for (i in 1:length(polyList)) {
			polyList[[i]] <- sf::st_as_sf(polyList[[i]])
		}
	}
	
	if (inherits(polyList, 'Spatial')) {
		polyList <- sf::st_as_sf(polyList)
	}
	
	if (!inherits(polyList, 'list')) {
		if (grepl('POINT', unique(as.character(sf::st_geometry_type(polyList))))) {
			polyCRS <- sf::st_crs(polyList)
			polyList <- occurrenceFormatting(polyList)
			polyList <- lapply(polyList, function(x) sf::st_as_sf(x, coords = 2:3, crs = polyCRS))
		}
	}

	if (!inherits(polyList, 'list')) {
		stop('Input must be a list of sf objects.')
	}

	if (!inherits(polyList[[1]], c('sf', 'sfc'))) {
		stop('Input must be a list of sf objects.')
	}
	
	cellType <- match.arg(cellType, c('hexagon', 'square'))
	
	if (!is.null(bb)) {
	    # broaden just a bit for clearer plotting
	    bb[1] <- bb[1] * 0.95
	    bb[2] <- bb[2] * 1.05
	    bb[3] <- bb[3] * 0.95
	    bb[4] <- bb[4] * 1.05
	}
		
	# coarse template
	# For resolution, use latitudinal and longitudinal breadth
	xrange <- abs(getExtentOfList(polyList)['xmax'] - getExtentOfList(polyList)['xmin']) / 100
	yrange <- abs(getExtentOfList(polyList)['ymax'] - getExtentOfList(polyList)['ymin']) / 100
	quickRes <- round((xrange + yrange) / 2, 0)
	
	polyCRS <- sf::st_crs(polyList[[1]])
			
	if (cellType == 'hexagon') {
	    
	    #get overall extent
	    masterExtent <- sf::st_as_sfc(getExtentOfList(polyList))
	    
	    if (!is.null(bb)) {
	        masterExtent <- getExtentOfList(polyList)
	        masterExtent['xmin'] <- bb[1]
	        masterExtent['xmax'] <- bb[2]
	        masterExtent['ymin'] <- bb[3]
	        masterExtent['ymax'] <- bb[4]
	        masterExtent <- sf::st_as_sfc(masterExtent)
	    }
		
		quickTemplate <- sf::st_make_grid(masterExtent, cellsize = c(quickRes, quickRes), square = FALSE, crs = sf::st_crs(masterExtent))
		quickCentroids <- sf::st_centroid(quickTemplate)
		quickTemplate <- sf::st_sf(quickTemplate, grid_id = 1:length(quickTemplate))

		taxonNames <- names(polyList)
		polyList <- lapply(polyList, function(x) sf::st_combine(sf::st_geometry(x)))
		polyList <- do.call(c, polyList)
		
		if (unique(as.character(sf::st_geometry_type(polyList[[1]]))) == 'MULTIPOLYGON') {
			quick <- sf::st_intersects(polyList, quickCentroids)
		} else {
			quick <- sf::st_intersects(polyList, quickTemplate)
		}
		names(quick) <- taxonNames
					
	} else {
		
	    #get overall extent
	    masterExtent <- getExtentOfList(polyList)
	    
	    if (!is.null(bb)) {
	        masterExtent['xmin'] <- bb[1]
	        masterExtent['xmax'] <- bb[2]
	        masterExtent['ymin'] <- bb[3]
	        masterExtent['ymax'] <- bb[4]
	    }
	    
		quickTemplate <- terra::rast(xmin = masterExtent['xmin'], xmax = masterExtent['xmax'], ymin = masterExtent['ymin'], ymax = masterExtent['ymax'], resolution = c(quickRes, quickRes), crs = ifelse(sf::st_is_longlat(polyList[[1]]), '', sf::st_crs(polyList[[1]])$input))
		
		quick <- pbapply::pblapply(polyList, function(x) {
		
			if (all(unique(as.character(sf::st_geometry_type(x))) == 'POINT')) {
			
				terra::cellFromXY(quickTemplate, sf::st_coordinates(x))

			} else {
				
				# rasterize polygon against grid (cells register if midpoint is within polygon)
				xx <- terra::rasterize(terra::vect(x), quickTemplate)

				which(!is.na(terra::values(xx)))
			}
		})
	}
	
	if (inherits(quickTemplate, 'SpatRaster')) {
		nGridCells <- terra::ncell(quickTemplate)
	} else if (inherits(quickTemplate, 'sf')) {
		nGridCells <- nrow(quickTemplate)
	}

	# create site by species matrix
	mat <- matrix(0, nrow = nGridCells, ncol = length(quick))
	colnames(mat) <- names(quick)
	for (i in 1:length(quick)) {
		mat[quick[[i]], names(quick)[i]] <- 1
	}
	
	# remove empty cells if hexagonal grid, and record cell richness
	if (inherits(quickTemplate, 'sf')) {
		quickTemplate <- quickTemplate[which(rowSums(mat) > 0),]
		quickTemplate$spRichness <- rowSums(mat)[which(rowSums(mat) > 0)]
	} else {
		terra::values(quickTemplate) <- rowSums(mat)
		quickTemplate[quickTemplate == 0] <- NA
	}	
	
	message('\n\tAn interactive coarse-grain map has been displayed.\n')
	message('\n\tPlease wait until plot is completed......', appendLF = FALSE)
	
	colramp <- function(n) viridisLite::turbo(n = n, begin = 0.1, end = 0.9)
	
	# add map for context
	if (!sf::st_is_longlat(polyCRS)) {
		wrld <- sf::st_transform(worldmap, crs = polyCRS)
	} else {
		wrld <- worldmap
	}

	wrld <- sf::st_cast(wrld, 'MULTILINESTRING')
	
	if (inherits(quickTemplate, 'sf')) {
		# plot centroid points rather than all hexagonal cells to speed up process
		centroids <- sf::st_centroid(sf::st_geometry(quickTemplate))
		centroids <- sf::st_sf(spRichness = quickTemplate[['spRichness']], geometry = centroids)
		plot(centroids['spRichness'], pch = 20, cex = 1, pal = colramp, key.pos = NULL, main = NULL, reset = FALSE)
		# plot(quickTemplate['spRichness'], pal = colramp, border = NA, key.pos = NULL, main = NULL, reset = FALSE)
	} else {
		terra::plot(quickTemplate, legend = FALSE, axes = FALSE, col = colramp(100))
	}
	
	grXY <- graphics::par("usr")
	clip <- sf::st_make_grid(sf::st_as_sf(rbind.data.frame(grXY[c(1,3)], grXY[c(2,4)]), coords = 1:2, crs = polyCRS), n = 1)
	wrld <- sf::st_intersection(wrld, clip)
	wrld <- sf::st_combine(wrld)				
	# graphics::clip(grXY[1], grXY[2], grXY[3], grXY[4]) # this ensures that world map is constrained to plot region

	graphics::plot(wrld, add = TRUE, lwd = 1, col = gray(0.5))
	
	#ggplot(quickTemplate) + geom_sf(aes(fill = spRichness), lwd = 0.5, show.legend = FALSE) + scale_fill_gradientn(colours = sf::sf.colors()) + geom_sf(data = wrld, fill = NA)
	
	message('done!\n')
	graphics::title(main = 'Define your extent polygon.')

	message('\tClick on the map to create a polygon that will define the extent of the epmGrid.')
	message('\tRight-clicking will close the polygon and terminate the interactive plot.\n\n')
	
	userPoly <- terra::draw(x = 'polygon', col = 'red', xpd = NA)
	userPoly <- terra::as.data.frame(userPoly, geom = "hex")
	userPoly$geometry <- structure(as.list(userPoly$geometry), class = "WKB")
	userPoly <- sf::st_as_sf(userPoly)
	sf::st_crs(userPoly) <- polyCRS
	userPoly <- sf::st_geometry(userPoly)
	
	# display call so user can use this extent in the future
	wkt <- sf::st_as_text(sf::st_geometry(userPoly))
	
	grDevices::dev.off()
	
	return(list(poly = userPoly, wkt = wkt))
}

