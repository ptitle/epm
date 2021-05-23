##' @title Interactively choose extent
##'
##' @description Given a list of polygons, sets up an interactive plot
##' 	to allow the user to draw the desired extent. This can be used to
##' 	define the extent in \link{createEPMgrid}.
##'
##' @param polyList a list of Simple Feature polygons.
##'
##' @param cellType either \code{hexagon} or \code{square}.
##' 
##' @param bb c(xmin, xmax, ymin, ymax) to limit the extent for 
##' the interactive plot.
##'
##' @param nThreads if > 1, then employ parallel computing. 
##'
##' @return A list with a polygon, and its WKT string 
##'
##' @author Pascal Title
##' 
##' @examples
##' \dontrun{
##' ex <- interactiveExtent(tamiasPolyList)
##'
##' # You can use this as the extent in createEPMgrid
##' grid <- createEPMgrid(tamiasPolyList, resolution = 50000, extent = ex$wkt)
##'
##'}
##' @export




interactiveExtent <- function(polyList, cellType = 'square', bb = NULL, nThreads = 1) {
	
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
	xrange <- abs(getExtentOfList(polyList, format = 'sf')['xmax'] - getExtentOfList(polyList, format = 'sf')['xmin']) / 100
	yrange <- abs(getExtentOfList(polyList, format = 'sf')['ymax'] - getExtentOfList(polyList, format = 'sf')['ymin']) / 100
	quickRes <- round((xrange + yrange) / 2, 0)
	
	# if projected, use 100km, if not, use 2 degrees
	# quickRes <- ifelse(sf::st_is_longlat(polyList[[1]]), 2, 100000)
	
	# prep for parallel computing
	if (nThreads == 1) {
		cl <- NULL
	} else if (nThreads > 1) {
		if (.Platform$OS.type != 'windows') {
			cl <- nThreads
		} else {
			cl <- parallel::makeCluster(nThreads)
			parallel::clusterExport(cl, c('method', 'gridTemplate', 'coverCutoff'))
		}
	}	
	
	if (cellType == 'hexagon') {
	    
	    #get overall extent
	    masterExtent <- sf::st_as_sfc(getExtentOfList(polyList, format = 'sf'))
	    
	    if (!is.null(bb)) {
	        masterExtent <- getExtentOfList(polyList, format = 'sf')
	        masterExtent['xmin'] <- bb[1]
	        masterExtent['xmax'] <- bb[2]
	        masterExtent['ymin'] <- bb[3]
	        masterExtent['ymax'] <- bb[4]
	        masterExtent <- sf::st_as_sfc(masterExtent)
	    }
		
		quickTemplate <- sf::st_make_grid(masterExtent, cellsize = c(quickRes, quickRes), square = FALSE, crs = sf::st_crs(masterExtent))
		quickCentroids <- sf::st_centroid(quickTemplate)
		quickTemplate <- sf::st_sf(quickTemplate, grid_id = 1:length(quickTemplate))
		
		quick <- pbapply::pblapply(polyList, function(x) polyToGridCells(x, method = 'centroid', quickTemplate, quickCentroids, subset = 0), cl = cl)
					
	} else {
		
	    #get overall extent
	    masterExtent <- getExtentOfList(polyList, format = 'sf')
	    
	    if (!is.null(bb)) {
	        masterExtent['xmin'] <- bb[1]
	        masterExtent['xmax'] <- bb[2]
	        masterExtent['ymin'] <- bb[3]
	        masterExtent['ymax'] <- bb[4]
	    }
	    
		quickTemplate <- terra::rast(xmin = masterExtent['xmin'], xmax = masterExtent['xmax'], ymin = masterExtent['ymin'], ymax = masterExtent['ymax'], resolution = c(quickRes, quickRes), crs = ifelse(sf::st_is_longlat(polyList[[1]]), '', sf::st_crs(polyList[[1]])$input))
		
		quick <- pbapply::pblapply(polyList, function(x) {
		# for (i in 1:length(polyList)) {
		#     x <- polyList[[i]]
		#     message(i)
		    
			xx <- terra::cells(quickTemplate, y = terra::vect(x), weights = TRUE)
			
			if (nrow(xx) > 0 & !all(is.na(xx[, 'cell']))) {
    			# if polygon extends beyond grid template, cell may be NaN
    			if (anyNA(xx[, 'cell'])) {
    				xx <- xx[!is.na(xx[, 'cell']),]
    			}
    					
    			# cells are returned regardless of how much they are covered by polygon
    			# using centroid method
    			quickCentroids <- terra::xyFromCell(quickTemplate, cell = xx[, 'cell'])
    			quickCentroids <- sf::st_as_sf(as.data.frame(quickCentroids), coords = 1:2, crs = sf::st_crs(x))	
    			
    			xx[, 'cell'][unlist(sf::st_intersects(x, quickCentroids))]
			}
    	}, cl = cl)
	}
	
	if (nThreads > 1 & .Platform$OS.type == 'windows') parallel::stopCluster(cl)
	
	# flip list from list of cells per species, to list of species per cells
	if (inherits(quickTemplate, 'SpatRaster')) {
		cellList <- vector('list', terra::ncell(quickTemplate))
	} else if (inherits(quickTemplate, 'sf')) {
		cellList <- vector('list', nrow(quickTemplate))
	}
	
	for (i in 1:length(quick)) {
		cellList[quick[[i]]] <- lapply(cellList[quick[[i]]], function(x) append(x, names(quick)[i]))
	}
	cellList <- lapply(cellList, unique)

	# add species richness in as attribute
	if (inherits(quickTemplate, 'SpatRaster')) {
		terra::values(quickTemplate) <- lengths(cellList)
		quickTemplate[quickTemplate == 0] <- NA
	} else if (inherits(quickTemplate, 'sf')) {
		quickTemplate$spRichness <- lengths(cellList)
		quickTemplate <- quickTemplate[which(lengths(cellList) > 0),]
	} 
	

	
	# add map for context
	if (!sf::st_is_longlat(polyList[[1]])) {
		wrld <- sf::st_transform(worldmap, crs = sf::st_crs(polyList[[1]]))
	} else {
		wrld <- worldmap
	}
	
	message('\n\tAn interactive coarse-grain map has been displayed.\n')
	message('\n\tPlease wait until plot is completed......', appendLF = FALSE)
	
	if (inherits(quickTemplate, 'sf')) {
		plot(quickTemplate['spRichness'], border = NA, key.pos = NULL, main = NULL, reset = FALSE)
	} else {
		terra::plot(quickTemplate, legend = FALSE, axes = FALSE, col = grDevices::colorRampPalette(c('blue', 'cyan', 'yellow', 'red'))(100))
	}

	graphics::plot(wrld, add = TRUE, lwd = 1, border = gray(0.5))
	
	#ggplot(quickTemplate) + geom_sf(aes(fill = spRichness), lwd = 0.5, show.legend = FALSE) + scale_fill_gradientn(colours = sf::sf.colors()) + geom_sf(data = wrld, fill = NA)
	
	message('done!\n')
	graphics::title(main = 'Define your extent polygon.')

	message('\tClick on the map to create a polygon that will define the extent of the epmGrid.')
	message('\tRight-clicking will close the polygon and terminate the interactive plot.\n\n')
	
	userPoly <- terra::draw(x = 'polygon', col = 'red', xpd = NA)
	userPoly <- terra::as.data.frame(userPoly, geom = "hex")
	userPoly$geometry <- structure(as.list(userPoly$geometry), class = "WKB")
	userPoly <- sf::st_as_sf(userPoly)
	sf::st_crs(userPoly) <- sf::st_crs(polyList[[1]])
	
	# display call so user can use this extent in the future
	wkt <- sf::st_as_text(sf::st_geometry(userPoly))
	
	grDevices::dev.off()
	
	return(list(poly = userPoly, wkt = wkt))
}

