##' @title Create epmGrid object
##' 
##' @description Takes a list of polygons and creates a epmGrid object.
##'
##'
##' @param polyList a list of polygon objects (sf or sp), named with taxon names.
##'
##' @param resolution vertical and horizontal spacing of grid cells, in units 
##'		of the polygons' projection.
##'
##' @param method approach used for gridding. Either \code{centroid} or \code{areaCutoff}. See details below.
##'
##' @param cellType either \code{hexagon} or \code{square}. See details below.
##' 
##' @param coverCutoff the percent that a species range must cover a grid cell to be considered present.
##'
##' @param retainSmallRanges boolean; should small ranged species be dropped or preserved.
##'		See details.
##'
##' @param extent if 'auto', then the maximal extent of the polygons will be used. 
##' 	If not 'auto', can be a SpatialPolygon or sf object, in which case the resulting epmGrid
##'		will be cropped and masked with respect to the polygon, or a spatial coordinates object, 
##' 	from which an extent object will be generated, or a numeric vector of length 4 
##' 	with minLong, maxLong, minLat, maxLat. If 'interactive', then an interactive plot
##'		will appear in which the user can draw the desired polygon extent. That extent will then be returned.
##'
##' @param checkValidity if \code{TRUE}, then check polygon validity and repair if needed, 
##' 		using sf::st_make_valid and the lwgeom package. 
##'
##' @param nThreads if > 1, then employ parallel computing. 
##' 
##' @param nGroups break up the grid into this many groups for processing. This will alleviate memory 
##' 		usage for datasets that are very high resolution (not relevant for large numbers of species)
##'
##'
##' @details 
##'		If \code{cellType = 'hexagon'}, then the grid is made of polygons via the sf package.
##' 	If \code{cellType = 'square'}, then the grid is a raster generated via the terra package.
##' 	Hexagonal cells have several advantages, including being able to be of different sizes (if the grid is in 
##' 	unprojected long/lat), and may be able to more naturally follow coastlines and non-linear features.
##'		However, the raster-based square cells will be much less memory intensive for high resolution datasets. 
##' 	Choice of grid type matters more for spatial resolution, than for number of species.
##'
##' 
##' 	In the polygon-to-grid conversion process, two approaches are implemented. 
##'		For \code{method = 'centroid'}, a range polygon registers in a cell if the polygon 
##' 	overlaps with the cell centroid. 
##' 	For \code{method = 'areaCutoff'}, a range polygon registers in a cell if it covers that cell by 
##' 	at least \code{coverCutoff} fraction of the cell.
##' 	If \code{retainSmallRanges = FALSE}, then species whose ranges are so small that no 
##' 	cell registers as present will be dropped. If \code{retainSmallRanges = TRUE}, then the 
##' 	cell that contains the majority of the the small polygon will be considered as present, even if it's a 
##' 	small percent of the cell.
##' 
##'		If \code{retainSmallRanges = TRUE}, and an extent is provided, then species may still be dropped if they
##' 	fall outside of that extent.
##'
##'		In interactive mode for defining the extent, the user can draw a bounding polygon on a 
##'		map. The drawn polygon will then be printed to the console so that the user can provide 
##'		that bounding polygon in future calls as the extent.
##'
##'		Any SpatialPolygon or SpatialPoints objects are converted to objects of class \code{sf}.
##'	
##' 
##' 	In interactive mode, the basemap is from \url{www.naturalearthdata.com}. 
##' 
##' @return an object of class \code{epmGrid}. If \code{extent = 'interactive'}, then a polygon is returned.
##' 
##' @author Pascal Title
##'
##' @examples
##' library(sf)
##' # example dataset: a list of 24 chipmunk distributions as polygons
##' head(tamiasPolyList)
##' 
##' # hexagonal grid
##' tamiasEPM <- createEPMgrid(tamiasPolyList, resolution = 50000, 
##' 	cellType = 'hexagon', method = 'centroid')
##' tamiasEPM
##'
##' # square grid
##' tamiasEPM2 <- createEPMgrid(tamiasPolyList, resolution = 50000, 
##' 	cellType = 'square', method = 'centroid')
##' tamiasEPM2
##' 
##' @export

createEPMgrid <- function(polyList, resolution = 50000, method = 'centroid', cellType = 'hexagon', coverCutoff = 0.1, retainSmallRanges = TRUE, extent = 'auto', checkValidity = FALSE, nThreads = 1, nGroups = 1) {
	
	# polyList <- tamiasPolyList; resolution = 50000; method = 'centroid'; cellType = 'hexagon'; coverCutoff = 0.1; retainSmallRanges = TRUE; extent = 'auto'; checkValidity = FALSE; nThreads = 1; nGroups = 1
	
	if (is.list(polyList)) {
		if (inherits(polyList[[1]], c('SpatialPolygons', 'SpatialPolygonsDataFrame'))) {
			# if class SpatialPolygons, convert to sf
			for (i in 1:length(polyList)) {
				polyList[[i]] <- sf::st_as_sf(polyList[[i]])
			}
		}
	} else if (!inherits(polyList[[1]], c('SpatialPolygons', 'SpatialPolygonsDataFrame', 'sf', 'sfc'))) {
		stop('polyList must be a list of SpatialPolygons or Simple Features.')
	}
	
	cellType <- match.arg(cellType, c('hexagon', 'square'))
	
	if (!cellType %in% c('hexagon', 'square')) {
		stop('cellType must be either hexagon or square.')
	}

	if (is.null(names(polyList))) {
		stop('List must be named with species names.')
	}
	
	if (anyDuplicated(names(polyList))) {
		stop('Some taxon names are duplicated.')
	}
	
	if (coverCutoff > 1) {
		stop('coverCutoff is a fraction.')
	}
	
	if (nThreads > parallel::detectCores()) {
		stop('nThreads is greater than what is available.')
	}
	
	# if lwgeom package available, and checking is requested, check sf polygon validity. 
	if (!requireNamespace('lwgeom', quietly = TRUE) & checkValidity) {
		stop('For polygon validity checking, the lwgeom package is required.')
	}
	
	if (requireNamespace('lwgeom', quietly = TRUE) & checkValidity) {
		for (i in 1:length(polyList)) {
			if (inherits(polyList[[i]], 'sf')) {
				if (!any(sf::st_is_valid(polyList[[i]]))) {
					message('\trepairing polygon ', i)
					polyList[[i]] <- sf::st_make_valid(polyList[[i]])
				}
			}
		}
	}

	# test that all have same CRS
	if (!any(sapply(polyList, function(x) identical(sf::st_crs(polyList[[1]]), sf::st_crs(x))))) {
		stop('proj4string or EPSG of all polygons must match.')
	}
	proj <- sf::st_crs(polyList[[1]])
	
	# if WKT string, then convert to sf polygon
	if (inherits(extent, 'character')) {
		if (grepl('POLYGON \\(', extent)) {
			extent <- sf::st_as_sfc(extent, crs = sf::st_crs(polyList[[1]]))
		}
	}

	if (inherits(extent, 'character')) {
		
		if (!extent %in% c('auto', 'interactive')) {
			stop("If extent is a character vector, it can only be 'auto' or 'interactive'.")
		}
		
		if (extent == 'auto') {
		
			#get overall extent
			masterExtent <- getExtentOfList(polyList, format = 'sf')
		}
	} else if (is.numeric(extent) & length(extent) == 4) {
		# use user-specified bounds
		masterExtent <- sf::st_bbox(polyList[[1]])
		masterExtent[[1]] <- extent[1]
		masterExtent[[2]] <- extent[3]
		masterExtent[[3]] <- extent[2]
		masterExtent[[4]] <- extent[4]
		
	} else if (inherits(extent, c('SpatialPolygons', 'SpatialPolygonsDataFrame', 'SpatialPoints', 'SpatialPointsDataFrame', 'sf', 'sfc'))) {
		
		if (inherits(extent, c('SpatialPolygons', 'SpatialPolygonsDataFrame', 'SpatialPoints', 'SpatialPointsDataFrame'))) {
			extent <- sf::st_as_sf(extent)
		}
			
		if (!is.na(sf::st_crs(extent))) {
			if (!identical(sf::st_crs(extent), proj)) {
				extent <- sf::st_transform(extent, crs = proj)
			}
		} else {
			sf::st_crs(extent) <- proj
		}
	  
		if (inherits(extent, 'sfc')) {
			extent <- sf::st_sf(extent)
		}
		
		# get extent from spatial object
		masterExtent <- sf::st_bbox(extent)

	} else if (inherits(extent, 'Extent')) {
		
		masterExtent <- sf::st_bbox(polyList[[1]])
		masterExtent[[1]] <- extent@xmin
		masterExtent[[2]] <- extent@ymin
		masterExtent[[3]] <- extent@xmax
		masterExtent[[4]] <- extent@ymax
		
	} else {
		stop("extent must be 'auto', a spatial object or a vector with minLong, maxLong, minLat, maxLat.")
	}
	
	# interactive extent: if this option is selected, a coarse richness grid
	# will be plotted, so that the user can designate an extent polygon
	wkt <- NULL
	if (inherits(extent, 'character')) {
		if (extent == 'interactive') {

			interactive <- interactiveExtent(polyList)
			extent <- interactive$poly
			wkt <- interactive$wkt
			
			return(wkt)
		}
	}
	
	
	# here, we take one of two routes:
	## if hexagonal grid, then use sf, if square grid, use terra
	
	uniqueSp <- sort(names(polyList))
	polyList <- polyList[uniqueSp]

	
	if (cellType == 'hexagon') {
		
		spGridList <- polyToHex(poly = polyList, method = method, coverCutoff = coverCutoff, extentVec = masterExtent, resolution = resolution, crs = proj, nGroups = nGroups, retainSmallRanges = retainSmallRanges, nThreads = nThreads)
		
		
	} else if (cellType == 'square') {
		
		spGridList <- polyToTerra(poly = polyList, method = method, coverCutoff = coverCutoff, extentVec = masterExtent, resolution = resolution, crs = proj, retainSmallRanges = retainSmallRanges, nThreads = nThreads)
		
	} else {
		stop('Cell type not supported.')
	}
	
	gridTemplate <- spGridList[[1]]
	spGridList <- spGridList[[2]]
	

	# update
	smallSp <- which(lengths(spGridList) == 0)
	
	# if small ranged species are not to be preserved, drop them
	# OR if some species are outside of proposed extent, drop them
	if (length(smallSp) > 0) {
		spGridList <- spGridList[ - smallSp]

		msg <- paste0('The following species are being dropped:\n\t', paste(names(smallSp), collapse='\n\t'))
		message(msg)
	}
	
	# flip list from list of cells per species, to list of species per cells
	if (inherits(gridTemplate, 'SpatRaster')) {
		cellList <- vector('list', terra::ncell(gridTemplate))
	} else if (inherits(gridTemplate, 'sf')) {
		cellList <- vector('list', nrow(gridTemplate))
	} else {
		stop('gridTemplate class mismatch.')
	}
	for (i in 1:length(spGridList)) {
		cellList[spGridList[[i]]] <- lapply(cellList[spGridList[[i]]], function(x) append(x, names(spGridList)[i]))
	}
	cellList <- lapply(cellList, unique)
	
	uniqueComm <- unique(cellList)
	uniqueComm2 <- sapply(uniqueComm, function(y) paste(y, collapse = '|'))
	cellList2 <- sapply(cellList, function(y) paste(y, collapse = '|'))

	cellCommVec <- integer(length = length(cellList))
	for (i in 1:length(uniqueComm)) {
		cellCommVec[which(cellList2 == uniqueComm2[i])] <- i
	}
	
	rm(cellList2, uniqueComm2)
	
	# set empty cells to NA, rather than NULL
	uniqueComm[sapply(uniqueComm, is.null)] <- NA

	# add unique community index and species richness in as attributes
	if (inherits(gridTemplate, 'SpatRaster')) {
		gridTemplate <- c(gridTemplate, terra::rast(gridTemplate))
		terra::values(gridTemplate) <- cbind(cellCommVec, lengths(cellList))
		names(gridTemplate) <- c('uniqueComm', 'spRichness')
		gridTemplate$spRichness[gridTemplate$spRichness == 0] <- NA
	} else if (inherits(gridTemplate, 'sf')) {
		gridTemplate$uniqueComm <- cellCommVec
		gridTemplate$spRichness <- lengths(cellList)
	} else {
		stop('gridTemplate class mismatch.')
	}
	
	# remove empty gridcells
	if (inherits(gridTemplate, 'sf')) {
		gridTemplate <- gridTemplate[which(lengths(cellList) > 0),]
		cellCommVec <- cellCommVec[which(lengths(cellList) > 0)]
	}

	# prepare output object
	obj <- vector('list', length = 7)
	names(obj) <- c('grid', 'speciesList', 'cellCommInd', 'geogSpecies', 'cellCount', 'data', 'phylo')

	obj[['grid']] <- gridTemplate
	obj[['speciesList']] <- uniqueComm
	obj[['cellCommInd']] <- cellCommVec
	obj[['geogSpecies']] <- uniqueSp
	obj[['cellCount']] <- lengths(spGridList)
	attr(obj, 'resolution') <- resolution
	if (inherits(gridTemplate, 'sf')) {
		attr(obj, 'crs') <- sf::st_crs(gridTemplate)$input
	} else {
		attr(obj, 'crs') <- terra::crs(gridTemplate, proj4 = TRUE)
	}
	if (inherits(gridTemplate, 'sf')) {
		attr(obj, 'projected') <- !sf::st_is_longlat(gridTemplate)
	} else {
		attr(obj, 'projected') <- !sf::st_is_longlat(proj)
	}
	attr(obj, 'gridType') <- cellType
	attr(obj, 'metric') <- 'spRichness'
	
	class(obj) <- 'epmGrid'
	return(obj)	
}


# Function to take a polygon and return occupied grid cells
# subset allows us to operate on a subset of grid cells at a time, if necessary
polyToGridCells <- function(poly, method, gridTemplate, gridCentroids, coverCutoff, subset) {

	if (subset == 0) {
		clusterInd <- 1:nrow(gridTemplate)
	} else {
		clusterInd <- which(gridTemplate$clusters == subset)
	}
	
	if (method == 'centroid') {
		tmp <- clusterInd[unlist(sf::st_intersects(poly, gridCentroids[clusterInd, ]))]
	}

	if (method == 'areaCutoff') {
		tmp <- unlist(sf::st_intersects(poly, gridTemplate[clusterInd, ]))
		tmpWithin <- clusterInd[tmp][lengths(sf::st_within(gridTemplate[clusterInd[tmp],], poly)) > 0]
		tmpCrosses <- setdiff(clusterInd[tmp], tmpWithin)
		# tmpWithin <- tmp[lengths(sf::st_within(gridTemplate[tmp,], poly)) > 0]
		# tmpCrosses <- setdiff(tmp, tmpWithin)
				
		if (length(tmpWithin) > 0 & length(tmpCrosses) > 0) {
			tmpCrosses <- tmpCrosses[as.vector(sf::st_area(sf::st_intersection(sf::st_union(poly), gridTemplate[tmpCrosses,])) / sf::st_area(gridTemplate[tmpCrosses,])) >= coverCutoff]
			tmp <- union(tmpWithin, tmpCrosses)
		}
	}
	
	return(tmp)
}



## IDEA: if hexagons, use sf, if square, use terra package. Creates same type of epm object, except that first slot contains a sf object or spatRaster object. Then, other functions will have to be able to handle either. 
## best of both worlds: hexagons and squares, and slow vs fast options.








## Hexagonal grid cells via the sf package. 
### if method == 'centroid', then species belongs to cell if its range intersects a cell's centroid.
### if method == 'areaCutoff', then species belongs to cell if it covers some specified percent of cell area.

## if retainSmallSpecies == TRUE, then the cell that has the greatest occupancy for a species (regardless of the amount) is marked as present. 

polyToHex <- function(poly, method, coverCutoff, extentVec, resolution, crs, nGroups, retainSmallRanges, nThreads) {
	
	# Generate template
	masterExtent <- sf::st_as_sfc(extentVec)
	gridTemplate <- sf::st_make_grid(masterExtent, cellsize = c(resolution, resolution), square = FALSE, crs = sf::st_crs(masterExtent))
	
	# if extent was polygon, then mask the grid template
	if (inherits(masterExtent, c('sf', 'sfc'))) {
		gridTemplate <- gridTemplate[lengths(sf::st_intersects(gridTemplate, masterExtent)) > 0,]
	}	
	
	gridTemplate <- sf::st_sf(gridTemplate, grid_id = 1:length(gridTemplate))
	
	gridCentroids <- sf::st_centroid(gridTemplate)
	
	# do kmeans clustering on centroid coords.
	if (nGroups > 1) {
		clusters <- stats::kmeans(sf::st_coordinates(gridCentroids), centers = nGroups)
		gridTemplate$clusters <- clusters$cluster
	}
	
	# prep for parallel computing
	if (nThreads == 1) {
		cl <- NULL
	} else if (nThreads > 1) {
		if (.Platform$OS.type != 'windows') {
			cl <- nThreads
		} else {
			cl <- parallel::makeCluster(nThreads)
			parallel::clusterExport(cl, c('polyToGridCells', 'method', 'gridTemplate', 'gridCentroids', 'coverCutoff'))
		}
	}
	
	if (nGroups > 1) {
		subsetsList <- list()
		for (i in 1:nGroups) {
			
			message('\tProcessing subgroup ', i, ' of ', nGroups, '...')
			
			# get grid cells per species
			subsetsList[[i]] <- pbapply::pblapply(poly, function(x) polyToGridCells(x, method = method, gridTemplate, gridCentroids, coverCutoff, subset = i), cl = cl)
		}
		
		spGridList <- lapply(1:length(subsetsList[[1]]), function(i) {
			unlist(lapply(subsetsList, function(x) x[[i]]))
		})
		names(spGridList) <- names(subsetsList[[1]])
		rm(subsetsList, clusters)	
	
	
	} else {
				
		spGridList <- pbapply::pblapply(poly, function(x) polyToGridCells(x, method = method, gridTemplate, gridCentroids, coverCutoff, subset = 0), cl = cl)
		
	}

	if (nThreads > 1 & .Platform$OS.type == 'windows') parallel::stopCluster(cl)
		
	# if we are retaining small ranged species that would otherwise be dropped,
	# then we will search for species that did not register in any grid cell
	smallSp <- which(lengths(spGridList) == 0)
	
	if (retainSmallRanges) {
		if (length(smallSp) > 0) {
			for (i in 1:length(smallSp)) {

				tmp <- unlist(sf::st_intersects(poly[[smallSp[i]]], gridTemplate))
				
				# of grid cells intersected by small-ranged species, keep the cell most occupied.
				# this is to avoid going from a species that would not appear in any cell, to a species occuring in multiple cells with 1% coverage.
				tmp <- tmp[which.max(as.vector(sf::st_area(sf::st_intersection(sf::st_union(poly[[smallSp[i]]]), gridTemplate[tmp,])) / sf::st_area(gridTemplate[tmp,])))]
				
				spGridList[[smallSp[i]]] <- tmp
			}
			rescued <- setdiff(names(smallSp), names(which(lengths(spGridList) == 0)))
			if (length(rescued) > 0) {
				msg <- paste0(length(rescued), ' small-ranged species ', ifelse(length(rescued) > 1, 'were', 'was'), ' preserved:\n\t', paste(rescued, collapse='\n\t'))
				message(msg)
			}
		}
	}
	
	return(list(gridTemplate, spGridList))
}	
	










# square gridcells via the terra package
 
### if method == 'centroid', then species belongs to cell if its range intersects a cell's centroid.
### if method == 'areaCutoff', then species belongs to cell if it covers some specified percent of cell area.

## if retainSmallSpecies == TRUE, then the cell that has the greatest occupancy for a species (regardless of the amount) is marked as present. 


polyToTerra <- function(poly, method, coverCutoff, extentVec, resolution, crs, retainSmallRanges, nThreads) {

	# Generate template
	gridTemplate <- terra::rast(xmin = extentVec['xmin'], xmax = extentVec['xmax'], ymin = extentVec['ymin'], ymax = extentVec['ymax'], resolution = c(resolution, resolution), crs = crs$input)

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


	spGridList <- pbapply::pblapply(poly, function(x) {
		
		xx <- terra::cells(gridTemplate, y = terra::vect(x), weights = TRUE)
				
		# cells are returned regardless of how much they are covered by polygon

		if (method == 'centroid') {
			centroids <- terra::xyFromCell(gridTemplate, cell = xx[, 'cell'])
			centroids <- sf::st_as_sf(as.data.frame(centroids), coords = 1:2, crs = sf::st_crs(x))	
			presenceCells <- xx[, 'cell'][unlist(sf::st_intersects(x, centroids))]

		} else if (method == 'areaCutoff') {

			presenceCells <- xx[xx[, 'weights'] >= coverCutoff, 'cell']
			# terra::xyFromCell(gridTemplate, cell = presenceCells) # for testing
			
		}		
		
		return(presenceCells)
	}, cl = cl)
	
	if (nThreads > 1 & .Platform$OS.type == 'windows') parallel::stopCluster(cl)
	
	## for testing
	# plot(st_geometry(poly[[i]]))
	# testGrid <- rast(gridTemplate)
	# testGrid[xx[, 'cell']] <- 1
	# plot(as.polygons(testGrid, dissolve = F), col = 'orange')	
	# plot(st_geometry(poly[[i]]), add=T)


	# if we are retaining small ranged species that would otherwise be dropped,
	# then we will search for species that did not register in any grid cell
	smallSp <- which(lengths(spGridList) == 0)
	
	if (retainSmallRanges) {
		if (length(smallSp) > 0) {
			
			gridTemplateHighRes <- terra::rast(xmin = extentVec['xmin'], xmax = extentVec['xmax'], ymin = extentVec['ymin'], ymax = extentVec['ymax'], resolution = c(resolution/10, resolution/10), crs = crs$input)
			
			for (i in 1:length(smallSp)) {
				
				xx <- terra::cells(gridTemplateHighRes, y = terra::vect(poly[[smallSp[i]]]), weights = TRUE)
				
				smallCells <- xx[which.max(xx[, 'weights']), 'cell']
				smallCells <- terra::cellFromXY(gridTemplate, terra::xyFromCell(gridTemplateHighRes, cell = smallCells))
				spGridList[[smallSp[i]]] <- smallCells

				# ## alternative approach using point sampling
				# xx <- terra::cells(gridTemplateHighRes, y = terra::vect(poly[[smallSp[i]]]), weights = TRUE)
				# if (nrow(xx) > 0) {
					# pts <- terra::xyFromCell(gridTemplateHighRes, cell = xx[, 'cell'])
				# } else {
					# # if range is too small to register, sample points within the polygon
					# pts <- sf::st_sample(poly[[smallSp[i]]], size = 100)
					# pts <- sf::st_coordinates(pts)
				# }
				# smallCells <- unique(terra::cellFromXY(gridTemplate, pts))
				# spGridList[[smallSp[i]]] <- smallCells		
			}
			
			rescued <- setdiff(names(smallSp), names(which(lengths(spGridList) == 0)))
			if (length(rescued) > 0) {
				msg <- paste0(length(rescued), ' small-ranged species ', ifelse(length(rescued) > 1, 'were', 'was'), ' preserved:\n\t', paste(rescued, collapse='\n\t'))
				message(msg)
			}
		}
	}
	
	return(list(gridTemplate, spGridList))
}

				
				
				
		



