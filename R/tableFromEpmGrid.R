##'@title Data table from epmGrid
##'
##'@description Given one or several epmGrid objects, sf objects, rasterLayers,
##'  SpatRasters, create a table of values and associated coordinate data.
##'
##'@param ... objects of class \code{epmGrid}, \code{sf}, \code{sp},
##'  \code{SpatRaster}, \code{RasterLayer} or \code{RasterStack}. All should
##'  have the same projection.
##'@param n number of cells to randomly subsample, no subsampling if \code{NULL}
##'@param minTaxCount integer; cells with at least this many taxa will be included. 
##'@param coords if NULL, then points are sampled as needed, otherwise, data
##'  will be extracted at these specified coordinates.
##'@param id boolean, should the grid cell index (of the first item in the inputs) be returned as well?
##'
##'@details A set of cells are identified in the input objects. If
##'  \code{n=NULL}, then all cells are used, otherwise cells are randomly
##'  subsampled. Values at those cells are then returned. This table
##'  construction can be particularly useful for subsequent statistical
##'  analyses.
##'
##'  Only cells with data in all inputs are returned. If n is greater than the
##'  number of cells with data, then fewer than n cells will be returned.
##'
##'  The first element provided should be a \code{epmGrid} object, and that will
##'  be the one used as a template for the sampled grid system.
##'
##'  If \code{coords} is provided, then data are extracted at those coordinates,
##'  and no subsetting of those points is done.
##'
##'@return data.frame with input variables, as well as \code{"x"} and
##'  \code{"y"}.
##'
##'@author Pascal Title
##'
##' @examples
##'
##' tamiasEPM
##' tamiasEPM <- addPhylo(tamiasEPM, tamiasTree)
##' tamiasEPM <- addTraits(tamiasEPM, tamiasTraits)
##' morphoDisp <- gridMetrics(tamiasEPM, metric='disparity')
##' meanPat <- gridMetrics(tamiasEPM, metric='meanPatristic')
##'
##' tableFromEpmGrid(tamiasEPM, morphoDisp, meanPat, n = 100, 
##' minTaxCount = 2)
##'
##' # this time request grid cell ID's, which would be useful
##' # for linking this table back to the grid system
##' tableFromEpmGrid(tamiasEPM, morphoDisp, meanPat, n = 100, 
##' minTaxCount = 2, id = TRUE)
##'
##' # from predetermined set of coordinates
##' pts <- sf::st_sample(tamiasEPM[[1]], size = 10)
##' tableFromEpmGrid(tamiasEPM, morphoDisp, meanPat, n = 100, 
##' minTaxCount = 1, coords = pts)
##'
##'@export

tableFromEpmGrid <- function(..., n = NULL, minTaxCount = 1, coords = NULL, id = FALSE) {
    
    # x <- list(tamiasEPM, morphoDisp, meanPat); n = 100; minTaxCount = 2; coords = NULL; id = TRUE
	# inputNames <- c('richness', 'morphdisp', 'meanpat')

	# capture input names
	inputNames <- as.list(substitute(list(...)))[-1L]
	
	x <- list(...)
	
	if (!all(sapply(x, inherits, c('epmGrid', 'sf', 'SpatialPolygonsDataFrame', 'RasterLayer', 'RasterStack', 'SpatRaster')))) {
		stop('Not all inputs are of an acceptable class.')
	}
	
	# coerce any RasterLayer or RasterStack objects to SpatRaster, and any sp objects to sf
	for (i in 1:length(x)) {
		if (inherits(x[[i]], c('RasterStack', 'RasterLayer'))) {
			x[[i]] <- as(x[[i]], 'SpatRaster')
		}
		
		if (inherits(x[[i]], 'SpatialPolygonsDataFrame')) {
			x[[i]] <- sf::st_as_sfc(x[[i]], crs = sf::st_crs(x[[i]]))
		}
	}
		
	# check projections
	xproj <- vector('list', length(x))
	for (i in 1:length(x)) {
		if (inherits(x[[i]], 'epmGrid')) {
			xproj[[i]] <- attributes(x[[i]])$crs
		} else {
			xproj[[i]] <- sf::st_crs(x[[i]])$input
		}
	}
	# if (length(unique(xproj)) > 1) {
		# stop('Not all inputs have the same projection.')
	# }

	if (is.null(coords)) {
		# first item will be used as a template
		if (inherits(x[[1]], 'epmGrid')) {
			if (minTaxCount > 1) {
				tooFewInd <- spCountIndex(x[[1]], count = 1:(minTaxCount - 1))
				if (inherits(x[[1]][[1]], 'sf')) {
					grid_multiSp <- x[[1]][[1]][- tooFewInd,]
				} else {
					grid_multiSp <- x[[1]][[1]][attributes(x[[1]])$metric]
					grid_multiSp[tooFewInd] <- NA
				}
				rm(tooFewInd)
			} else {
				if (inherits(x[[1]][[1]], 'sf')) {
					grid_multiSp <- x[[1]][[1]]
				} else {
					grid_multiSp <- x[[1]][[1]][attributes(x[[1]])$metric]
				}
			}
			if (inherits(x[[1]][[1]], 'sf')) {
				templateCentroids <- sf::st_centroid(sf::st_geometry(grid_multiSp))
			} else {
				templateCentroids <- sf::st_geometry(sf::st_as_sf(terra::as.points(grid_multiSp, values = FALSE, na.rm = TRUE)))
			}
			rm(grid_multiSp)
		} else if (inherits(x[[1]], 'sf')) {
			templateCentroids <- sf::st_centroid(sf::st_geometry(x[[1]]))
		} else if (inherits(x[[1]], 'SpatRaster')) {
			rasterSum <- terra::app(x[[1]], fun = sum)
			templateCentroids <- sf::st_geometry(sf::st_as_sf(terra::as.points(rasterSum, values = FALSE, na.rm = TRUE)))
		} else if (inherits(x[[1]], 'SpatialPolygonsDataFrame')) {
			tmp <- sf::st_as_sfc(x[[1]], crs = sf::st_crs(x[[1]]))
			templateCentroids <- sf::st_centroid(tmp)
			rm(tmp)
		}
		
		# for each other item, determine which cells intersect with template, and which of those have data
		resList <- vector('list', length(x))
	
		for (i in 2:length(x)) {
			
			if (inherits(x[[i]], 'epmGrid')) {
				if (minTaxCount > 1) {				
				    tooFewInd <- spCountIndex(x[[i]], count = 1:(minTaxCount - 1))
					if (inherits(x[[i]][[1]], 'sf')) {
						grid_multiSp <- x[[i]][[1]][- tooFewInd,]
					} else {
						grid_multiSp <- x[[i]][[1]][attributes(x[[i]])$metric]
						grid_multiSp[tooFewInd] <- NA
					}
					rm(tooFewInd)
				} else {
					if (inherits(x[[i]][[1]], 'sf')) {
						grid_multiSp <- x[[i]][[1]]
					} else {
						grid_multiSp <- x[[i]][[1]][attributes(x[[i]])$metric]
					}
				}
				
				# get the points that intersect with templateCentroids and that have data
				if (inherits(grid_multiSp, 'sf')) {
					tmp <- sf::st_join(sf::st_as_sf(templateCentroids), grid_multiSp, join = sf::st_within)
					resList[[i]] <- which(!is.na(tmp[[attributes(x[[i]])$metric]]))
				} else {
					tmp <- terra::extract(grid_multiSp, sf::st_coordinates(templateCentroids))
					resList[[i]] <- which(!is.na(tmp))
				}	
			
			} else if (inherits(x[[i]], 'sf')) {
				datCol <- setdiff(colnames(x[[i]]), attributes(x[[i]])$sf_column)
				if (length(datCol) > 1) {
					stop(paste('List item', i, 'contains more than one attribute.'))
				}
				tmp <- sf::st_join(sf::st_as_sf(templateCentroids), x[[i]], join = sf::st_within)
				resList[[i]] <- which(!is.na(tmp[[datCol]]))
	
			} else if (inherits(x[[i]], 'SpatRaster')) {
				rasterSum <- terra::app(x[[i]], fun = sum)
				tmp <- terra::extract(rasterSum, sf::st_coordinates(templateCentroids))
				resList[[i]] <- which(!is.na(tmp))
						
			} else if (inherits(x[[i]], 'SpatialPolygonsDataFrame')) {
				if (ncol(x[[i]]@data) > 1) {
					stop(paste('List item', i, 'contains more than one attribute.'))
				}
				tmp <- sf::st_as_sfc(x[[i]], crs = sf::st_crs(x[[i]]))
				datCol <- setdiff(colnames(tmp), attributes(tmp)$sf_column)
				tmp <- sf::st_join(sf::st_as_sf(templateCentroids), tmp, join = sf::st_within)
				resList[[i]] <- which(!is.na(tmp[[datCol]]))
	
			} else {
				stop(paste('Class of item', i, 'not expected.'))
			}
		}
		
		# each list element contains the indices of template coords where element i has data.
		# so we will reduce this to the intersection of all elements to get the coordinates where all items have data.
		resList <- resList[ - 1] # empty
		goodCells <- Reduce(intersect, resList)
		
		# get subset if requested
		if (!is.null(n)) { 
			if (n < length(goodCells)) {
				goodCells <- sample(goodCells, size = n, replace = FALSE)
			}
		}
		
		gridTemplate <- templateCentroids[goodCells]
		gridTemplate <- sf::st_geometry(gridTemplate)
		
		
	} else {
		
		# coordinates have been provided.
		## Will accept either a matrix of longitude and latitude or a sf points object
		if (inherits(coords, c('sf', 'sfc'))) {
			coords <- sf::st_coordinates(coords)
		}
		if (ncol(coords) != 2 | !is.numeric(coords)) {
			stop('Coords expected to be a 2-column matrix of coordinates.')
		}
		
		gridTemplate <- sf::st_geometry(sf::st_as_sf(as.data.frame(coords), coords = 1:2, crs = xproj[[1]]))
	}
	
	# extract the information from each item
	df <- as.data.frame(matrix(nrow = length(gridTemplate), ncol = length(x) + 2))
	colnames(df) <- c('x', 'y', inputNames)

	for (i in 1:length(x)) {

		if (inherits(x[[i]], 'epmGrid')) {
			if (minTaxCount > 1) {				
			    tooFewInd <- spCountIndex(x[[i]], count = 1:(minTaxCount - 1))
				if (inherits(x[[i]][[1]], 'sf')) {
					grid_multiSp <- x[[i]][[1]][- tooFewInd,]
					grid_multiSp <- grid_multiSp[attributes(x[[i]])$metric]
				} else {
					grid_multiSp <- x[[i]][[1]][attributes(x[[i]])$metric]
					grid_multiSp[tooFewInd] <- NA
				}
				rm(tooFewInd)
			} else {
				if (inherits(x[[i]][[1]], 'sf')) {
					grid_multiSp <- x[[i]][[1]][attributes(x[[i]])$metric]
				} else {
					grid_multiSp <- x[[i]][[1]][attributes(x[[i]])$metric]
				}
			}
			
			if (inherits(grid_multiSp, 'sf')) {
				tmp <- sf::st_join(sf::st_as_sf(gridTemplate), grid_multiSp, join = sf::st_within)
				df[, i + 2] <- sf::st_drop_geometry(tmp)[attributes(x[[i]])$metric]
			} else {
				df[, i + 2] <- terra::extract(grid_multiSp, sf::st_coordinates(gridTemplate))
			}
		
		} else if (inherits(x[[i]], 'sf')) {
			datCol <- setdiff(colnames(x[[i]]), attributes(x[[i]])$sf_column)
			if (length(datCol) > 1) {
				stop(paste('List item', i, 'contains more than one attribute.'))
			}
			tmp <- sf::st_join(sf::st_as_sf(gridTemplate), x[[i]], join = sf::st_within)
			df[, i + 2] <- sf::st_drop_geometry(tmp)[datCol]

		} else if (inherits(x[[i]], 'SpatRaster')) {
			for (j in 1:terra::nlyr(x[[i]])) {
				df[, i + 2] <- terra::extract(x[[i]][[j]], sf::st_coordinates(gridTemplate))
			}
					
		} else if (inherits(x[[i]], 'SpatialPolygonsDataFrame')) {
			if (ncol(x[[i]]@data) > 1) {
				stop(paste('List item', i, 'contains more than one attribute.'))
			}
			tmp <- sf::st_as_sfc(x[[i]], crs = sf::st_crs(x[[i]]))
			datCol <- setdiff(colnames(tmp), attributes(tmp)$sf_column)
			tmp <- sf::st_join(sf::st_as_sf(templateCentroids), tmp, join = sf::st_within)
			resList[[i]] <- which(!is.na(tmp[[datCol]]))
	
		} else {
			stop(paste('Class of item', i, 'not expected.'))
		}
	}
	
	# fill in grid coordinates
	df[, 1:2] <- sf::st_coordinates(gridTemplate)
		
	# avoid identical column names
	if (anyDuplicated(colnames(df)) > 0) {
		dupNames <- names(which(table(colnames(df)) > 1) == TRUE)
		for (i in 1:length(dupNames)) {
			colnames(df)[which(colnames(df) == dupNames[i])] <- paste0(dupNames[i], '.', 1:length(which(colnames(df) == dupNames[i])))
			
		}
	}
	
	# add in grid cell index if requested
	if (id) {
		
		if (inherits(x[[1]], 'epmGrid')) {
			df_grid <- x[[1]][[1]]
		} else {
			df_grid <- x[[1]]
		}
		
		if (inherits(df_grid, 'sf')) {
			# convert the xy to spatial sf object
			xxPts <- sf::st_geometry(sf::st_as_sf(df, coords = c('x', 'y'), crs = sf::st_crs(df_grid)))
			ptCheck <- sf::st_intersects(xxPts, df_grid)
			df$gridCellID <- unlist(ptCheck)

		} else if (inherits(df_grid, 'SpatRaster')) {
			df$gridCellID <- terra::cellFromXY(df_grid, df[, c('x', 'y')])
		}
		
		df <- df[, c('x', 'y', 'gridCellID', setdiff(colnames(df), c('x', 'y', 'gridCellID')))]
	}	

	return(df)
}














