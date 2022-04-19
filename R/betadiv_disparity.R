##' @title Map change in morphological disparity
##'
##' @description Change in morphological disparity is calculating across a
##'   moving window of neighboring grid cells. To implement a custom function,
##' 	  see \code{\link{customBetaDiv}}.
##'
##' @param x object of class \code{epmGrid}.
##' @param radius Radius of the moving window in map units.
##' @param slow if TRUE, use an alternate implementation that has a smaller
##'   memory footprint but that is likely to be much slower. Most useful for
##'   high spatial resolution.
##' @param nThreads number of threads for parallelization

##'
##' @details For each gridcell neighborhood (defined by the radius), we
##' calculate the proportion of the full disparity contained in those grid
##' cells, and then take the standard deviation of those proportions across the
##' gridcell neighborhood. This way, the returned values reflect how much
##' disparity (relative to the overall total disparity) changes across a moving
##' window.
##'
##' If the R package spdep is installed, this function should run more quickly.
##'
##'
##' @return Returns a sf polygons object (if hex grid) or a SpatRaster object
##'   (if square grid).
##'
##' @author Pascal Title
##'
##' @references
##'
##' Foote M. 1993. Contributions of individual taxa to overall morphological
##' disparity. Paleobiology. 19:403–419.
##'
##'
##' @examples
##' \donttest{
##' tamiasEPM
##'
##' tamiasEPM <- addTraits(tamiasEPM, tamiasTraits)
##'
##' z <- betadiv_disparity(tamiasEPM, radius = 150000)
##'
##' plot(z)
##'
##' # using square grid epmGrid
##' tamiasEPM2 <- createEPMgrid(tamiasPolyList, resolution = 50000,
##' 	cellType = 'square', method = 'centroid')
##' tamiasEPM2 <- addTraits(tamiasEPM2, tamiasTraits)
##' z2 <- betadiv_disparity(tamiasEPM2, radius = 150000)
##'
##' terra::plot(z2, col = sf::sf.colors(100))
##'
##' }
##' @export


betadiv_disparity <- function(x, radius, slow = FALSE, nThreads = 1) {
	# x <- tamiasEPM
	# x <- addTraits(x, tamiasTraits)
	# radius = 150000; slow = FALSE; nThreads = 1
	# radius is in map units
		
	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	if (radius <= attributes(x)$resolution) {
		stop(paste0('The radius must at greater than ', attributes(x)$resolution, '.'))
	}
	
	if (nThreads > parallel::detectCores()) {
		stop('nThreads is greater than what is available.')
	}
	
	if (nThreads > 1 & .Platform$OS.type == 'windows') {
		stop('Parallel processing for Windows OS is currently non-functional.')
	}
	
	# There appear to be some potential numerical precision issues with defining neighbors.
	# Therefore, if radius is a multiple of resolution, add 1/100 of the resolution. 
	if (radius %% attributes(x)$resolution == 0) {
		radius <- radius + attributes(x)$resolution * 0.01
	}
	
	# ----------------------------------------------------------
	# Subset appropriately depending on dataset of interest
	
	if (is.null(x[['data']])) {
		stop('epmGrid object does not contain any associated trait data!')
	}
	
	if (is.vector(x[['data']])) {
		stop('This function is intended for multivariate shape data.')
	}
	
	# prune epmGrid object down to species shared with data
	if (is.vector(x[['data']])) {
		if (!identical(sort(x[['geogSpecies']]), sort(names(x[['data']])))) {
			x[['speciesList']] <- intersectList(x[['speciesList']], names(x[['data']]))
		}
	} else {
		if (!identical(sort(x[['geogSpecies']]), sort(rownames(x[['data']])))) {
			x[['speciesList']] <- intersectList(x[['speciesList']], rownames(x[['data']]))
		}
	}
		
	
	# ----------------------------------------------------------
	# Calculate partial disparity of all species in dataset
	partialDisp <- getSpPartialDisparities(x[['data']])
	
	# Here, total disparity = sum(partialDisp)
	
	# for each cell neighborhood, we will calculate the proportion of the full disparity contained in those grid cells, and then take the standard deviation. This way, the returned value reflects how much relative disparity changes across a moving window.
	
	if (inherits(x[[1]], 'sf')) {
		
		# Generate list of neighborhoods

		if (requireNamespace('spdep', quietly = TRUE)) {
			nb <- spdep::dnearneigh(sf::st_centroid(sf::st_geometry(x[[1]])), d1 = 0, d2 = radius)
		} else {
			nb <- sf::st_is_within_distance(sf::st_centroid(sf::st_geometry(x[[1]])), sf::st_centroid(sf::st_geometry(x[[1]])), dist = radius)
			
			# remove focal cell from set of neighbors
			for (i in 1:length(nb)) {
				nb[[i]] <- setdiff(nb[[i]], i)
			}
		}	

		# average neighborhood size
		message('\tgridcell neighborhoods: median ', stats::median(lengths(nb)), ', range ', min(lengths(nb)), ' - ', max(lengths(nb)))
		# hist(lengths(nb))
		
		
		# prep for parallel computing
		if (nThreads == 1) {
			cl <- NULL
		} else if (nThreads > 1) {
			if (.Platform$OS.type != 'windows') {
				cl <- nThreads
			} else {
				cl <- parallel::makeCluster(nThreads)
				parallel::clusterExport(cl, c('x', 'nb', 'getSpPartialDisparities'))
			}
		}	
				
		cellVals <- pbapply::pblapply(1:nrow(x[[1]]), function(k) {
			
			focalCell <- k
			nbCells <- nb[[k]]
			
			# # visual verification of focal cell + neighborhood
			# plot(st_buffer(st_centroid(st_geometry(x[[1]][focalCell,])), dist = radius * 1.5), border = 'white')
			# plot(st_geometry(x[[1]]), add = TRUE, border = 'gray95')
			# plot(st_geometry(x[[1]])[nbCells,], add = TRUE)
			# plot(st_geometry(x[[1]])[focalCell, ], col = 'red', add = TRUE)
			# plot(st_buffer(st_centroid(st_geometry(x[[1]][focalCell,])), dist = radius), add = TRUE, border = 'blue')
	
			# convert to cell indices
			focalCell <- x[['cellCommInd']][focalCell]
			nbCells <- x[['cellCommInd']][nbCells]
					
			# get species communities
			sites <- x[['speciesList']][c(focalCell, nbCells)]
						
			return(sd(sapply(sites, function(y) sum(partialDisp[y])) / sum(partialDisp)))
		}, cl = cl)
		
	} else if (inherits(x[[1]], 'SpatRaster')) {
		
		if (!slow) {
			
			# Generate list of neighborhoods
			datCells <- which(terra::values(x[[1]]['spRichness']) > 0)
			gridCentroids <- terra::xyFromCell(x[[1]], cell = datCells)
			gridCentroids <- sf::st_as_sf(as.data.frame(gridCentroids), coords = 1:2, crs = attributes(x)$crs)
			gridCentroids$cellInd <- datCells

			if (requireNamespace('spdep', quietly = TRUE)) {
				nb <- spdep::dnearneigh(gridCentroids, d1 = 0, d2 = radius)
			} else {
				nb <- sf::st_is_within_distance(gridCentroids, gridCentroids, dist = radius)
				
				# remove focal cell from set of neighbors
				for (i in 1:length(nb)) {
					nb[[i]] <- setdiff(nb[[i]], i)
				}
			}	
			
			# average neighborhood size
			message('\tgridcell neighborhoods: median ', stats::median(lengths(nb)), ', range ', min(lengths(nb)), ' - ', max(lengths(nb)), ' cells')
			
			if (nThreads == 1) {
				cl <- NULL
			} else if (nThreads > 1) {
				if (.Platform$OS.type != 'windows') {
					cl <- nThreads
				} else {
					cl <- parallel::makeCluster(nThreads)
					parallel::clusterExport(cl, c('x', 'datCells', 'nb', 'getSpPartialDisparities'))
				}
			}			

			cellVals <- rep(NA, terra::ncell(x[[1]]))
			cellVals[datCells] <- pbapply::pblapply(1:nrow(gridCentroids), function(k) {
				
				nbCells <- sf::st_drop_geometry(gridCentroids[c(k, nb[[k]]), 'cellInd'])[,1]
				nbCells <- x[['cellCommInd']][nbCells]
				sites <- x[['speciesList']][nbCells]
				
				# # visual verification of focal cell + neighborhood
				# plot(gridCentroids[nb[[i]],])
				# plot(gridCentroids[i,], pch = 3, col = 'red', cex = 2, add = TRUE)
				
				# plot(xyFromCell(x[[1]], nbCells))
				# points(xyFromCell(x[[1]], focalCell), col = 'red', cex = 2, pch = 3)
												
							
				return(sd(sapply(sites, function(y) sum(partialDisp[y])) / sum(partialDisp)))
			}, cl = cl)
			
		} else {
					
			# prep for parallel computing
			## makeCluster seems to have issues with SpatRaster
			
			if (nThreads == 1) {
				cl <- NULL
			} else if (nThreads > 1) {
				if (.Platform$OS.type != 'windows') {
					cl <- nThreads
				} else {
					cl <- parallel::makeCluster(nThreads)
					parallel::clusterExport(cl, c('x', 'getSpPartialDisparities'))
				}
			}			
		
			cellVals <- pbapply::pblapply(1:terra::ncell(x[[1]]), function(k) {
				
				focalCell <- k
				if (!anyNA((x[['speciesList']][[x[['cellCommInd']][focalCell]]]))) {
					
					# get neighborhood for focal cell
					focalxy <- sf::st_as_sf(as.data.frame(terra::xyFromCell(x[[1]], focalCell)), coords = 1:2, crs = terra::crs(x[[1]]))
					focalCircle <- sf::st_buffer(focalxy, dist = radius * 2)
					nbCells <- terra::cells(x[[1]], terra::vect(focalCircle))[, 'cell']
					gridCentroids <- sf::st_as_sf(as.data.frame(terra::xyFromCell(x[[1]], nbCells)), coords = 1:2, crs = terra::crs(x[[1]]))

					if (requireNamespace('spdep', quietly = TRUE)) {
						nb <- spdep::dnearneigh(gridCentroids, d1 = 0, d2 = radius)
					} else {
						nb <- sf::st_is_within_distance(gridCentroids, gridCentroids, dist = radius)
						
						# remove focal cell from set of neighbors
						for (i in 1:length(nb)) {
							nb[[i]] <- setdiff(nb[[i]], i)
						}
					}	

					nbCells <- nbCells[nb[[which(nbCells == focalCell)]]]
					
					# without dnearneigh step					
					# focalxy <- sf::st_as_sf(as.data.frame(terra::xyFromCell(x[[1]], focalCell)), coords = 1:2, crs = terra::crs(x[[1]]))
					# focalCircle <- sf::st_buffer(focalxy, dist = radius)
					# nbCells <- setdiff(terra::cells(x[[1]], terra::vect(focalCircle))[, 'cell'], focalCell)
			
					# convert to cell indices
					focalCell <- x[['cellCommInd']][focalCell]
					nbCells <- x[['cellCommInd']][nbCells]
							
					# get species communities
					sites <- x[['speciesList']][c(focalCell, nbCells)]
					sites <- sites[!sapply(sites, anyNA)]
						
					return(sd(sapply(sites, function(y) sum(partialDisp[y])) / sum(partialDisp)))

				} else {
					return(NA)
				}
			}, cl = cl)
		}
	} else {
		stop('Grid format not recognized.')
	}
	
	if (nThreads > 1 & .Platform$OS.type == 'windows') parallel::stopCluster(cl)
	
	cellVals <- unlist(cellVals)

	
	metricName <- paste0('disparityTurnover')
	
	if (inherits(x[[1]], 'sf')) {	
		ret <- x[[1]]
		ret[metricName] <- cellVals
		ret <- ret[metricName]

	} else {
		ret <- terra::rast(x[[1]][[1]])
		names(ret) <- metricName
		terra::values(ret) <- cellVals		
	}
	return(ret)
}





		

