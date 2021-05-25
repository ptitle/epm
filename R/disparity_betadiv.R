##' @title Map change in morphological disparity
##'
##' @description Change in morphological partial disparity is calculating across a moving window
##' of neighboring grid cells. 
##' 
##' @param x object of class \code{epmGrid}.
##' @param radius Radius of the moving window in map units.
##' @param slow if TRUE, use an alternate implementation that has a lower memory footprint 
##' 	but that is likely to be much slower. Most useful for high spatial resolution.
##' @param nThreads number of threads for parallelization

##'
##' @details
##' 	Disparity is calculated for each cell + its neighboring cells. Partial disparity is then calculated
##' 	for each species, and mean partial disparity of species in neighboring cells (but not the focal cell)
##' 	is calculated and assigned to the focal cell. The value at each grid cell is therefore a measure of how 
##' 	much neighbors contribute to the overall local disparity. 
##' 
##'
##' 
##' @return Returns a sf polygons object with mean partial disparity.
##' 
##' @author Pascal Title
##'
##' @references
##' 
##' Foote M. 1993. Contributions of individual taxa to overall morphological disparity. Paleobiology. 19:403–419.
##' 
##' 
##' @examples
##' \donttest{
##' tamiasEPM
##'
##' tamiasEPM <- addTraits(tamiasEPM, tamiasTraits)
##'
##' z <- disparity_betadiv(tamiasEPM, radius = 150000)
##' 
##' plot(z)
##' 
##' # using square grid epmGrid
##' tamiasEPM2 <- createEPMgrid(tamiasPolyList, resolution = 50000, 
##' 	cellType = 'square', method = 'centroid')
##' tamiasEPM2 <- addTraits(tamiasEPM2, tamiasTraits)
##' z2 <- disparity_betadiv(tamiasEPM2, radius = 150000)
##'
##' terra::plot(z2, col = sf::sf.colors(100), range = c(0,1))
##' 
##' }
##' @export


disparity_betadiv <- function(x, radius, slow = FALSE, nThreads = 1) {
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
	# Here, total disparity = disparity of species in focal cell + neighboring cells.
	# For each cell, calculate the summed partial disparity of the species in neighboring cells (excluding the focal cell)
	# take ratio of summed neighbor partial disparities to total disparity

	if (inherits(x[[1]], 'sf')) {
		
		# Generate list of neighborhoods

		# nb <- sf::st_is_within_distance(sf::st_centroid(sf::st_geometry(x[[1]])), x[[1]], dist = radius)
		nb <- spdep::dnearneigh(sf::st_centroid(sf::st_geometry(x[[1]])), d1 = 0, d2 = radius)

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
				
		cellVals <- pbapply::pblapply(1:nrow(x[[1]]), function(i) {
			
			focalCell <- i
			nbCells <- nb[[i]]
			
			# # visual verification of focal cell + neighborhood
			# plot(st_buffer(st_centroid(st_geometry(x[[1]][focalCell,])), dist = radius * 1.5), border = 'white')
			# plot(st_geometry(x[[1]]), add = TRUE, border = 'gray95')
			# plot(st_geometry(x[[1]])[nbCells,], add = TRUE)
			# plot(st_geometry(x[[1]])[focalCell, ], col = 'red', add = TRUE)
			# plot(st_buffer(st_centroid(st_geometry(x[[1]][focalCell,])), dist = radius), add = TRUE, border = 'blue')
	
			# convert to cell indices
			focalCell <- x[['cellCommInd']][focalCell]
			nbCells <- x[['cellCommInd']][nbCells]
					
			# get species
			focalSp <- x[['speciesList']][[focalCell]]
			nbSp <- unique(unlist(x[['speciesList']][nbCells]))
			
			pd <- getSpPartialDisparities(x[['data']][union(focalSp, nbSp), ])
			# sum(getSpPartialDisparities(x[['data']][union(focalSp, nbSp), ])) - sum(diag(cov(x[['data']][union(focalSp, nbSp), ])))
			
			if (all(is.nan(pd))) {
				return(0)
			} else {
				return(abs((sum(pd[nbSp]) / sum(pd)) - (sum(pd[focalSp]) / sum(pd))))
			}
		}, cl = cl)
		
	} else if (inherits(x[[1]], 'SpatRaster')) {
		
		if (!slow) {
			
			# Generate list of neighborhoods
			datCells <- which(terra::values(x[[1]]['spRichness']) > 0)
			gridCentroids <- terra::xyFromCell(x[[1]], cell = datCells)
			gridCentroids <- sf::st_as_sf(as.data.frame(gridCentroids), coords = 1:2, crs = attributes(x)$crs)
			gridCentroids$cellInd <- datCells
			nb <- spdep::dnearneigh(gridCentroids, d1 = 0, d2 = radius)
			
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

			cellVals <- pbapply::pblapply(1:nrow(gridCentroids), function(i) {
				
				focalCell <- datCells[i]
				nbCells <- nb[[i]]
				nbCells <- sf::st_drop_geometry(gridCentroids[nbCells, 'cellInd'])[,1]
				
				# # visual verification of focal cell + neighborhood
				# plot(gridCentroids[nb[[i]],])
				# plot(gridCentroids[i,], pch = 3, col = 'red', cex = 2, add = TRUE)
				
				# plot(xyFromCell(x[[1]], nbCells))
				# points(xyFromCell(x[[1]], focalCell), col = 'red', cex = 2, pch = 3)
						
				# convert to cell indices
				focalCell <- x[['cellCommInd']][focalCell]
				nbCells <- x[['cellCommInd']][nbCells]
						
				# get species
				focalSp <- x[['speciesList']][[focalCell]]
				nbSp <- unique(unlist(x[['speciesList']][nbCells]))
				
				pd <- getSpPartialDisparities(x[['data']][union(focalSp, nbSp), ])
				# sum(getSpPartialDisparities(x[['data']][union(focalSp, nbSp), ])) - sum(diag(cov(x[['data']][union(focalSp, nbSp), ])))
				
				if (all(is.nan(pd))) {
					return(0)
				} else {
					return(abs((sum(pd[nbSp]) / sum(pd)) - (sum(pd[focalSp]) / sum(pd))))
				}
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
		
			cellVals <- pbapply::pblapply(1:terra::ncell(x[[1]]), function(i) {
				
				focalCell <- i
				if (!anyNA((x[['speciesList']][[x[['cellCommInd']][focalCell]]]))) {
					
					# get neighborhood for focal cell
					focalxy <- sf::st_as_sf(as.data.frame(terra::xyFromCell(x[[1]], focalCell)), coords = 1:2, crs = terra::crs(x[[1]]))
					focalCircle <- sf::st_buffer(focalxy, dist = radius * 2)
					nbCells <- terra::cells(x[[1]], terra::vect(focalCircle))[, 'cell']
					gridCentroids <- sf::st_as_sf(as.data.frame(terra::xyFromCell(x[[1]], nbCells)), coords = 1:2, crs = terra::crs(x[[1]]))
					nb <- spdep::dnearneigh(gridCentroids, d1 = 0, d2 = radius)
					nbCells <- nbCells[nb[[which(nbCells == focalCell)]]]
					
					# without dnearneigh step					
					# focalxy <- sf::st_as_sf(as.data.frame(terra::xyFromCell(x[[1]], focalCell)), coords = 1:2, crs = terra::crs(x[[1]]))
					# focalCircle <- sf::st_buffer(focalxy, dist = radius)
					# nbCells <- setdiff(terra::cells(x[[1]], terra::vect(focalCircle))[, 'cell'], focalCell)
			
					# convert to cell indices
					focalCell <- x[['cellCommInd']][focalCell]
					nbCells <- x[['cellCommInd']][nbCells]
							
					# get species
					focalSp <- x[['speciesList']][[focalCell]]
					nbSp <- unique(unlist(x[['speciesList']][nbCells]))
					nbSp <- nbSp[!is.na(nbSp)]
					
					
					pd <- getSpPartialDisparities(x[['data']][union(focalSp, nbSp), ])
					# sum(getSpPartialDisparities(x[['data']][union(focalSp, nbSp), ])) - sum(diag(cov(x[['data']][union(focalSp, nbSp), ])))
					
					if (all(is.nan(pd))) {
						return(0)
					} else {
						return(abs((sum(pd[nbSp]) / sum(pd)) - (sum(pd[focalSp]) / sum(pd))))
					}
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

	
	# values > 0 indicate that nb taxa contribute unique aspects to disparity, 
	# values < 0 indicate that focal taxa contribute unique aspects to disparity
	# values = 0 indicate that contributions of nb taxa and focal taxa are the same. 
	# absolute val version: ranges 0 to 1:
		# 0 = same contribution to disparity in nb vs focal taxa
		# as values get bigger, greater and greater difference in contribution to total disparity
		# a value of 0.95 would mean that nb taxa contribute 95% of the disparity, whereas the focal taxon contributes 5%, or vice versa. 
	
	metricName <- paste0('disparityTurnover')
	
	if (inherits(x[[1]], 'sf')) {	
		ret <- x[[1]]
		ret[metricName] <- cellVals
		ret <- ret[metricName]

	} else {
		ret <- terra::rast(x[[1]][[1]])
		names(ret) <- metricName
		if (slow) {
			terra::values(ret) <- cellVals		
		} else {
			ret[which(terra::values(x[[1]]['spRichness']) > 0)] <- cellVals	
		}
	}
	return(ret)
}


