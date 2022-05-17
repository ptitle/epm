##'@title Map turnover in species communities
##'
##'@description Multisite taxonomic community dissimilarity is calculated for
##'  each cell within a circular moving window of neighboring cells. To implement 
##' a custom function, see \code{\link{customBetaDiv}}.
##'
##'@param x object of class \code{epmGrid}.
##'@param radius Radius of the moving window in map units.
##'@param component which component of beta diversity to use, can be
##'  \code{"turnover"}, \code{"nestedness"} or \code{"full"}
##'@param focalCoord vector of x and y coordinate, see details
##'@param slow if TRUE, use an alternate implementation that has a smaller
##'  memory footprint but that is likely to be much slower. Most useful for high
##'  spatial resolution.
##'@param nThreads number of threads for parallelization
##'
##'@details For each cell, multisite dissimilarity is calculated from the focal
##'cell and its neighbors. If \code{focalCoord} is specified, then instead of
##'multisite dissimilarity within a moving window of gridcells, pairwise
##'dissimilarity is calculated from the cell at the focal coordinates, to all
##'other cells.
##'
##'All metrics are based on Sorensen dissimilarity and range from 0 to 1. \cr
##'For each metric, the following components can be specified. These components
##'are additive, such that the full metric = turnover + nestedness. 
##'\itemize{
##'    \item{turnover}: species turnover without the influence of richness
##'differences 
##'    \item{nestedness}: species turnover due to differences in
##'richness 
##     \item{full}: the combined turnover due to both differences in
##'richness and pure turnover 
##'}
##'
##'If the R package spdep is installed, this function should run more quickly.
##'
##'
##'
##'@return Returns a grid with multi-site community dissimilarity for each cell.
##'
##'@author Pascal Title
##'
##'@references
##'
##'Baselga, A. The relationship between species replacement, dissimilarity
##'derived from nestedness, and nestedness. Global Ecology and Biogeography 21
##'(2012): 1223–1232.
##'
##'
##' @examples
##' \donttest{
##' tamiasEPM
##'
##' tamiasEPM <- addPhylo(tamiasEPM, tamiasTree)
##' tamiasEPM <- addTraits(tamiasEPM, tamiasTraits)
##'
##' # taxonomic turnover
##' beta_taxonomic_turnover <- betadiv_taxonomic(tamiasEPM, radius = 70000,
##' 		component = 'turnover')
##' beta_taxonomic_nestedness <- betadiv_taxonomic(tamiasEPM, radius = 70000,
##'			component = 'nestedness')
##' beta_taxonomic_full <- betadiv_taxonomic(tamiasEPM, radius = 70000,
##'			component = 'full')
##'
##'
##'
##' oldpar <- par(mfrow = c(1, 3))
##' plot(beta_taxonomic_turnover, reset = FALSE, key.pos = NULL)
##' plot(beta_taxonomic_nestedness, reset = FALSE, key.pos = NULL)
##' plot(beta_taxonomic_full, reset = FALSE, key.pos = NULL)
##'
##'
##'
##' # using square grid epmGrid
##' tamiasEPM2 <- createEPMgrid(tamiasPolyList, resolution = 50000,
##' 	cellType = 'square', method = 'centroid')
##'
##' beta_taxonomic_full <- betadiv_taxonomic(tamiasEPM2, radius = 70000,
##' 		component = 'full')
##' beta_taxonomic_full_slow <- betadiv_taxonomic(tamiasEPM2, radius = 70000,
##' 		component = 'full', slow = TRUE)
##'
##' par(mfrow=c(1,2))
##' terra::plot(beta_taxonomic_full, col = sf::sf.colors(100))
##' terra::plot(beta_taxonomic_full_slow, col = sf::sf.colors(100))
##'
##' # dissimilarity from a focal cell
##' focalBeta <- betadiv_taxonomic(tamiasEPM, radius = 70000,
##'			component = 'full', focalCoord = c(-1413764, 573610.8))
##' plot(focalBeta, reset = FALSE)
##' points(-1413764, 573610.8, pch = 3, col = 'white')
##'
##' par(oldpar)
##' }
##'@export


betadiv_taxonomic <- function(x, radius, component = 'full', focalCoord = NULL, slow = FALSE, nThreads = 1) {
	# radius is in map units
	
	# x <- tamiasEPM; radius <- 70000; component = 'full'; slow = FALSE; nThreads = 1; focalCoord = NULL
	# focalCoord <- c(-1413764, 573610.8)
	
	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	if (radius <= attributes(x)$resolution) {
		stop(paste0('The radius must at least be greater than ', attributes(x)$resolution, '.'))
	}
	
	component <- match.arg(component, choices = c('turnover', 'nestedness', 'full'))
	if (!component %in% c('turnover', 'nestedness', 'full')) {
		stop('component must be turnover, nestedness or full.')
	}
	
	# There appear to be some potential numerical precision issues with defining neighbors.
	# Therefore, if radius is a multiple of resolution, add 1/100 of the resolution. 
	if (radius %% attributes(x)$resolution == 0) {
		radius <- radius + attributes(x)$resolution * 0.01
	}
	
	# ----------------------------------------------------------
	# if focal coordinate, then we will calculate pairwise dissimilarity from the focal cell to all other cells.
	
	if (!is.null(focalCoord)) {
				
		# convert NA communities to "empty" (easier to handle in Rcpp)
		isNAind <- which(sapply(x[['speciesList']], anyNA) == TRUE)
		if (length(isNAind) > 0) {
			for (i in 1:length(isNAind)) {
				x[['speciesList']][[isNAind[i]]] <- 'empty'
			}
		}	
	
		pairwiseD <- calcPairwiseTaxonomicSorensen(x[['speciesList']], component = component)
		pairwiseD[pairwiseD < 0] <- NA
		
		# which community does the focal coordinate represent?
		focalSp <- extractFromEpmGrid(x, focalCoord)
		focalInd <- which(sapply(x[['speciesList']], function(y) identical(focalSp, y)) == TRUE)
		if (length(focalInd) != 1) stop('focalInd != 1.')
		
		if (inherits(x[[1]], 'sf')) {
			cellVals <- pbapply::pbsapply(1:nrow(x[[1]]), function(k) pairwiseD[focalInd, x[['cellCommInd']][k]])
		} else if (inherits(x[[1]], 'SpatRaster')) {
			cellVals <- pbapply::pbsapply(1:terra::ncell(x[[1]]), function(k) pairwiseD[focalInd, x[['cellCommInd']][k]])
		} else {
			stop('epm format not recognized.')
		}
		
	} else {	
		# ----------------------------------------------------------
		# Generate list of neighborhoods
	
		if (inherits(x[[1]], 'sf')) {
			
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
			message('\tgridcell neighborhoods: median ', stats::median(lengths(nb)), ', range ', min(lengths(nb)), ' - ', max(lengths(nb)), ' cells')
			
			# prep for parallel computing
			if (nThreads == 1) {
				cl <- NULL
			} else if (nThreads > 1) {
				if (.Platform$OS.type != 'windows') {
					cl <- nThreads
				} else {
					cl <- parallel::makeCluster(nThreads)
					parallel::clusterExport(cl, c('x', 'nb', 'pairwiseD'))
				}
			}	
			
			cellVals <- pbapply::pbsapply(1:nrow(x[[1]]), function(k) {
			
				cells <- x[['cellCommInd']][c(k, nb[[k]])]
				allsites <- x[['speciesList']][cells]
				
				# visual verification of focal cell + neighborhood
				# plot(sf::st_geometry(x[[1]])[nb[[k]]], border = NA)
				# plot(sf::st_geometry(x[[1]]), add = TRUE, border = 'gray95')
				# plot(sf::st_geometry(x[[1]])[nb[[k]]], add = TRUE)
				# plot(sf::st_geometry(x[[1]])[k], col = 'red', add = TRUE)
				# plot(sf::st_buffer(sf::st_centroid(sf::st_geometry(x[[1]])[k]), dist = radius), add = TRUE, border = 'blue')
				
				# return(beta.multi(generateOccurrenceMatrix(x, c(k, nb[[k]])), index.family = 'sor')$beta.SOR)
				if (length(unique(allsites)) > 1) {
					calcTaxonomicMultisiteBeta(allsites, component)
				} else {
					0
				}
			}, cl = cl)
				
			# identical(cellVals, cellVals2)
			# which((cellVals == cellVals2) == FALSE)
						
		} else if (inherits(x[[1]], 'SpatRaster')) {
	
			if (!slow) {
				
				# convert grid cells to points, and get neighborhoods
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
				
				message('\tgridcell neighborhoods: median ', stats::median(lengths(nb)), ', range ', min(lengths(nb)), ' - ', max(lengths(nb)), ' cells')
				
				# prep for parallel computing
				if (nThreads == 1) {
					cl <- NULL
				} else if (nThreads > 1) {
					if (.Platform$OS.type != 'windows') {
						cl <- nThreads
					} else {
						cl <- parallel::makeCluster(nThreads)
						parallel::clusterExport(cl, c('x', 'nb', 'gridCentroids', 'datCells', 'calcTaxonomicMultisiteBeta', 'component'))
					}
				}	
				
				cellVals <- rep(NA, terra::ncell(x[[1]]))
				cellVals[datCells] <- pbapply::pbsapply(1:nrow(gridCentroids), function(k) {
				# for (i in 1:nrow(x[[1]])) {
			
					# nbCells <- sf::st_drop_geometry(gridCentroids[c(datCells[k], nb[[k]]), 'cellInd'])[,1]
					nbCells <- sf::st_drop_geometry(gridCentroids[c(k, nb[[k]]), 'cellInd'])[,1]
					nbCells <- x[['cellCommInd']][nbCells]
					allsites <- x[['speciesList']][nbCells]
					
					if (length(unique(allsites)) > 1) {
						calcTaxonomicMultisiteBeta(allsites, component)
					} else {
						0
					}
				}, cl = cl)
					
			} else {
	
				# prep for parallel computing
				if (nThreads == 1) {
					cl <- NULL
				} else if (nThreads > 1) {
					if (.Platform$OS.type != 'windows') {
						cl <- nThreads
					} else {
						cl <- parallel::makeCluster(nThreads)
						parallel::clusterExport(cl, c('x', 'radius', 'pairwiseD'))
					}
				}	
				
				cellVals <- pbapply::pbsapply(1:terra::ncell(x[[1]]), function(k) {
				# for (k in 1:terra::ncell(x[[1]])) {
					
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
						
						allsites <- x[['speciesList']][c(focalCell, nbCells)]
						allsites <- allsites[!sapply(allsites, anyNA)]
						
						if (length(unique(allsites)) > 1) {
							calcTaxonomicMultisiteBeta(allsites, component)
						} else {
							0
						}						
					} else {
						NA
					}
				}, cl = cl)
			}		
		} else {
			stop('Grid format not recognized.')
		}
		
		if (nThreads > 1 & .Platform$OS.type == 'windows') parallel::stopCluster(cl)

	}
	
	cellVals <- unlist(cellVals)
	cellVals[is.nan(cellVals)] <- NA 

	metricName <- paste0('taxonomic_', component)

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


calcTaxonomicMultisiteBeta <- function(siteList, component) {
	
	minDiffs <- matrix(0, nrow = length(siteList), ncol = length(siteList))
	maxDiffs <- minDiffs
	for (i in 1:length(siteList)) {
		for (j in 1:length(siteList)) {
			if (i < j) {
				# message(i, ' -- ', j)
				diff1 <- length(setdiff(siteList[[i]], siteList[[j]]))
				diff2 <- length(setdiff(siteList[[j]], siteList[[i]]))
				minDiffs[i,j] <- min(c(diff1, diff2))
				maxDiffs[i,j] <- max(c(diff1, diff2))
			}
		}
	}
		
	minDiffs <- sum(minDiffs)
	maxDiffs <- sum(maxDiffs)
		
	if (component == 'full') {
		(minDiffs + maxDiffs) / (2 * (sum(lengths(siteList)) - length(unique(unlist(siteList)))) + minDiffs + maxDiffs)
	} else if (component == 'nestedness') {
		minDiffs / ((sum(lengths(siteList)) - length(unique(unlist(siteList)))) + minDiffs)
	} else {
		((minDiffs + maxDiffs) / (2 * (sum(lengths(siteList)) - length(unique(unlist(siteList)))) + minDiffs + maxDiffs)) - (minDiffs / ((sum(lengths(siteList)) - length(unique(unlist(siteList)))) + minDiffs))
	}
}
