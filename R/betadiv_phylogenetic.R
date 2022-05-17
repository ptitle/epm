##'@title Map phylogenetic turnover in species communities
##'
##'@description Multisite phylogenetic community dissimilarity is calculated for
##'  each cell within a circular moving window of neighboring cells. To implement a 
##'  custom function, see \code{\link{customBetaDiv}}.
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
##'@details For each cell, multisite dissimilarity is calculated for the focal
##'  cell and its neighbors. If \code{focalCoord} is specified, then instead of
##'  multisite dissimilarity within a moving window of gridcells, pairwise
##'  dissimilarity is calculated from the cell at the focal coordinates, to all
##'  other cells.
##'
##'  All metrics are based on Sorensen dissimilarity and range from 0 to 1: For
##'  each metric, the following components can be specified. These components
##'  are additive, such that the full metric = turnover + nestedness. \itemize{
##'  \item{turnover}: species turnover without the influence of richness
##'  differences \item{nestedness}: species turnover due to differences in
##'  richness \item{full}: the combined turnover due to both differences in
##'  richness and pure turnover }
##'
##'  If the R package spdep is installed, this function should run more quickly.
##'
##'
##'
##'@return Returns a sf polygons object (if hex grid) or a SpatRaster object (if
##'  square grid) with multisite community dissimilarity for each grid cell.
##'
##'@author Pascal Title
##'
##'@references
##'
##'Baselga, A. The relationship between species replacement, dissimilarity
##'derived from nestedness, and nestedness. Global Ecology and Biogeography 21
##'(2012): 1223–1232.
##'
##'Leprieur, F, Albouy, C, De Bortoli, J, Cowman, PF, Bellwood, DR & Mouillot,
##'D. Quantifying Phylogenetic Beta Diversity: Distinguishing between "True"
##'Turnover of Lineages and Phylogenetic Diversity Gradients. PLoS ONE 7 (2012):
##'e42760–12.
##'
##'
##' @examples
##' \donttest{
##' tamiasEPM
##'
##' tamiasEPM <- addPhylo(tamiasEPM, tamiasTree)
##'
##' # phylogenetic turnover
##' beta_phylo_turnover <- betadiv_phylogenetic(tamiasEPM, radius = 70000,
##' 		component = 'turnover')
##' beta_phylo_nestedness <- betadiv_phylogenetic(tamiasEPM, radius = 70000,
##' 		component = 'nestedness')
##' beta_phylo_full <- betadiv_phylogenetic(tamiasEPM, radius = 70000,
##' 		component = 'full')
##'
##' oldpar <- par(mfrow=c(1,3))
##' plot(beta_phylo_turnover, reset = FALSE, key.pos = NULL)
##' plot(beta_phylo_nestedness, reset = FALSE, key.pos = NULL)
##' plot(beta_phylo_full, reset = FALSE, key.pos = NULL)
##'
##' # using square grid epmGrid
##' tamiasEPM2 <- createEPMgrid(tamiasPolyList, resolution = 50000,
##' 	cellType = 'square', method = 'centroid')
##' tamiasEPM2 <- addPhylo(tamiasEPM2, tamiasTree)
##'
##' beta_phylo_full <- betadiv_phylogenetic(tamiasEPM2, radius = 70000,
##' 		component = 'full')
##' beta_phylo_full_slow <- betadiv_phylogenetic(tamiasEPM2, radius = 70000,
##' 		component = 'full', slow = TRUE)
##'
##' par(mfrow = c(1,2))
##' terra::plot(beta_phylo_full, col = sf::sf.colors(100))
##' terra::plot(beta_phylo_full_slow, col = sf::sf.colors(100))
##'
##' # dissimilarity from a focal cell
##' focalBeta <- betadiv_phylogenetic(tamiasEPM, radius = 70000,
##'			component = 'full', focalCoord = c(-1413764, 573610.8))
##' plot(focalBeta, reset = FALSE)
##' points(-1413764, 573610.8, pch = 3, col = 'white')
##'
##'	par(oldpar)
##' }
##'@export

betadiv_phylogenetic <- function(x, radius, component = 'full', focalCoord = NULL, slow = FALSE, nThreads = 1) {
	# radius is in map units
		
	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	if (radius <= attributes(x)$resolution) {
		stop(paste0('The radius must at greater than ', attributes(x)$resolution, '.'))
	}
	
	component <- match.arg(component, choices = c('turnover', 'nestedness', 'full'))
	if (!component %in% c('turnover', 'nestedness', 'full')) {
		stop('component must be turnover, nestedness or full.')
	}
	
	if (inherits(x[['phylo']], 'phylo')) {
		x[['phylo']] <- list(x[['phylo']])
		class(x[['phylo']]) <- 'multiPhylo'
	}
	
	# There appear to be some potential numerical precision issues with defining neighbors.
	# Therefore, if radius is a multiple of resolution, add 1/100 of the resolution. 
	if (radius %% attributes(x)$resolution == 0) {
		radius <- radius + attributes(x)$resolution * 0.01
	}

	# ----------------------------------------------------------
	# Subset appropriately depending on dataset of interest
	
	if (is.null(x[['phylo']])) {
		stop('epmGrid object does not contain a phylo object!')
	}
	
	# prune epmGrid object down to species shared with phylogeny
	if (!identical(sort(x[['geogSpecies']]), sort(x[['phylo']][[1]]$tip.label))) {
		x[['speciesList']] <- intersectList(x[['speciesList']], x[['phylo']][[1]]$tip.label)
	}

	# do calculations for each tree in epmGrid object
		
	retList <- vector('list', length(x[['phylo']]))
	
	for (i in 1:length(x[['phylo']])) {
		
		if (length(x[['phylo']]) > 1) {
			message('\tCalculating metric for tree ', i, '...')
		}
		
		tree <- x[['phylo']][[i]]
	
		# ----------------------------------------------------------
		# if focal coordinate, then we will calculate pairwise dissimilarity from the focal cell to all other cells.
		
		if (!is.null(focalCoord)) {
			
			# convert NA communities to "empty" (easier to handle in Rcpp)
			isNAind <- which(sapply(x[['speciesList']], anyNA) == TRUE)
			if (length(isNAind) > 0) {
				for (j in 1:length(isNAind)) {
					x[['speciesList']][[isNAind[j]]] <- 'empty'
				}
			}	
		
			pairwiseD <- calcPairwisePhylosor2(x[['speciesList']], tree, component = component)
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
		
			spEdges <- ape::nodepath(tree)
			names(spEdges) <- tree$tip.label
	
			# ----------------------------------------------------------
			# Generate list of neighborhoods
			
			if (inherits(x[[1]], 'sf')) {
	
				if (requireNamespace('spdep', quietly = TRUE)) {
					nb <- spdep::dnearneigh(sf::st_centroid(sf::st_geometry(x[[1]])), d1 = 0, d2 = radius)
				} else {
					nb <- sf::st_is_within_distance(sf::st_centroid(sf::st_geometry(x[[1]])), sf::st_centroid(sf::st_geometry(x[[1]])), dist = radius)
					
					# remove focal cell from set of neighbors
					for (j in 1:length(nb)) {
						nb[[j]] <- setdiff(nb[[j]], j)
					}
				}	
				
				# average neighborhood size
				message('\tgridcell neighborhoods: median ', stats::median(lengths(nb)), ', range ', min(lengths(nb)), ' - ', max(lengths(nb)), ' cells')
				# hist(lengths(nb))
						
				# prep for parallel computing
				if (nThreads == 1) {
					cl <- NULL
				} else if (nThreads > 1) {
					if (.Platform$OS.type != 'windows') {
						cl <- nThreads
					} else {
						cl <- parallel::makeCluster(nThreads)
						parallel::clusterExport(cl, c('x', 'nb', 'tree', 'getEdgeIndices', 'calcPhyloMultiSite', 'spEdges', 'component'))
					}
				}	
				
				cellVals <- pbapply::pbsapply(1:nrow(x[[1]]), function(k) {
				# for (k in 1:nrow(x[[1]])) {
			
					cells <- x[['cellCommInd']][c(k, nb[[k]])]
					allsites <- x[['speciesList']][cells]
		
					if (length(unique(allsites)) > 1) {
						# phylo.beta.multi(generateOccurrenceMatrix(x, c(k, nb[[k]])), tree)$phylo.beta.SOR
						calcPhyloMultiSite(allsites, tree, spEdges, component)
					} else {
						0
					}
				
				}, cl = cl)
				
				# identical(cellVals, cellVals1)
				# which((cellVals == cellVals1) == FALSE)
				# range(abs(cellVals - cellVals1), na.rm = TRUE)
				# tail(order(abs(cellVals - cellVals1)))
						
			} else if (inherits(x[[1]], 'SpatRaster')) {
				
				# write two versions, one low memory, one high memory
				
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
						for (j in 1:length(nb)) {
							nb[[j]] <- setdiff(nb[[j]], j)
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
							parallel::clusterExport(cl, c('x', 'datCells', 'tree', 'gridCentroids', 'nb', 'getEdgeIndices', 'calcPhyloMultiSite', 'spEdges', 'component'))
						}
					}	
					
					cellVals <- rep(NA, terra::ncell(x[[1]]))
					cellVals[datCells] <- pbapply::pbsapply(1:nrow(gridCentroids), function(k) {
					# for (k in 1:nrow(x[[1]])) {
		
						nbCells <- sf::st_drop_geometry(gridCentroids[c(k, nb[[k]]), 'cellInd'])[,1]
						nbCells <- x[['cellCommInd']][nbCells]
						allsites <- x[['speciesList']][nbCells]
					
						# message(k)
						# visual verification of focal cell + neighborhood
						# plot(st_geometry(x[[1]])[nbCells,], border = 'white')
						# plot(st_geometry(x[[1]]), add = TRUE, border = 'gray95')
						# plot(st_geometry(x[[1]])[nbCells,], add = TRUE)
						# plot(st_geometry(x[[1]])[focalCell, ], col = 'red', add = TRUE)
						# plot(st_buffer(st_centroid(st_geometry(x[[1]][focalCell,])), dist = radius), add = TRUE, border = 'blue')
						
						if (length(unique(allsites)) > 1) {
							calcPhyloMultiSite(allsites, tree, spEdges, component)
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
							parallel::clusterExport(cl, c('x', 'getEdgeIndices', 'calcPhyloMultiSite', 'spEdges', 'tree', 'component'))
						}
					}	
		
					cellVals <- pbapply::pbsapply(1:terra::ncell(x[[1]]), function(k) {
					# for (k in 1:terra::ncell(x[[1]])) {
						focalCell <- k
						# message(k)
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
								for (j in 1:length(nb)) {
									nb[[j]] <- setdiff(nb[[j]], j)
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
								calcPhyloMultiSite(allsites, tree, spEdges, component)
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
	
		metricName <- paste0('phylogenetic_', component)
		
		if (inherits(x[[1]], 'sf')) {
			ret <- x[[1]]
			ret[metricName] <- cellVals
			ret <- ret[metricName]
		} else {
			ret <- terra::rast(x[[1]][[1]])
			names(ret) <- metricName
			terra::values(ret) <- cellVals
		}
		
		retList[[i]] <- ret
	}
	
	if (length(retList) == 1) {
		retList <- retList[[1]]
	}
	
	return(retList)
}



getEdgeIndices <- function(spEdges, taxa, tree) {
	setdiff(match(unique(unlist(spEdges[taxa])), tree$edge[,2]), NA)
}

	
calcPhyloMultiSite <- function(siteList, tree, spEdges, component) {
	
	# total branch lengths across sites
	sumSi <- sum(tree$edge.length[unlist(lapply(siteList, function(y) getEdgeIndices(spEdges, y, tree)))])
	
	# common branch length across sites
	St <- sum(tree$edge.length[getEdgeIndices(spEdges, unique(unlist(siteList)), tree)])
	
	minDiffs <- matrix(0, nrow = length(siteList), ncol = length(siteList))
	maxDiffs <- minDiffs
	for (i in 1:length(siteList)) {
		for (j in 1:length(siteList)) {
			if (i < j) {
				edges_i <- getEdgeIndices(spEdges, siteList[[i]], tree)
				edges_j <- getEdgeIndices(spEdges, siteList[[j]], tree)
				pd_ij <- sum(tree$edge.length[setdiff(edges_i, edges_j)])
				pd_ji <- sum(tree$edge.length[setdiff(edges_j, edges_i)])
				minDiffs[i,j] <- min(c(pd_ij, pd_ji))
				maxDiffs[i,j] <- max(c(pd_ij, pd_ji))
			}
		}
	}

	minDiffs <- sum(minDiffs)
	maxDiffs <- sum(maxDiffs)
	
	if (component == 'full') {
		(minDiffs + maxDiffs) / (2 * (sumSi - St) + minDiffs + maxDiffs)	
	} else if (component == 'nestedness') {
		minDiffs / ((sumSi - St) + minDiffs)
	} else {
		((minDiffs + maxDiffs) / (2 * (sumSi - St) + minDiffs + maxDiffs)) - (minDiffs / ((sumSi - St) + minDiffs))
	}
}
	

## TESTING ##
# verify distances with betapart package to ensure calculations are correct

# # sourceCpp('~/Dropbox/ecoPhyloMapper/src/betaDiversity.cpp')
# testWithBetaPart <- function(x, metric, component) {
	# if (!metric %in% c('taxonomic', 'phylogenetic') | !component %in% c('turnover', 'nestedness', 'full')) {
		# stop()
	# }
	
	# commMat <- matrix(0, nrow = length(x[['speciesList']]), ncol = length(na.omit(unique(unlist(x[['speciesList']])))))
	# colnames(commMat) <- na.omit(unique(unlist(x[['speciesList']])))
	# for (i in 1:length(x[['speciesList']])) {
		# if (!anyNA(x[['speciesList']][[i]])) {
			# commMat[i,] <- as.numeric(colnames(commMat) %in% x[['speciesList']][[i]])
		# }
	# }
	# if ('empty' %in% colnames(commMat)) {
		# commMat[, 'empty'] <- 0
	# }
	
	# emptyInd <- which(rowSums(commMat) == 0)
	
	# if (metric == 'phylogenetic') {
		# beta.mat <- betapart::phylo.beta.pair(commMat, x[['phylo']])
		# if (component == 'turnover') {
			# beta.mat <- as.matrix(beta.mat$phylo.beta.sim)
		# }
		# if (component == 'nestedness') {
			# beta.mat <- as.matrix(beta.mat$phylo.beta.sne)
		# }
		# if (component == 'full') {
			# beta.mat <- as.matrix(beta.mat$phylo.beta.sor)
		# }
	# } else if (metric == 'taxonomic') {
		# beta.mat <- betapart::beta.pair(commMat)
		# if (component == 'turnover') {
			# beta.mat <- as.matrix(beta.mat$beta.sim)
		# }
		# if (component == 'nestedness') {
			# beta.mat <- as.matrix(beta.mat$beta.sne)
		# }
		# if (component == 'full') {
			# beta.mat <- as.matrix(beta.mat$beta.sor)
		# }			
	# }
	
	# dtest <- matrix(nrow = length(x[['speciesList']]), ncol = length(x[['speciesList']]))
	# for (i in 1:length(x[['speciesList']])) {
		# for (j in 1:length(x[['speciesList']])) {
			# if (i >= j) {
				# dtest[i,j] <- beta.mat[i,j]
				# dtest[j,i] <- dtest[i,j]
			# }
		# }
	# }
	# dtest[is.nan(dtest)] <- NA
	# dtest[emptyInd,] <- NA
	# dtest[, emptyInd] <- NA
	# return(dtest)
# }



# # taxonomic

# pairwiseD <- calcPairwiseTaxonomicSorensen(x[['speciesList']], component = 'turnover'); pairwiseD[pairwiseD < 0] <- NA
# betapartD <- testWithBetaPart(x, metric = 'taxonomic', component = 'turnover')
# pairwiseD[1:10,1:10]
# betapartD[1:10,1:10]
# identical(pairwiseD, betapartD)

# pairwiseD <- calcPairwiseTaxonomicSorensen(x[['speciesList']], component = 'nestedness'); pairwiseD[pairwiseD < 0] <- NA
# betapartD <- testWithBetaPart(x, metric = 'taxonomic', component = 'nestedness')
# pairwiseD[1:10,1:10]
# betapartD[1:10,1:10]
# identical(pairwiseD, betapartD)

# pairwiseD <- calcPairwiseTaxonomicSorensen(x[['speciesList']], component = 'full'); pairwiseD[pairwiseD < 0] <- NA
# betapartD <- testWithBetaPart(x, metric = 'taxonomic', component = 'full')
# pairwiseD[1:10,1:10]
# betapartD[1:10,1:10]
# identical(pairwiseD, betapartD)



# # phylogenetic
# pairwiseD <- calcPairwisePhylosor2(x[['speciesList']], x[['phylo']], component ='turnover'); pairwiseD[pairwiseD < 0] <- NA
# betapartD <- testWithBetaPart(x, metric = 'phylogenetic', component = 'turnover')
# pairwiseD[1:10,1:10]
# betapartD[1:10,1:10]
# identical(pairwiseD, betapartD)
# identical(round(pairwiseD, 6), round(betapartD, 6))

# pairwiseD <- calcPairwisePhylosor(x[['speciesList']], x[['phylo']], component ='nestedness'); pairwiseD[pairwiseD < 0] <- NA
# betapartD <- testWithBetaPart(x, metric = 'phylogenetic', component = 'nestedness')
# pairwiseD[1:10,1:10]
# betapartD[1:10,1:10]
# identical(pairwiseD, betapartD)
# identical(round(pairwiseD, 6), round(betapartD, 6))

# pairwiseD <- calcPairwisePhylosor(x[['speciesList']], x[['phylo']], component ='full'); pairwiseD[pairwiseD < 0] <- NA
# betapartD <- testWithBetaPart(x, metric = 'phylogenetic', component = 'full')
# pairwiseD[1:10,1:10]
# betapartD[1:10,1:10]
# identical(pairwiseD, betapartD)
# identical(round(pairwiseD, 6), round(betapartD, 6))




# test with morphdist

# # morphDist <- function(x) {
	# dtest <- matrix(nrow = length(x[['speciesList']]), ncol = length(x[['speciesList']]))

	# d <- as.matrix(dist(x[['data']]))
	# diag(d) <- NA

	# for (i in 1:length(x[['speciesList']])) {
		# for (j in 1:length(x[['speciesList']])) {
			# comm1 <- x[['speciesList']][[i]]
			# comm2 <- x[['speciesList']][[j]]
			# if (i >= j) {
				# if (!anyNA(c(comm1, comm2)) & !('empty' %in% c(comm1, comm2))) {
					# if (identical(comm1, comm2)) {
						# dtest[i,j] <- 0
						# dtest[j,i] <- 0
					# } else {
						# dtest[i,j] <- mean(d[comm1, comm2], na.rm=TRUE)			
						# dtest[j,i] <- dtest[i,j]
					# }
				# }
			# }
		# }
	# }
	# return(dtest)
# }



