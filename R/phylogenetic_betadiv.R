##' @title Map turnover in species communities
##'
##' @description Mean phylogenetic community dissimilarity is 
##' 	calculated for each cell within a circular moving window of neighboring cells. 
##' 
##' @param x object of class \code{epmGrid}.
##' @param radius Radius of the moving window in map units.
##' @param component which component of beta diversity to use, can be \code{"turnover"}, 
##' 	\code{"nestedness"} or \code{"full"}
##' @param slow if TRUE, use an alternate implementation that has a lower memory footprint 
##' 	but that is likely to be much slower. Most useful for high spatial resolution.
##' @param nThreads number of threads for parallelization
##'
##' @details
##' 	For each cell, mean dissimilarity is calculated from the focal cell to each of its neighbors.
##' 
##' 	All metrics are based on Sorensen dissimilarity and range from 0 to 1:
##'		For each metric, the following components can be specified. These components are additive,
##' 	such that the full metric = turnover + nestedness. 
##' 	\itemize{
##'			\item{turnover}: species turnover without the influence of richness differences
##'			\item{nestedness}: species turnover due to differences in richness
##'			\item{full}: the combined turnover due to both differences in richness and pure turnover
##'		}
##'
##'
##'
##' 
##' 
##' @return Returns a sf polygons object with mean community dissimilarity for each grid cell.
##' 
##' @author Pascal Title
##'
##' @references
##' 
##' Baselga, A. The relationship between species replacement, dissimilarity derived from nestedness, 
##' and nestedness. Global Ecology and Biogeography 21 (2012): 1223–1232.
##' 
##' Laffan, SW, et al. Range-weighted metrics of species and phylogenetic turnover can better 
##' resolve biogeographic transition zones. Methods in Ecology and Evolution 7 (2016): 580-588.
##'
##' Leprieur, F, Albouy, C, De Bortoli, J, Cowman, PF, Bellwood, DR & Mouillot, D. Quantifying 
##' Phylogenetic Beta Diversity: Distinguishing between "True" Turnover of Lineages and Phylogenetic 
##' Diversity Gradients. PLoS ONE 7 (2012): e42760–12.
##'
##' Rosauer, D, Laffan, SW, Crisp, MD, Donnellan, SC, Cook, LG. Phylogenetic endemism: a new approach 
##' for identifying geographical concentrations of evolutionary history. Molecular Ecology 
##' 18 (2009): 4061-4072.
##' 
##' @examples
##' \donttest{
##' tamiasEPM
##'
##' tamiasEPM <- addPhylo(tamiasEPM, tamiasTree)
##'
##' # phylogenetic turnover
##' beta_phylo_turnover <- phylogenetic_betadiv(tamiasEPM, radius = 70000, 
##' 		component = 'turnover')
##' beta_phylo_nestedness <- phylogenetic_betadiv(tamiasEPM, radius = 70000, 
##' 		component = 'nestedness')
##' beta_phylo_full <- phylogenetic_betadiv(tamiasEPM, radius = 70000, 
##' 		component = 'full')
##' 
##' par(mfrow=c(1,3))
##' plot(beta_phylo_turnover, reset = FALSE, key.pos = NULL)
##' plot(beta_phylo_nestedness, reset = FALSE, key.pos = NULL)
##' plot(beta_phylo_full, reset = FALSE, key.pos = NULL)
##'
##' # using square grid epmGrid
##' tamiasEPM2 <- createEPMgrid(tamiasPolyList, resolution = 50000, 
##' 	cellType = 'square', method = 'centroid')
##' tamiasEPM2 <- addPhylo(tamiasEPM2, tamiasTree)
##' 
##' beta_phylo_full <- phylogenetic_betadiv(tamiasEPM2, radius = 70000, 
##' 		component = 'full')
##' beta_phylo_full_slow <- phylogenetic_betadiv(tamiasEPM2, radius = 70000, 
##' 		component = 'full', slow = TRUE)
##'
##' par(mfrow=c(1,2))
##' terra::plot(beta_phylo_full, col = sf::sf.colors(100))
##' terra::plot(beta_phylo_full_slow, col = sf::sf.colors(100))
##' 
##'
##' 
##' }
##' @export

phylogenetic_betadiv <- function(x, radius, component = 'full', slow = FALSE, nThreads = 1) {
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
	
	# ----------------------------------------------------------
	# Subset appropriately depending on dataset of interest
	
	if (is.null(x[['phylo']])) {
		stop('epmGrid object does not contain a phylo object!')
	}
	
	# prune epmGrid object down to species shared with phylogeny
	if (!identical(sort(x[['geogSpecies']]), sort(x[['phylo']]$tip.label))) {
		x[['speciesList']] <- intersectList(x[['speciesList']], x[['phylo']]$tip.label)
	}
	
	
	# convert NA communities to "empty" (easier to handle in Rcpp)
	isNAind <- which(sapply(x[['speciesList']], anyNA) == TRUE)
	if (length(isNAind) > 0) {
		for (i in 1:length(isNAind)) {
			x[['speciesList']][[isNAind[i]]] <- 'empty'
		}
	}	
		
	pairwiseD <- calcPairwisePhylosor2(x[['speciesList']], x[['phylo']], component = component)

	pairwiseD[pairwiseD < 0] <- NA
	
	
	# ----------------------------------------------------------
	# Generate list of neighborhoods
	
	if (inherits(x[[1]], 'sf')) {
		# nb <- sf::st_is_within_distance(sf::st_centroid(sf::st_geometry(x[[1]])), x[[1]], dist = radius)
		nb <- spdep::dnearneigh(sf::st_centroid(sf::st_geometry(x[[1]])), d1 = 0, d2 = radius)
		
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
				parallel::clusterExport(cl, c('x', 'nb', 'pairwiseD'))
			}
		}	
		

		cellVals <- pbapply::pblapply(1:nrow(x[[1]]), function(i) {
		# for (i in 1:nrow(x[[1]])) {
	
			focalCell <- i
			nbCells <- nb[[i]]
		
			# convert to cell indices
			focalCell <- x[['cellCommInd']][focalCell]
			nbCells <- x[['cellCommInd']][nbCells]
		
			
			# visual verification of focal cell + neighborhood
			# plot(st_geometry(x[[1]])[nbCells,], border = 'white')
			# plot(st_geometry(x[[1]]), add = TRUE, border = 'gray95')
			# plot(st_geometry(x[[1]])[nbCells,], add = TRUE)
			# plot(st_geometry(x[[1]])[focalCell, ], col = 'red', add = TRUE)
			# plot(st_buffer(st_centroid(st_geometry(x[[1]][focalCell,])), dist = radius), add = TRUE, border = 'blue')
			
			return(mean(pairwiseD[focalCell, nbCells], na.rm = TRUE))
	
		}, cl = cl)
		
	} else if (inherits(x[[1]], 'SpatRaster')) {
		
		# write two versions, one low memory, one high memory
		
		if (!slow) {
			# convert grid cells to points, and get neighborhoods
			datCells <- which(terra::values(x[[1]]['spRichness']) > 0)
			gridCentroids <- terra::xyFromCell(x[[1]], cell = datCells)
			gridCentroids <- sf::st_as_sf(as.data.frame(gridCentroids), coords = 1:2, crs = attributes(x)$crs)
			gridCentroids$cellInd <- datCells
			nb <- spdep::dnearneigh(gridCentroids, d1 = 0, d2 = radius)
			
			message('\tgridcell neighborhoods: median ', stats::median(lengths(nb)), ', range ', min(lengths(nb)), ' - ', max(lengths(nb)), ' cells')
			
			# prep for parallel computing
			if (nThreads == 1) {
				cl <- NULL
			} else if (nThreads > 1) {
				if (.Platform$OS.type != 'windows') {
					cl <- nThreads
				} else {
					cl <- parallel::makeCluster(nThreads)
					parallel::clusterExport(cl, c('x', 'datCells', 'nb', 'pairwiseD'))
				}
			}	
			
			cellVals <- pbapply::pblapply(1:nrow(gridCentroids), function(i) {
			# for (i in 1:nrow(x[[1]])) {
		
				focalCell <- datCells[i]
				nbCells <- nb[[i]]
				nbCells <- sf::st_drop_geometry(gridCentroids[nbCells, 'cellInd'])[,1]
			
				# convert to cell indices
				focalCell <- x[['cellCommInd']][focalCell]
				nbCells <- x[['cellCommInd']][nbCells]
			
				
				# visual verification of focal cell + neighborhood
				# plot(st_geometry(x[[1]])[nbCells,], border = 'white')
				# plot(st_geometry(x[[1]]), add = TRUE, border = 'gray95')
				# plot(st_geometry(x[[1]])[nbCells,], add = TRUE)
				# plot(st_geometry(x[[1]])[focalCell, ], col = 'red', add = TRUE)
				# plot(st_buffer(st_centroid(st_geometry(x[[1]][focalCell,])), dist = radius), add = TRUE, border = 'blue')
				
				return(mean(pairwiseD[focalCell, nbCells], na.rm = TRUE))
		
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
					
					return(mean(pairwiseD[focalCell, nbCells], na.rm=TRUE))
							
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
	cellVals[is.nan(cellVals)] <- NA 

	metricName <- paste0('phylogenetic_', component)
	
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



