##'@title Custom beta diversity metrics
##'
##'@description Define your own function for summarizing information across a
##'  moving window of grid cells.
##'
##'@param x object of class \code{epmGrid}
##'@param fun a function to apply to grid cell neighborhoods (see details)
##'@param radius Radius of the moving window in map units.
##'@param minTaxCount the minimum number of taxa needed to apply the function.
##'  For instance, should the function be applied to gridcells with just 1
##'  taxon?
##'@param focalCoord vector of x and y coordinate, see details
##'@param metricName the name you would like to attach to the output
##'
##'@details This function will identify the neighbors of every cell and will
##' apply the specified function to those sets of cell neighborhoods. 
##'
##'  The custom function should have just one input: a list of taxon names, where
##'  the list will represent a set of grid cells (focal cell + neighboring cells).
##' 
##'  However, if a set of focal coordinates is provided, then rather than apply the
##'  function to each neighborhood of cells, the function should have two inputs: the
##'  focal cell and another cell, and that function will be applied to every pair 
##'  defined by the focal cell and another cell. See examples.
##' 
##'  Within the function call, the trait data already attached to the epmGrid object 
##'  must be referred to as dat, and the phylogenetic tree already attached to the 
##'  epmGrid must be referred to as phylo.\cr 
##'
##'  If the input epmGrid object contains a set of trees, then this function will
##' 	 be applied, using each tree in turn, and will return a list of results. This
##'  list can then be passed to \code{\link{summarizeEpmGridList}} to be summarized.
##'  
##'  See examples below.
##'
##'@return object of class \code{epmGrid}, or list of \code{epmGrid} objects
##'
##'@author Pascal Title
##'
##' @examples
##' \donttest{
##' tamiasEPM <- addPhylo(tamiasEPM, tamiasTree)
##' tamiasEPM <- addTraits(tamiasEPM, tamiasTraits)
##' 
##' # An example using a multivariate dataset
##' ## For each focal cell + neighbors, calculate the morphological range
##' ## per grid cell and return the standard deviation.
##' f <- function(cellList) {
##' 	vec <- numeric(length(cellList))
##' 	for (i in 1:length(cellList)) {
##' 		vec[[i]] <- max(dist(dat[cellList[[i]], ]))
##' 	}
##' 	return(sd(vec, na.rm = TRUE))	
##' }
##' 
##' xx <- customBetaDiv(tamiasEPM, fun = f, radius = 70000, minTaxCount = 2, metricName = 'maxdist')
##' 
##' 
##' # An example using only the phylogeny.
##' ## Calculate standard deviation of phylogenetic diversity across cell neighborhood.
##' f <- function(cellList) {
##' 	vec <- numeric(length(cellList))
##' 	for (i in 1:length(cellList)) {
##' 		vec[[i]] <- faithPD(phylo, cellList[[i]])
##' 	}
##' 	return(sd(vec, na.rm = TRUE))	
##' }
##' 
##' xx <- customBetaDiv(tamiasEPM, fun = f, radius = 70000, minTaxCount = 1, metricName = 'faithPD')
##' 
##' 
##' 
##' # an example that involves both morphological and phylogenetic data
##' ## nonsensical, but for illustrative purposes:
##' ## ratio of Faith's phylogenetic diversity to morphological range
##' ## the standard deviation of this measure across grid cells 
##' ## in a neighborhood.
##' f <- function(cellList) {
##' 	vec <- numeric(length(cellList))
##' 	for (i in 1:length(cellList)) {
##' 		vec[[i]] <- faithPD(phylo, cellList[[i]]) / 
##' 		max(dist(dat[cellList[[i]], ]))
##' 	}
##' 	return(sd(vec, na.rm = TRUE))	
##' }
##' 
##' xx <- customBetaDiv(tamiasEPM, fun = f, radius = 70000, minTaxCount = 2, 
##'   metricName = 'ratio_PD_maxdist')
##' 
##'
##' 
##' 
##' # from a focal coordinate to all other sites
##' ## Here, the function has 2 inputs.
##' ## Example: calculate the per grid cell mean and take the distance.
##' f <- function(focalCell, otherCell) {
##' 	x1 <- colMeans(dat[focalCell, ])
##' 	x2 <- colMeans(dat[otherCell, ])
##' 	return(as.matrix(dist(rbind(x1, x2)))[1,2])
##' }
##' 
##' xx <- customBetaDiv(tamiasEPM, fun = f, radius = 70000, minTaxCount = 1,
##'  focalCoord = c(-1413764, 573610.8), metricName = 'meandist')
##' 
##' 
##' 
##' # Example involving a set of trees
##' tamiasEPM <- addPhylo(tamiasEPM, tamiasTreeSet, replace = TRUE)
##' 
##' ## Calculate standard deviation of phylogenetic diversity across cell
##' ## neighborhood.
##' f <- function(cellList) {
##' 	vec <- numeric(length(cellList))
##' 	for (i in 1:length(cellList)) {
##' 		vec[[i]] <- faithPD(phylo, cellList[[i]])
##' 	}
##' 	return(sd(vec, na.rm = TRUE))	
##' }
##' 
##' # This time, a list of sf objects will be returned, one for each input tree.
##' xx <- customBetaDiv(tamiasEPM, fun = f, radius = 70000, minTaxCount = 1,
##'  metricName = 'faithPD')
##' 
##' 
##' 
##' # also works with square grid cells
##' tamiasEPM2 <- createEPMgrid(tamiasPolyList, resolution = 50000,
##' 	cellType = 'square', method = 'centroid')
##' tamiasEPM2 <- addPhylo(tamiasEPM2, tamiasTree)
##' tamiasEPM2 <- addTraits(tamiasEPM2, tamiasTraits)
##' 
##' 
##' f <- function(cellList) {
##' 	vec <- numeric(length(cellList))
##' 	for (i in 1:length(cellList)) {
##' 		vec[[i]] <- faithPD(phylo, cellList[[i]]) / 
##' 		max(dist(dat[cellList[[i]], ]))
##' 	}
##' 	return(sd(vec, na.rm = TRUE))	
##' }
##' 
##' xx <- customBetaDiv(tamiasEPM2, fun = f, radius = 70000, minTaxCount = 2,
##'  metricName = 'ratio_PD_maxdist')
##' }
##' 
##'
##'@export



# where phenotypic or other species-specific data = dat
# where phylogeny is phylo
# where the only input variable can be the list of grid cell taxa that include a focal cell and its neighbors.
customBetaDiv <- function(x, fun, radius, minTaxCount = 1, focalCoord = NULL, metricName = 'custom_metric') {

	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	if (inherits(x[['phylo']], 'phylo')) {
		x[['phylo']] <- list(x[['phylo']])
		class(x[['phylo']]) <- 'multiPhylo'
	}
	
	# what is being called in the function fun?
	useDat <- FALSE
	usePhylo <- FALSE
	if (any(grepl('dat', as.character(body(fun))))) {
		useDat <- TRUE
	}
	if (any(grepl('phylo', as.character(body(fun))))) {
		usePhylo <- TRUE
	}
	if (length(names(formals(fun))) != 1 & is.null(focalCoord)) {
		stop('function should only have one input argument.')
	}

	if (length(names(formals(fun))) != 2 & !is.null(focalCoord)) {
		stop('function should two input arguments for pairwise calculations.')
	}
		
	if (useDat) {
		if (is.null(x[['data']])) {
			stop('epmGrid object does not contain trait data!')
		}
	}
	
	if (usePhylo) {
		if (is.null(x[['phylo']])) {
			stop('epmGrid object does not contain a phylo object!')
		}
	}
				
	# prune epmGrid
	if (useDat & !usePhylo) {
		## if only traits are involved	
		# prune epmGrid object down to species shared with data
		if (is.vector(x[['data']])) {
			x[['speciesList']] <- intersectList(x[['speciesList']], names(x[['data']]))
		} else {
			x[['speciesList']] <- intersectList(x[['speciesList']], rownames(x[['data']]))
		}
	}
	
	if (!useDat & usePhylo) {
	 	# prune epmGrid object down to species shared with phylogeny
		x[['speciesList']] <- intersectList(x[['speciesList']], x[['phylo']][[1]]$tip.label)
	}
	
	if (useDat & usePhylo) {
		if (is.vector(x[['data']])) {
			sharedTaxa <- intersect(names(x[['data']]), x[['phylo']][[1]]$tip.label)
			x[['data']] <- x[['data']][sharedTaxa]
		} else {
			sharedTaxa <- intersect(rownames(x[['data']]), x[['phylo']][[1]]$tip.label)
			x[['data']] <- x[['data']][sharedTaxa,]
		}
		x[['phylo']] <- lapply(x[['phylo']], ape::keep.tip, sharedTaxa)
		x[['speciesList']] <- intersectList(x[['speciesList']], sharedTaxa)
	}

	# There appear to be some potential numerical precision issues with defining neighbors.
	# Therefore, if radius is a multiple of resolution, add 1/100 of the resolution. 
	if (radius %% attributes(x)$resolution == 0) {
		radius <- radius + attributes(x)$resolution * 0.01
	}

	nRuns <- 1
	if (usePhylo & length(x[['phylo']]) > 1) {
		nRuns <- length(x[['phylo']])
		op <- getOption("pboptions")
		opb <- pbapply::pboptions(type = 'none')
	}
	
	retList <- vector('list', length = nRuns)
	
	if (is.null(focalCoord)) {
		# generate neighborhoods
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
	
		} else if (inherits(x[[1]], 'SpatRaster')) {
				
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
		}
		
		message('\tgridcell neighborhoods: median ', stats::median(lengths(nb)), ', range ', min(lengths(nb)), ' - ', max(lengths(nb)), ' cells')
	
	}

	if (nRuns > 1) pb <- utils::txtProgressBar(max = length(x[['phylo']]), char = '+', style = 3)
	for (i in 1:nRuns) {
		
		if (nRuns > 1) utils::setTxtProgressBar(pb, i)
		
		# create custom environment for function and add in trait and phylo data
		.epm.env <- new.env()
		assign('dat', x[['data']], envir = .epm.env)
		assign('phylo', x[['phylo']][[i]], envir = .epm.env)
		# ls(envir = .epm.env)
		environment(fun) <- .epm.env

		if (!is.null(focalCoord)) {
					
			# which community does the focal coordinate represent?
			focalSp <- extractFromEpmGrid(x, focalCoord)
			focalInd <- which(sapply(x[['speciesList']], function(y) identical(focalSp, y)) == TRUE)
			if (length(focalInd) != 1) stop('focalInd != 1.')
						
			pairwiseMat <- matrix(nrow = length(x[['speciesList']]), ncol = length(x[['speciesList']]))
			for (k in 1:length(x[['speciesList']])) {
				pairwiseMat[focalInd, k] <- fun(focalSp, x[['speciesList']][[k]])
			}

		
			if (inherits(x[[1]], 'sf')) {
				cellVals <- pbapply::pbsapply(1:nrow(x[[1]]), function(k) pairwiseMat[focalInd, x[['cellCommInd']][k]])
			} else if (inherits(x[[1]], 'SpatRaster')) {
				cellVals <- pbapply::pbsapply(1:terra::ncell(x[[1]]), function(k) pairwiseMat[focalInd, x[['cellCommInd']][k]])
			} else {
				stop('epm format not recognized.')
			}

		} else {
		
			if (inherits(x[[1]], 'sf')) {
							
				cellVals <- pbapply::pbsapply(1:nrow(x[[1]]), function(k) {
				# for (k in 1:nrow(x[[1]])) {
		
					cells <- x[['cellCommInd']][c(k, nb[[k]])]
					allsites <- x[['speciesList']][cells]
					allsites <- allsites[lengths(allsites) >= minTaxCount]
					if (length(allsites) > 1) {
						fun(allsites)
					} else {
						NA
					}	
				})							
			} else if (inherits(x[[1]], 'SpatRaster')) {

				cellVals <- rep(NA, terra::ncell(x[[1]]))
				cellVals[datCells] <- pbapply::pbsapply(1:nrow(gridCentroids), function(k) {
				# for (k in 1:nrow(gridCentroids)) {
			
					nbCells <- sf::st_drop_geometry(gridCentroids[c(k, nb[[k]]), 'cellInd'])[,1]
					nbCells <- x[['cellCommInd']][nbCells]
					allsites <- x[['speciesList']][nbCells]
				
					allsites <- allsites[lengths(allsites) >= minTaxCount]
	
					if (length(allsites) > 1) {
						fun(allsites)
					} else {
						NA
					}				
				})
			}
		}

		cellVals <- unlist(cellVals)
		cellVals[is.nan(cellVals)] <- NA 

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
	
	if (nRuns > 1) {
		close(pb)
		opb <- pbapply::pboptions(op)
	}

	
	if (length(retList) == 1) {
		retList <- retList[[1]]
	}
	
	return(retList)
}



