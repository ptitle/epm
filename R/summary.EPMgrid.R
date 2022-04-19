##' @title epmGrid summary
##'
##' @description Generates a summary of a epmGrid object.
##'
##' @param object object of class \code{epmGrid}
##' @param ... further arguments passed to \code{\link{summary}}
##'
##' @details
##' Summary information includes
##' 
##' @return A list containing the summary information is returned invisibly.
##'
##' @author Pascal Title
##' 
##' @examples
##' summary(tamiasEPM)
##' attr <- summary(tamiasEPM)
##' attr
##' 
##' @rdname summary
##' @aliases summary.epmGrid
##' @export

summary.epmGrid <- function(object, ...) {
	
	if (!inherits(object, 'epmGrid')) {
		stop('Object must be of class epmGrid.')
	}
	
	# if data present in object, then report info
	if (inherits(object[['data']], c('numeric', 'matrix', 'data.frame'))) {
		if (is.vector(object[['data']])) {
			data <- length(intersect(object[['geogSpecies']], names(object[['data']])))
		} else {
			data <- length(intersect(object[['geogSpecies']], rownames(object[['data']])))
		}
	} else {
		data <- NA
	}
	
	# if phylogeny present in object, then report info
	if (inherits(object[['phylo']], 'phylo')) {
		object[['phylo']] <- list(object[['phylo']])
		class(object[['phylo']]) <- 'multiPhylo'
	}

	if (inherits(object[['phylo']], 'multiPhylo')) {
		phylo <- length(intersect(object[['geogSpecies']], object[['phylo']][[1]]$tip.label))
		nTrees <- length(object[['phylo']])
	} else {
		phylo <- NA
		nTrees <- 0
	}
	
	metric <- attributes(object)$metric

	if (inherits(object[[1]], 'sf')) {
		ncells <- nrow(object[['grid']])
		gridExtent <- sf::st_bbox(object[['grid']])
	} else if (inherits(object[[1]], 'SpatRaster')) {
		ncells <- terra::ncell(object[['grid']])
		gridExtent <- as.vector(terra::ext(object[['grid']]))
	} else {
		stop('Grid format not recognized.')
	}

	gridType <- attributes(object)$gridType
	resolution <- attributes(object)$resolution
	proj <- attributes(object)$crs
	isProj <- attributes(object)$projected
	lengthUniqueSp <- length(object[['geogSpecies']])
	minSp <- min(lengths(object[['speciesList']]))
	maxSp <- max(lengths(object[['speciesList']]))
	
	print.epmGrid(object)
		
	obj <- list(
				ncells = ncells, 
				gridType = gridType,
				extent = gridExtent, 
				resolution = resolution, 
				projected = isProj,
				crs = proj, 
				numberUniqueSpecies = lengthUniqueSp, 
				minSp = minSp,
				maxSp = maxSp,
				overlapWithPhylogeny = phylo,
				nTrees = nTrees)
				
	invisible(obj)
}
