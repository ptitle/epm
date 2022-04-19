##' @export

print.epmGrid <- function(x, ...) {
	
	if (!inherits(x, 'epmGrid')) {
		stop('Object must be of class epmGrid.')
	}
	
	# if data present in object, then report info
	if (inherits(x[['data']], c('numeric', 'matrix', 'data.frame'))) {
		if (is.vector(x[['data']])) {
			data <- length(intersect(x[['geogSpecies']], names(x[['data']])))
		} else {
			data <- length(intersect(x[['geogSpecies']], rownames(x[['data']])))
		}
	} else {
		data <- NA
	}
	
	# if phylogeny present in object, then report info
	if (inherits(x[['phylo']], 'phylo')) {
		x[['phylo']] <- list(x[['phylo']])
		class(x[['phylo']]) <- 'multiPhylo'
	}
	
	if (inherits(x[['phylo']], 'multiPhylo')) {
		phylo <- length(intersect(x[['geogSpecies']], x[['phylo']][[1]]$tip.label))
		nTrees <- length(x[['phylo']])
	} else {
		phylo <- NA
		nTrees <- 0
	}
	
	metric <- attributes(x)$metric

	if (inherits(x[[1]], 'sf')) {
		ncells <- nrow(x[['grid']])
		gridExtent <- sf::st_bbox(x[['grid']])
	} else if (inherits(x[[1]], 'SpatRaster')) {
		ncells <- terra::ncell(x[['grid']])
		gridExtent <- as.vector(terra::ext(x[['grid']]))
	} else {
		stop('Grid format not recognized.')
	}
	
	resolution <- attributes(x)$resolution
	proj <- attributes(x)$crs
	isProj <- attributes(x)$projected
	lengthUniqueSp <- length(x[['geogSpecies']])
	minSp <- min(lengths(x[['speciesList']]))
	maxSp <- max(lengths(x[['speciesList']]))
	
	cat('\n\tSummary of epm object:\n\n')
	cat('\tmetric:', metric, '\n')
	cat('\tgrid type: ', attributes(x)$gridType, '\n')
	cat('\tnumber of grid cells:', ncells, '\n')
	cat('\tgrid resolution:', resolution, 'by', resolution, '\n')
	cat('\tprojected:', isProj, '\n')
	cat('\tcrs:', proj, '\n\n')
	cat(paste0('\tnumber of unique species: ', lengthUniqueSp, ' (richness range: ', minSp, ' - ', maxSp, ')'), '\n')
	cat('\tdata present:', ifelse(is.na(data), 'No', 'Yes'), '\n')
	if (!is.na(data)) {
		cat('\tnumber of species shared between data and grid:', data, '\n')
	}

	cat('\tphylogeny present:', ifelse(is.na(phylo), 'No', paste0('Yes (', nTrees, ifelse(nTrees > 1, ' trees)', ' tree)'))), '\n')
	if (!is.na(phylo)) {
		cat('\tnumber of species shared between phylogeny and grid:', phylo)
	}
	
	cat('\n')
	
}