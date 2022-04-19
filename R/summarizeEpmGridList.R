##' @title Summarize lists of epmGrid objects
##'
##' @description If a diversity metric was calculated for an epmGrid object
##' that contained a phylogenetic distribution, then a list of resulting
##' epmGrid objects was returned. This function will take that list, and 
##' apply a summary statistic, returning a single epmGrid object. 
##' If the input list is from \code{\link{betadiv_phylogenetic}}, then that
##' list of sf or SpatRaster objects can also be summarized with this function.
##'
##' @param x a list of objects of class \code{epmGrid} or \code{sf} or \code{SpatRaster}. 
##'
##'@param fun a function to apply to grid cells across the list x. 
##' 
##'
##' @details It is assumed that across the objects in list x,
##' the only difference is the values for the grid cells. 
##'
##'
##' @return a single object of class \code{epmGrid} or \code{sf} or \code{SpatRaster}.  
##'
##' @author Pascal Title
##' 
##' @examples
##'
##' tamiasEPM <- addPhylo(tamiasEPM, tamiasTreeSet)
##' tamiasEPM <- addTraits(tamiasEPM, tamiasTraits)
##' 
##' x <- gridMetrics(tamiasEPM, metric='meanPatristicNN')
##' z <- summarizeEpmGridList(x, fun = var)
##' 
##' # using a custom function
##' f <- function(y) sum(y) / length(y)
##' 
##' z <- summarizeEpmGridList(x, fun = f)
##' 
##' # works with square grid epmGrids too
##' tamiasEPM2 <- createEPMgrid(tamiasPolyList, resolution = 50000, 
##' 	cellType = 'square', method = 'centroid')
##' tamiasEPM2 <- addPhylo(tamiasEPM2, tamiasTreeSet)
##' tamiasEPM2 <- addTraits(tamiasEPM2, tamiasTraits)
##' 
##' x <- gridMetrics(tamiasEPM2, metric='meanPatristicNN')
##' 
##' z <- summarizeEpmGridList(x, fun = median)
##' 
##' \donttest{
##' # With phylogenetic distribution
##'
##' tamiasEPM <- addPhylo(tamiasEPM, tamiasTreeSet, replace = TRUE)
##' beta_phylo_turnover <- betadiv_phylogenetic(tamiasEPM, radius = 70000,
##' 	component = 'turnover')
##' 
##' z <- summarizeEpmGridList(beta_phylo_turnover)
##' 
##' tamiasEPM2 <- createEPMgrid(tamiasPolyList, resolution = 50000,
##' 	cellType = 'square', method = 'centroid')
##' tamiasEPM2 <- addPhylo(tamiasEPM2, tamiasTreeSet)
##' 
##' beta_phylo_turnover <- betadiv_phylogenetic(tamiasEPM2, radius = 70000,
##' 		component = 'turnover')
##' 
##' z <- summarizeEpmGridList(beta_phylo_turnover, fun = median)
##'	
##'
##' }
##' @export


summarizeEpmGridList <- function(x, fun = mean) {
	
	# check: we expect that all elements in list are (1) epmGrid with the same metric, (2) a list of sf objects with a common data column, or (3) a set of SpatRasters
	if (!all(sapply(x, inherits, c('epmGrid', 'sf', 'SpatRaster')))) {
		stop('Not all list elements are of class epmGrid.')
	}



	# for list of epmGrid objects
	
	if (all(sapply(x, inherits, 'epmGrid'))) {
	
		if (length(unique(sapply(x, function(y) attributes(y)$metric))) != 1) {
			stop('Not all epmGrid objects are for the same metric.')
		}
		
		metricCol <- attributes(x[[1]])$metric
		
		if (inherits(x[[1]][[1]], 'sf')) {
			
			xTables <- lapply(x, function(y) sf::st_drop_geometry(y[[1]]))
			vals <- matrix(nrow = nrow(x[[1]][[1]]), ncol = length(x))
			for (i in 1:nrow(x[[1]][[1]])) {
				vals[i,] <- vapply(xTables, function(y) y[i, metricCol], FUN.VALUE = numeric(1))
			}
			
			newvals <- apply(vals, 1, function(y) fun(y))
			
			ret <- x[[1]]
			ret[[1]][, metricCol] <- newvals
			
		} else if (inherits(x[[1]][[1]], 'SpatRaster')) {
					
			# combine all into a single multilayer SpatRaster
			xx <- do.call(c, lapply(x, function(y) y[[1]][metricCol]))
			ret <- x[[1]]
			terra::values(ret[[1]][[metricCol]]) <- as.numeric(terra::values(terra::app(xx, fun = fun)))
			
		} else {
			stop('Grid system format not recognized.')
		}
		
	
	# for list of sf objects
	} else if (all(sapply(x, inherits, 'sf'))) {
	
		if (length(unique(lapply(x, colnames))) != 1) {
			stop('sf objects do not have the same column names.')
		}
		
		xTables <- lapply(x, function(y) sf::st_drop_geometry(y))
		
		if (!all(sapply(xTables, ncol) == 1)) {
			stop('sf objects should just have 1 non-geometry column to summarize.')
		}
		
		xTables <- lapply(xTables, function(y) y[, 1])
		xTables <- do.call(rbind, xTables)
		
		newvals <- apply(xTables, 2, function(y) fun(y))
		
		ret <- x[[1]]
		ret[, colnames(sf::st_drop_geometry(x[[1]]))] <- newvals
	
	
	 # for list of SpatRaster objects
	 } else if (all(sapply(x, inherits, 'SpatRaster'))) {
	 	
	 	# there should be just one layer per SpatRaster
	 	if (length(unique(sapply(x, names))) != 1) {
	 		stop('Not all inputs have the same layer.')
	 	}
	 	
 		# combine all into a single multilayer SpatRaster
		xx <- do.call(c, lapply(x, function(y) y))
		ret <- terra::rast(x[[1]])
		terra::values(ret) <- as.numeric(terra::values(terra::app(xx, fun = fun)))

	 } else {
	 	stop('Input is either a mix of classes or is not supported.')
	 }
			
	return(ret)
}