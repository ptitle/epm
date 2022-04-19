##' @title Interactively identify species in epmGrid
##'
##' @description Plots a epmGrid object and allows you to click on the plot
##' 	to return the species found in the cell you clicked on. 
##'
##' @param x object of class \code{epmGrid} or \code{sf}
##' @param returnCell boolean; if FALSE, then species names are returned,
##' 	if TRUE, then cell indices are returned. 
##' @param ... additional arguments passed to sf::plot
##'
##' @details This is a wrapper function for the \code{identify} function
##' in base graphics. This is primarily intended as a useful function for 
##' data exploration and spot-checking. 
##' 
##' @return A vector of species names or cell indices. 
##'
##' @author Pascal Title
##' 
##' @rdname identify
##' @aliases identify.epmGrid
##' @export
##' @export identify.epmGrid

identify.epmGrid <- function(x, returnCell = FALSE, ...) {
	
	colramp <- function(n) viridisLite::turbo(n = n, begin = 0.1, end = 0.9)
	
	if (inherits(x, 'epmGrid')) {
		plotMetric <- attributes(x)$metric
		
		if (inherits(x[[1]], 'sf')) {
			
			plot(x[[1]][plotMetric], lwd = 0.5, reset = FALSE, pal = colramp, ...)
			cell <- identify(sf::st_coordinates(sf::st_centroid(sf::st_geometry(x[[1]]))), n = 1, labels = '')
			
			plot(sf::st_centroid(sf::st_geometry(x[[1]])[cell]), col = 'orange', add = TRUE, pch = 3)
			
		} else if (inherits(x[[1]], 'SpatRaster')) {
			
			plot.epmGrid(x, use_tmap = FALSE)
			cell <- terra::click(x[[1]][plotMetric], n = 1, cell = TRUE)$cell
			
			graphics::points(terra::xyFromCell(x[[1]], cell), col = 'orange', pch = 3)
			
		} else {
			stop('Grid format not recognized.')
		}
		
		if (returnCell) {
			return(cell)
		} else {
			return(x[['speciesList']][[x[['cellCommInd']][cell]]])
		}

	} else if (inherits(x, 'sf')) {
		if (ncol(x) > 2) {
			plotMetric <- utils::tail(setdiff(colnames(x), 'sf_column'), 1)
		} else {
			plotMetric <- setdiff(colnames(x), 'sf_column')
		}
		plot(x[plotMetric], lwd = 0.5, reset = FALSE)
		cell <- identify(sf::st_coordinates(sf::st_centroid(sf::st_geometry(x))), n = 1, labels = '')
		plot(sf::st_centroid(sf::st_geometry(x)[cell,]), col = 'orange', add = TRUE, pch = 3)
		return(cell)
	}

}








	
