##' @title Retrieve coordinates from epmGrid
##'
##' @description Return the centroid coordinates for a specified set of grid cells.
##'
##' @param x object of class \code{epmGrid}
##' @param sites locations of sites, see details. 
##'
##' @details 
##' 	Sites can be cell indices as a numeric vector, or you can specify 
##' 	\code{sites = 'all'} to get all grid cells. If the epmGrid object
##' 	is hexagon-based, then all grid cells that are occupied are returned.
##' 	If the epmGrid is square-based, then all grid cells, occupied or empty,
##' 	are returned.
##'
##' @return matrix with x and y coordinates.
##'
##' @author Pascal Title
##'
##' @examples
##' tamiasEPM
##'
##' # from cell indices
##' cells <- c(2703, 90, 3112, 179)
##' coordsFromEpmGrid(tamiasEPM, cells)
##'
##' # for all grid cells
##' dim(coordsFromEpmGrid(tamiasEPM, 'all'))
##' 
##' 
##' @export

coordsFromEpmGrid <- function(x, sites) {
	
	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	if (inherits(x[[1]], 'sf')) {
	
		# if sites is a string (presumably 'all'), then we will use all sites
		if (is.character(sites)) {		
			sites <- 1:nrow(x[[1]])
		}
	
		coords <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(x[[1]][sites,])))[, 1:2]
		
	} else {
		
		# if sites is a string (presumably 'all'), then we will use all sites
		if (is.character(sites)) {
			sites <- 1:terra::ncell(x[[1]])
		}
		
		coords <- terra::xyFromCell(x[[1]], sites)
	}

	colnames(coords) <- c('x', 'y')
	rownames(coords) <- NULL
	return(coords)
}