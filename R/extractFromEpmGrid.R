##' @title Extract from epmGrid
##'
##' @description Return species from intersection between spatial points or
##'   polygons and a epmGrid object.
##'
##' @param x object of class \code{epmGrid}

##' @param spatial coordinates as either a spatial points object (sp or sf), a
##'   matrix/dataframe with two columns, a numeric vector of c(long, lat), or as
##'   a spatial polygon object (sp or sf).
##'
##' @param returnCells boolean, if \code{TRUE}, cell indices are returned rather
##'   than taxa
##' @param collapse boolean; if \code{TRUE}, then a vector of unique species is
##'   returned, pooled from all cells, if \code{FALSE}, then list is returned
##'   with species from every cell as intersected by \code{spatial}.
##'
##' @details If \code{spatial} is a spatial object, it will be transformed to
##'   the same projection as \code{x} if needed. If \code{spatial} is not a
##'   spatial object, it is assumed to be in the same projection as \code{x}.
##'
##' @return A vector of species if \code{collapse = TRUE}, or a list of species
##'   by cell if \code{collapse = FALSE}. If \code{returnCells = TRUE}, a vector
##'   of cell indices that correspond to the rows in the epmGrid sf object.
##'
##' @author Pascal Title
##'
##' @examples
##' library(sf)
##' # get the projection of the epmGrid object
##' proj <- summary(tamiasEPM)$crs
##' # define some points
##' pts <- rbind.data.frame(
##' 		c(-120.5, 38.82),
##' 		c(-84.02, 42.75),
##' 		c(-117.95, 55.53))
##' colnames(pts) <- c('x', 'y')
##' ptsSF <- st_as_sf(pts, coords = 1:2, crs = "epsg:4326")
##' pts <- st_coordinates(st_transform(ptsSF, crs = proj))
##'
##' # extract with table of coordinates
##' extractFromEpmGrid(tamiasEPM, pts)
##'
##' # extract with spatial points object
##' extractFromEpmGrid(tamiasEPM, ptsSF)
##'
##' # extract with spatial polygon
##' hull <- st_convex_hull(st_union(ptsSF))
##' extractFromEpmGrid(tamiasEPM, hull)
##'
##'
##' # returns each cell's contents
##' extractFromEpmGrid(tamiasEPM, hull, collapse=FALSE)
##'
##' # collapses results to unique set of species
##' extractFromEpmGrid(tamiasEPM, hull, collapse=TRUE)
##'
##' @export


extractFromEpmGrid <- function(x, spatial, returnCells = FALSE, collapse=TRUE) {
	
	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	if (inherits(spatial, c('SpatialPolygons', 'SpatialPolygonsDataFrame', 'SpatialPoints', 'SpatialPointsDataFrame', 'sf', 'sfc'))) {
		
		if (inherits(spatial, c('SpatialPolygons', 'SpatialPolygonsDataFrame', 'SpatialPoints', 'SpatialPointsDataFrame'))) {
			spatial <- sf::st_as_sf(spatial)
		}
		
		if (!is.na(sf::st_crs(spatial))) {
			if (!identical(sf::st_crs(spatial), sf::st_crs(x[[1]]))) {
				spatial <- sf::st_transform(spatial, crs = sf::st_crs(x[[1]]))
			}
		} else {
			stop('spatial must have projection information.')
		}
			
	} else if (is.numeric(spatial) & is.vector(spatial)) {
		spatial <- as.data.frame(matrix(spatial, nrow = 1, ncol = 2))
		spatial <- as.data.frame(spatial)
		colnames(spatial) <- c('x', 'y')
		spatial <- sf::st_as_sf(spatial, coords = 1:2, crs = sf::st_crs(x[[1]]))	
	} else if (is.matrix(spatial) | is.data.frame(spatial) & !inherits(spatial, 'sf')) {
		spatial <- as.data.frame(spatial)
		colnames(spatial) <- c('x', 'y')
		spatial <- sf::st_as_sf(spatial, coords = 1:2, crs = sf::st_crs(x[[1]]))	
	} else {
		stop('spatial argument not understood.')
	}
			
			
	spatial <- sf::st_geometry(spatial)
	
	if (inherits(x[[1]], 'sf')) {
		cells <- unlist(sf::st_intersects(spatial, x[[1]]))
	} else if (inherits(x[[1]], 'SpatRaster')) {
		if (unique(as.character(sf::st_geometry_type(spatial))) == 'POINT') {
			cells <- terra::cellFromXY(x[[1]], sf::st_coordinates(spatial))
		} else if (unique(as.character(sf::st_geometry_type(spatial))) %in% c('MULTIPOLYGON', 'POLYGON', 'GEOMETRY')) {
			cells <- terra::cells(x[[1]], terra::vect(spatial))[, 2]
		} else {
			stop('Spatial type not recognized.')
		}
	} else {
		stop('EPM grid format not recognized.')
	}
	
	if (returnCells) {
		if (collapse) {
			return(sort(unique(cells)))
		} else {
			return(sort(cells))
		}		
	} else {
			
		# extract species from epmGrid
		res <- sapply(x[['cellCommInd']][cells], function(y) x[['speciesList']][[y]], simplify = FALSE)
		names(res) <- paste0('cell', cells)
		
		if (collapse) {
			return(sort(unique(unlist(res))))
		} else {
			return(res)
		}
	}
}	
		
		
		
