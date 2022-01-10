##' @title Convert epmGrid to community matrix
##'
##' @description Given specific sites, convert epmGrid to 
##' 	phylocomm matrix, with sites as rows, and species as columns
##'
##' @param x object of class \code{epmGrid}
##' @param sites locations of sites, see details. 
##'
##' @details If sites are site coordinates, 
##' 	then dataframe or matrix with two columns; 
##' 	if sites are cell indices, then numeric vector;
##' 	if \code{sites = 'all'}, then all cells will be returned as sites. 
##'
##' @return community matrix, with sites as rows and species as columns
##'
##' @author Pascal Title
##'
##' @examples
##' tamiasEPM
##'
##' # from cell indices
##' cells <- c(2703, 90, 3112, 179)
##' epmToPhyloComm(tamiasEPM, cells)
##' 
##' # from coordinates
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
##' epmToPhyloComm(tamiasEPM, pts)
##' 
##' @export

epmToPhyloComm <- function(x, sites) {
	
	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	# if sites are in table format, then must be coordinates
	if (inherits(sites, c('matrix', 'data.frame'))) {
		sites <- extractFromEpmGrid(x, spatial = sites, returnCells = TRUE, collapse = TRUE)

	# if sites is a string (presumably 'all'), then we will use all sites
	} else if (is.character(sites)) {
		if (inherits(x[[1]], 'SpatRaster')) {
			sites <- 1:terra::ncell(x[[1]])	
		} else if (inherits(x[[1]], 'sf')) {	
			sites <- 1:nrow(x[[1]])
		}
	}
	
	# extract relevant cells from empGrid
	dat <- sapply(x[['cellCommInd']][sites], function(y) x[['speciesList']][[y]], simplify = FALSE)
	uniqueSp <- sort(unique(unlist(dat)))
	
	if (all(is.null(uniqueSp))) {
		# build phylocomm matrix
		resMat <- matrix(nrow = length(dat), ncol = length(uniqueSp))
		rownames(resMat) <- paste0('site', 1:length(dat))
		message('\tNo species present at sites.')

	} else {
	
		uniqueSp <- uniqueSp[stats::complete.cases(uniqueSp)]
		
		# build phylocomm matrix
		resMat <- matrix(nrow = length(dat), ncol = length(uniqueSp))
		rownames(resMat) <- paste0('site', 1:length(dat))
		colnames(resMat) <- uniqueSp
		
		for (i in 1:nrow(resMat)) {
			resMat[i,] <- as.numeric(uniqueSp %in% dat[[i]])
		}
	}
	
	return(resMat)
}