##' @title Read a epmGrid object
##'
##' @description Load a saved epmGrid object. 
##'
##' @param filename filename, with extension \code{rds}
##'
##' @details This function will read in epmGrid objects that
##' were saved with \code{\link{write.epmGrid}}. 
##'
##' @return object of class epmGrid
##'
##' @author Pascal Title
##'
##' @examples
##' #save
##' write.epmGrid(tamiasEPM, paste0(tempdir(), '/tamiasEPM'))
##' 
##' # read back in
##' tamiasEPM <- read.epmGrid(paste0(tempdir(), '/tamiasEPM.rds'))
##'
##' # delete the file
##' unlink(paste0(tempdir(), '/tamiasEPM.rds'))
##'
##' @export

read.epmGrid <- function(filename) {
	
	if (!grepl('\\.rds$', filename)) {
		filename <- paste0(filename, '.rds')
	}
	
	x <- readRDS(filename)
	
	if (inherits(x[[1]], 'PackedSpatRaster')) {
		x[[1]] <- terra::rast(x[[1]])
	}

	return(x)	
}



