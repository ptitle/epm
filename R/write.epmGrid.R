##' @title Save epmGrid object
##'
##' @description Write a epmGrid object to disk. 
##'
##' @param x object of class \code{epmGrid}
##' @param filename filename with no extension
##'
##' @details This function writes a .rds file with xz compression. 
##' This file can be read back in with \code{\link{read.epmGrid}}.
##'
##' @return Nothing is returned, but object is written to disk.
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
##'
##' @export

write.epmGrid <- function(x, filename) {
	
	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	if (!grepl('\\.rds$', filename)) {
		filename <- paste0(filename, '.rds')
	}
	
	if (inherits(x[[1]], 'SpatRaster')) {
	    x[[1]] <- terra::wrap(x[[1]])
	}
	
	saveRDS(x, file = filename, compress = 'xz')	
}



