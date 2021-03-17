##' @title Drop species from epmGrid
##'
##' @description Removes particular species from a epmGrid object.
##'
##' @param x object of class \code{epmGrid}
##' @param sp a character vector of species names to be dropped.
##'
##' @details If species in \code{sp} are not in \code{x}, they will be ignored.
##'
##' @return new \code{epmGrid} object.
##'
##' @author Pascal Title
##' 
##' @examples
##' tamiasEPM
##'
##' new <- dropSpecies(tamiasEPM, sp = c('Tamias_alpinus', 'Tamias_bulleri'))
##'
##' setdiff(tamiasEPM[['geogSpecies']], new[['geogSpecies']])
##'	
##' 
##' @export



dropSpecies <- function(x, sp) {
	
	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	if (any(sp %in% x[['geogSpecies']])) {
	
		fullList <- expandSpeciesCellList(x)
		
		for (i in 1:length(fullList)) {		
			if (any(sp %in% fullList[[i]])) {
				tmp <- sort(setdiff(fullList[[i]], sp))
				if (length(tmp) > 0) {
					fullList[[i]] <- tmp
				} else {
					fullList[[i]] <- NA
				}
			}		
		}
		
		return(rebuildEpmGrid(x, fullList))
		
	} else {
		return(x)
	}
}