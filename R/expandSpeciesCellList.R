##' Function to expand condensed species list to full set of cells

##' @title Expand species list
##'
##' @description The epmGrid object contains an accounting of species per
##' cell in a condensed format. This function returns a complete list of species
##' per cell. 
##'
##' @param x object of class \code{epmGrid}
##'
##' @return list of species for each cell.
##'
##' @author Pascal Title
##' 
##' @examples
##' tamiasEPM
##'	head(expandSpeciesCellList(tamiasEPM))
##' 
##' @export



expandSpeciesCellList <- function(x) {
	
	if (!'epmGrid' %in% class(x)) {
		stop('x must be of class epmGrid.')
	}	
	
	return(sapply(x[['cellCommInd']], function(y) x[['speciesList']][[y]]))
	
}