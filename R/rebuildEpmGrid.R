##' @title Rebuild epmGrid
##'
##' @description Given a list of species in each cell, all components of the
##' 	epmGrid object are reindexed and regenerated.
##'
##' @param x object of class \code{epmGrid}
##' @param spCellList list in which each element is a character vector of species names.
##'
##' @details This function is used internally by \code{\link{dropSpecies}}. 
##'
##' @return new \code{epmGrid} object.
##'
##' @author Pascal Title
##'	
##' 
##' @export


rebuildEpmGrid <- function(x, spCellList) {
		
	# convert spCellList into condensed version: 
	# 	list of unique cell communities
	#	vector of indices to map these back to grid cells
	if (inherits(x[[1]], 'sf')) {
		cellCommVec <- integer(length = nrow(x[[1]]))
	} else if (inherits(x[[1]], 'SpatRaster')) {
		cellCommVec <- integer(length = terra::ncell(x[[1]]))
	} else {
		stop('Grid format not recognized.')
	}
	uniqueComm <- unique(spCellList)
	fullList2 <- sapply(spCellList, function(y) paste(y, collapse = '|'))
	uniqueComm2 <- sapply(uniqueComm, function(y) paste(y, collapse = '|'))
	for (i in 1:length(uniqueComm2)) {
		cellCommVec[which(fullList2 == uniqueComm2[i])] <- i
	}
	
	uniqueSp <- sort(unique(unlist(uniqueComm)))
	
	spCellCount <- countCells(convertNAtoEmpty(spCellList), uniqueSp)
	names(spCellCount) <- uniqueSp
	
	gridVals <- lengths(spCellList)
	gridVals[which(sapply(spCellList, anyNA))] <- NA

	if (inherits(x[[1]], 'sf')) {
		x[[1]]['spRichness'] <- gridVals
	} else {
		x[[1]][['spRichness']] <- gridVals
	}
	
	x[['speciesList']] <- uniqueComm
	x[['cellCommInd']] <- cellCommVec
	x[['geogSpecies']] <- uniqueSp
	x[['cellCount']] <- spCellCount

	return(x)
}


convertNAtoEmpty <- function(spCellList) {
	ind <- which(sapply(spCellList, anyNA) == TRUE)
	spCellList[ind] <- rep(list('empty'), length(ind))
	return(spCellList)
}
