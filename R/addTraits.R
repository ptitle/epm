##' @title addTraits
##'
##' @description Add univariate or multivariate trait data to an epmGrid object.
##'
##' @param x object of class \code{epmGrid}
##' @param data named numeric vector, matrix or dataframe with rownames
##'   corresponding to species in \code{x} or pairwise matrix with row and
##'   column names corresponding to species in \code{x}. If pairwise matrix, the
##'   upper triangle of the matrix will be used for calculations.
##' @param replace boolean; if data is already a part of \code{x}, should it be
##'   replaced?
##' @param verbose if TRUE, list out all species that are dropped/excluded,
##'   rather than counts.
##'
##' @details If any species in \code{data} are not found in the epmGrid
##'   geographical data, then those species will be dropped from \code{data},
##'   and a warning will be issued.
##'
##' @return object of class \code{epmGrid}, with trait data as the list element
##'   named \code{data}.
##'
##' @author Pascal Title
##'
##' @examples
##' tamiasEPM
##' tamiasTraits
##'
##' addTraits(tamiasEPM, tamiasTraits)
##'
##' @export


addTraits <- function(x, data, replace = FALSE, verbose = FALSE) {
	
	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	if (!inherits(data, c('numeric', 'matrix', 'data.frame'))) {
		stop('data must be either a numeric vector, matrix or dataframe.')
	}
	
	if (inherits(x[['data']], c('numeric', 'matrix', 'data.frame')) & !replace) {
		stop('Data already present. If data are to be replaced, set replace = TRUE')
	}
	
	# drop species from trait vector if missing from grid
	if (is.null(dim(data))) {
		if (is.null(names(data))) {
			stop('Data must have names.')
		}
		data <- data[!is.na(data)]
		traitSpecies <- intersect(x$geogSpecies, names(data))
		inGeogNotData <- setdiff(x$geogSpecies, names(data))
		inDataNotGeog <- setdiff(names(data), x$geogSpecies)
		if (length(traitSpecies) == 0) {
			stop('There are no common species in geographic and trait data.')
		}
		x[['data']] <- data[traitSpecies]
	} else if (class(data) %in% c('matrix','data.frame')) {
		if (is.null(rownames(data))) {
			stop('Data must have rownames.')
		}
		traitSpecies <- intersect(x$geogSpecies, rownames(data))
		inGeogNotData <- setdiff(x$geogSpecies, rownames(data))
		inDataNotGeog <- setdiff(rownames(data), x$geogSpecies)
		if (length(traitSpecies) == 0) {
			stop('There are no common species in geographic and trait data.')
		}
		
		if (identical(rownames(data), colnames(data))) {
			# pairwise matrix
			x[['data']] <- as.data.frame(data[traitSpecies, traitSpecies], stringsAsFactors = FALSE)
		} else {
			x[['data']] <- as.data.frame(data[traitSpecies,], stringsAsFactors = FALSE)
		}
	}
	
	if (length(inDataNotGeog) > 0) {
	    if (verbose) {
		    msg <- paste0('The following species were dropped from the trait data because they lacked geographic data:\n\t', paste(inDataNotGeog, collapse='\n\t'))
	    } else {
	        msg <- paste0(length(inDataNotGeog), ' species ', ifelse(length(inDataNotGeog) > 1, 'were', 'was'), ' dropped from the trait data because ', ifelse(length(inDataNotGeog) > 1, 'they lack', 'it lacks'), ' geographic data.\n')
	    }
		warning(msg)
	}
	
	return(x)
}

