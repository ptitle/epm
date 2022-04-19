##' @title addPhylo
##'
##' @description Add a phylogeny to epmGrid object.
##'
##' @param x object of class \code{epmGrid}
##' @param tree a phylogeny of class \code{phylo}, or a set of trees of class \code{multiPhylo}
##' @param replace boolean; if a tree is already a part of \code{x}, should it
##'   be replaced?
##' @param verbose if TRUE, list out all species that are dropped/excluded,
##'   rather than counts.
##'
##' @details If any species in the phylogeny are not found in the epmGrid
##'   geographical data, then those species will be dropped from the phylogeny,
##'   and a warning will be issued.
##' 
##' 	 If providing a set of trees as a multiPhylo object, it is expected that all
##' 	  trees have the same tips. 
##'
##' @return object of class \code{epmGrid}, with a \code{phylo} object as the
##'   list element named \code{phylo}.
##'
##' @author Pascal Title
##'
##' @examples
##' tamiasEPM
##' tamiasTree
##'
##' addPhylo(tamiasEPM, tamiasTree)
##'
##' @export

addPhylo <- function(x, tree, replace = FALSE, verbose = FALSE) {
	
	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	if (!inherits(tree, c('phylo', 'multiPhylo'))) {
		stop('tree must be a phylo or multiPhylo object.')
	}
	
	if (inherits(x[['phylo']], c('phylo', 'multiPhylo')) & !replace) {
		stop('Phylogeny already present. If phylogeny is to be replaced, set replace = TRUE')
	}
	
	# to simplify the code, if input is a single tree, we will convert to a single-tree multiPhylo
	# object. Then no need for code duplication.
	if (inherits(tree, 'phylo')) {
		tree <- list(tree)
		class(tree) <- 'multiPhylo'
	}
	
	# if multiPhylo, then all trees must have the same tips. 
	if (inherits(tree, 'multiPhylo')) {
		alltips <- lapply(tree, function(x) x$tip.label)
		if (!all(sapply(alltips, function(y) length(intersect(alltips[[1]], y)) / length(alltips[[1]])) == 1)) {
			stop('All trees in multiPhylo are expected to have the same tips.')
		}
	}
	
	# if needed, prune species in phylogeny down to species with geog data
	if (length(intersect(tree[[1]]$tip.label, x[['geogSpecies']])) == 0) {
		stop('No matching taxa between geographic data and phylogeny.')
	}		

	inPhyloNotGeog <- setdiff(tree[[1]]$tip.label, x[['geogSpecies']])
	inGeogNotPhylo <- setdiff(x[['geogSpecies']], tree[[1]]$tip.label)
	tree <- lapply(tree, function(y) ape::drop.tip(y, inPhyloNotGeog))
	class(tree) <- 'multiPhylo'
		
	if (length(tree) == 1) {
		tree <- tree[[1]]
	}

	x[['phylo']] <- tree
	
	if (length(inPhyloNotGeog) > 0) {
	    if (verbose) {
		    msg <- paste0('The following species were pruned from the phylogeny because they lack geographic data:\n\t', paste(inPhyloNotGeog, collapse='\n\t'))
    	} else {
    	    msg <- paste0(length(inPhyloNotGeog), ' species ', ifelse(length(inPhyloNotGeog) > 1, 'were', 'was'), ' pruned from the phylogeny because ', ifelse(length(inPhyloNotGeog) > 1, 'they lack', 'it lacks'), ' geographic data.\n')
    	}
	    warning(msg)
	}
	
	return(x)
}

