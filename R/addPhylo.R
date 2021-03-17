##' @title addPhylo
##'
##' @description Add a phylogeny to epmGrid object.
##'
##' @param x object of class \code{epmGrid}
##' @param tree a phylogeny of class \code{phylo}
##' @param replace boolean; if a tree is already a part of \code{x},
##' should it be replaced?
##'
##' @details If any species in the phylogeny are not found in the epmGrid
##' geographical data, then those species will be dropped from the phylogeny, and
##' a warning will be issued. 
##'
##' @return object of class \code{epmGrid}, with a \code{phylo}
##' object as the list element named \code{phylo}. 
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

addPhylo <- function(x, tree, replace = FALSE) {
	
	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	if (!inherits(tree, 'phylo')) {
		stop('tree must be a phylo object.')
	}
	
	if (inherits(x[['phylo']], 'phylo') & !replace) {
		stop('Phylogeny already present. If phylogeny is to be replaced, set replace = TRUE')
	}
	
	# if needed, prune species in phylogeny down to species with geog data
	if (length(intersect(tree$tip.label, x[['geogSpecies']])) == 0) {
		stop('No matching taxa between geographic data and phylogeny.')
	}
	inPhyloNotGeog <- setdiff(tree$tip.label, x[['geogSpecies']])
	inGeogNotPhylo <- setdiff(x[['geogSpecies']], tree$tip.label)
	tree <- ape::drop.tip(tree, inPhyloNotGeog)
	x[['phylo']] <- tree
	
	if (length(inPhyloNotGeog) > 0) {
		msg <- paste0('The following species were pruned from the phylogeny because they lack geographic data:\n\t', paste(inPhyloNotGeog, collapse='\n\t'))
		warning(msg)
	}
	
	return(x)
}

