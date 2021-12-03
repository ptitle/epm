##' @title Subset epmGrid to shared taxa
##'
##' @description An EpmGrid object may contain more taxa with morphological
##' 	data than taxa with phylogenetic information, or vice versa. This 
##' 	function subsets all epmGrid components to the set of taxa shared 
##' 	across geographic, phenotypic and phylogenetic datasets. This might
##' 	desirable to ensure that all diversity metrics are based on the same 
##' 	set of taxa.
##'
##' @param x object of class \code{epmGrid}
##'
##'
##' @return new \code{epmGrid} object.
##'
##' @author Pascal Title
##' 
##' @examples
##' tamiasEPM
##' # randomly drop a few species for demonstration
##' tamiasEPM <- addPhylo(tamiasEPM, ape::drop.tip(tamiasTree, sample(tamiasTree$tip.label, 5)))
##' tamiasEPM <- addTraits(tamiasEPM, tamiasTraits[-(3:5),])
##'
##' new <- reduceToCommonTaxa(tamiasEPM)
##'
##' tamiasEPM
##' new
##'	
##' 
##' @export


reduceToCommonTaxa <- function(x) {
	
	if (!inherits(x, 'epmGrid')) {
		stop('Object must be of class epmGrid')
	}
	
	if (is.null(x[['data']]) & !inherits(x[['phylo']], 'phylo')) {
		stop('This epmGrid only contains geographic data.')
	}
	
	set1 <- x$geogSpecies

	if (!is.null(x[['data']])) {
		if (is.vector(x[['data']])) {
			set2 <- names(x[['data']])
		} else {
			set2 <- rownames(x[['data']])
		}
	} else {
		set2 <- set1
	}
	
	if (inherits(x[['phylo']], 'phylo')) {
		set3 <- x[['phylo']]$tip.label
	} else {
		set3 <- set1
	}
	
	# These are the taxa shared across all components
	commontaxa <- Reduce(intersect, list(set1, set2, set3))
	
	# Now subset each component to these taxa
	newEpm <- dropSpecies(x, setdiff(set1, commontaxa))
	
	if (!is.null(x[['data']])) {
		if (is.vector(x[['data']])) {
			newEpm[['data']] <- x[['data']][commontaxa]
		} else {
			newEpm[['data']] <- x[['data']][commontaxa, ]
		}
	}
	
	if (inherits(x[['phylo']], 'phylo')) {
		newEpm[['phylo']] <- ape::keep.tip(x[['phylo']], commontaxa)
	}
	
	return(newEpm)
	
}