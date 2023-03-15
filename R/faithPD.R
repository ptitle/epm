##' @title Calculate Faith's Phylogenetic Diversity
##'
##' @description Calculates Faith's PD for a specific set of tips
##'
##' @param phy phylogeny of class \code{phylo}
##' @param tips tip names to be included
##'
##' @return numeric value of summed phylogenetic diversity
##'
##' @details Returns the sum of total branch lengths that unite a 
##' set of species. The root is always included in these calculations. 
##' If tip is just one species, then the root-to-tip distance is returned.
##' 
##' @author Pascal Title
##'
##' @references
##' Faith D.P. (1992) Conservation evaluation and phylogenetic diversity. 
##' Biological Conservation, 61, 1-10.
##' 
##' @examples
##' tamiasTree
##'	faithPD(tamiasTree, c('Tamias_minimus', 'Tamias_speciosus'))
##' 
##' @export



faithPD <- function(phy, tips) {
	
	spEdges <- ape::nodepath(phy)
	names(spEdges) <- phy$tip.label
	spEdges <- spEdges[tips]
	spEdges <- lapply(spEdges, function(x) setdiff(match(x, phy$edge[,2]), NA))
	PDedges <- unique(unlist(spEdges))
	return(sum(phy$edge.length[PDedges]))	
}


# # for testing purposes, compare to picante::pd
# picantePD <- function(phy, tips) {
	# mat <- matrix(0, nrow = 1, ncol = length(phy$tip.label))
	# colnames(mat) <- phy$tip.label
	# mat[1, colnames(mat) %in% tips] <- 1

	# picante::pd(mat, phy, include.root = TRUE)$PD
# }






