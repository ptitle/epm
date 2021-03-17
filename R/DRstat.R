##' @title Calculate the DR statistic
##'
##' @description Calculates the tip-specific DR statistic for speciation rates
##'
##' @param tree phylogeny of class \code{phylo}
##'
##' @return named numeric vector of speciation rates
##'
##' @author Pascal Title
##'
##' @references
##' Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O. (2012). 
##' 	The global diversity of birds in space and time. Nature, 491, 444–448.
##' 
##' @examples
##' tamiasTree
##'	DRstat(tamiasTree)
##' 
##' @export

DRstat <- function(tree) {
	
	# ape::nodepath returns the nodes that make up the path, including tip nodes
	# We will drop the first node so that each listed node is the tip-side node of each branch
	spEdges <- ape::nodepath(tree)
	spEdges <- lapply(spEdges, function(x) rev(x[- 1]))
	# spEdges <- lapply(spEdges, function(x) sapply(x, function(y) tree$edge.length[tree$edge[, 2] == y]))
	edge2 <- tree$edge[, 2]
	bl <- tree$edge.length
	spEdges <- lapply(spEdges, function(x) rev(bl[edge2 %in% x]))
	
	rates <- sapply(spEdges, function(x) sum(sapply(1:length(x), function(y) x[y] * (1 / (2 ^ (y - 1))))))
	rates <- rates ^ (-1)
	names(rates) <- tree$tip.label
	
	return(rates)
}
