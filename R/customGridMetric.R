##'@title Custom grid metrics
##'
##'@description Define your own function for summarizing information across grid
##'  cells.
##'
##'@param x object of class \code{epmGrid}
##'@param fun a function to apply to all grid cells (see details)
##'@param column If a univariate morphological metric is specified, and the data
##'  in \code{x} are multivariate, which trait should be used? This can also
##'  specify which subset of columns a multivariate metric should be applied to.
##'@param minTaxCount the minimum number of taxa needed to apply the function.
##'  For instance, should the function be applied to gridcells with just 1
##'  taxon?
##'@param metricName the name you would like to attach to the output
##'
##'@details This function allows you to not be limited to the diversity metrics
##'  available via the \code{\link{gridMetrics}} function. \cr
##'
##'  The custom function should have just one input: a vector of taxon names
##'  that will then be used to subset the trait or phylogenetic data. Within the
##'  function call, the trait data already attached to the epmGrid object must
##'  be referred to as dat, and the phylogenetic tree already attached to the
##'  epmGrid must be referred to as phylo.\cr 
##'
##'  If the input epmGrid object contains a set of trees, then this function will
##' 	 be applied, using each tree in turn, and will return a list of results. This
##'  list can then be passed to \code{\link{summarizeEpmGridList}} to be summarized.
##'  
##'  See examples below.
##'
##'@return object of class \code{epmGrid}, or list of \code{epmGrid} objects
##'
##'@author Pascal Title
##'
##' @examples
##' tamiasEPM
##' tamiasEPM <- addPhylo(tamiasEPM, tamiasTree)
##' tamiasEPM <- addTraits(tamiasEPM, tamiasTraits)
##'
##' # In the following examples, notice that any mention of the trait data or 
##' ## phylogeny that are already attached to the epmGrid object are referred 
##' ## to as dat and phylo.
##'
##' # example: calculate morphological disparity 
##' ## (already implemented in gridMetrics)
##' f <- function(cells) {
##' 	sum(diag(cov(dat[cells,])))
##' }
##'
##' # to calculate disparity, we need at least 2 taxa
##'
##' xx <- customGridMetric(tamiasEPM, fun = f, minTaxCount = 2, 
##' metricName = 'disparity')
##'
##' # In the example above, gridcells with 1 species are left as NA. 
##' ## But if we wanted those gridcells to have a value of 0 rather than NA, 
##' ## we could do the following:
##' f <- function(sp) {
##' 	if (length(sp) == 1) {
##' 		0
##' 	} else {
##' 		sum(diag(cov(dat[sp,])))
##' 	}
##' }
##'
##' # and change minTaxCount to 1
##' xx <- customGridMetric(tamiasEPM, fun = f, minTaxCount = 1, 
##' metricName = 'disparity')
##'
##'
##' # phylogenetic example: mean patristic distance
##' ## this example doesn't actually involve the phylogeny internally,
##' ## we can just supply what is needed to the function
##' patdist <- cophenetic(tamiasEPM[['phylo']])
##' patdist[upper.tri(patdist, diag = TRUE)] <- NA
##' f <- function(cells) {
##' 	mean(patdist[cells, cells], na.rm = TRUE)
##' }
##'
##' xx <- customGridMetric(tamiasEPM, fun = f, minTaxCount = 1, 
##' metricName = 'mean patristic')
##'
##' # an example that involves both morphological and phylogenetic data
##' ## nonsensical, but for illustrative purposes:
##' ## ratio of Faith's phylogenetic diversity to morphological range
##' f <- function(cells) {
##' 	faithPD(phylo, cells) / max(dist(dat[cells, ]))
##' }
##'
##' xx <- customGridMetric(tamiasEPM, fun = f, minTaxCount = 2, 
##' metricName = 'PD_range_ratio')
##'
##'
##' # Example involving a set of trees
##' tamiasEPM <- addPhylo(tamiasEPM, tamiasTreeSet, replace = TRUE)
##' 
##' # get crown clade age of clade containing taxa present in grid cell
##' f <- function(sp) {
##' 	ape::branching.times(phylo)[as.character(ape::getMRCA(phylo, sp))]
##' }
##' 
##' xx <- customGridMetric(tamiasEPM, fun = f, minTaxCount = 2, metric = 'nodeAge')
##' 
##' 
##'
##'@export



# where phenotypic or other species-specific data = dat
# where phylogeny is phylo
# where the only input variable can be the grid cell taxa
customGridMetric <- function(x, fun, column = NULL, minTaxCount = 1, metricName = 'custom_metric') {

	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	if (inherits(x[['phylo']], 'phylo')) {
		x[['phylo']] <- list(x[['phylo']])
		class(x[['phylo']]) <- 'multiPhylo'
	}
	
	# what is being called in the function fun?
	useDat <- FALSE
	usePhylo <- FALSE
	if (any(grepl('dat', as.character(body(fun))))) {
		useDat <- TRUE
	}
	if (any(grepl('phylo', as.character(body(fun))))) {
		usePhylo <- TRUE
	}
	if (length(names(formals(fun))) != 1) {
		stop('function should only have one input argument.')
	}
		
	if (useDat) {
		if (is.null(x[['data']])) {
			stop('epmGrid object does not contain trait data!')
		}
	}
	
	if (usePhylo) {
		if (is.null(x[['phylo']])) {
			stop('epmGrid object does not contain a phylo object!')
		}
	}
	
	pairwise <- FALSE
	
	# if data is pairwise matrix, then set flag appropriately
	if (inherits(x[['data']], c('matrix', 'data.frame'))) {
		if (identical(rownames(x[['data']]), colnames(x[['data']]))) {
			message('\t...detected pairwise distance matrix...\n') 
			column <- NULL
			pairwise <- TRUE
			# make the diagonal and lower triangle NA
			diag(x[['data']]) <- NA
			# x[['data']][lower.tri(x[['data']], diag = TRUE)] <- NA
			if (all(is.na(x[['data']][upper.tri(x[['data']])]))) {
				stop('There are no values in the upper triangle of the pairwise matrix.')
			}
		}
	}
	
	# if column is defined but data are univariate, disable column. 
	if (is.vector(x[['data']]) & !is.null(column)) {
		column <- NULL
	}
		
	if (!is.null(column) & inherits(x[['data']], c('matrix', 'data.frame'))) {
		if (is.character(column)) {
			if (any(!column %in% colnames(x[['data']]))) {
				stop('column not a valid column name of the data.')
			}
		} else if (is.numeric(column)) {
			if (any(!all(column %in% 1:ncol(x[['data']])))) {
				stop('Requested data column indices are out of range.')
			}
		}
	}
	
	# if (inherits(x[['data']], c('matrix', 'data.frame')) | is.vector(x[['data']])) {
		# if (inherits(x[['data']], c('matrix', 'data.frame'))) {
			# metricType <- 'multiVar'
			# if (!is.null(column) & length(column) == 1) {
				# metricType <- 'uniVar'
			# }
		# } else {
			# metricType <- 'uniVar'
		# }
	# } else {
		# metricType <- 'none'
	# }
		
	# if a subset of data columns are requested, subset the data table
	if (!is.null(column) & inherits(x[['data']], c('matrix', 'data.frame'))) {
		if (length(column) > 1) {
			x[['data']] <- x[['data']][, column]
		} else {
			x[['data']] <- setNames(x[['data']][, column], rownames(x[['data']]))
		}
		column <- NULL
	}
		
	# prune epmGrid
	if (useDat & !usePhylo) {
		## if only traits are involved	
		# prune epmGrid object down to species shared with data
		if (is.vector(x[['data']])) {
			x[['speciesList']] <- intersectList(x[['speciesList']], names(x[['data']]))
		} else {
			x[['speciesList']] <- intersectList(x[['speciesList']], rownames(x[['data']]))
		}
	}
	
	if (!useDat & usePhylo) {
	 	# prune epmGrid object down to species shared with phylogeny
		x[['speciesList']] <- intersectList(x[['speciesList']], x[['phylo']][[1]]$tip.label)
	}
	
	if (useDat & usePhylo) {
		if (is.vector(x[['data']])) {
			sharedTaxa <- intersect(names(x[['data']]), x[['phylo']][[1]]$tip.label)
			x[['data']] <- x[['data']][sharedTaxa]
		} else {
			sharedTaxa <- intersect(rownames(x[['data']]), x[['phylo']][[1]]$tip.label)
			x[['data']] <- x[['data']][sharedTaxa,]
		}
		x[['phylo']] <- lapply(x[['phylo']], ape::keep.tip, sharedTaxa)
		x[['speciesList']] <- intersectList(x[['speciesList']], sharedTaxa)
	}

	# get unique gridcell communities
	uniqueComm <- x[['speciesList']]
	
	############################
	
	if (usePhylo) {
		if (length(x[['phylo']]) == 1) {
			
			# create custom environment for function and add in trait and phylo data
			.epm.env <- new.env()
			assign('dat', x[['data']], envir = .epm.env)
			assign('phylo', x[['phylo']][[1]], envir = .epm.env)
			# ls(envir = .epm.env)
			environment(fun) <- .epm.env
						
			ret <- calcCustomGridMetric(x, uniqueComm, fun, minTaxCount, metricName)

		} else {
			# multiPhylo
			op <- getOption("pboptions")
			opb <- pbapply::pboptions(type = 'none')
			pb <- utils::txtProgressBar(max = length(x[['phylo']]), char = '+', style = 3)
			ret <- vector('list', length(x[['phylo']]))
			for (i in 1:length(x[['phylo']])) {
				utils::setTxtProgressBar(pb, i)
				
				.epm.env <- new.env()
				assign('dat', x[['data']], envir = .epm.env)
				assign('phylo', x[['phylo']][[i]], envir = .epm.env)
				# ls(envir = .epm.env)
				environment(fun) <- .epm.env
				
				ret[[i]] <- calcCustomGridMetric(x, uniqueComm, fun, minTaxCount, metricName)
			}
			close(pb)
			opb <- pbapply::pboptions(op)

		}
		
	} else {
		
		# create custom environment for function and add in trait and phylo data
		.epm.env <- new.env()
		assign('dat', x[['data']], envir = .epm.env)
		assign('phylo', x[['phylo']][[1]], envir = .epm.env)
		# ls(envir = .epm.env)
		environment(fun) <- .epm.env
					
		ret <- calcCustomGridMetric(x, uniqueComm, fun, minTaxCount, metricName)
		
	}
	
	return(ret)
}
	
	
	
calcCustomGridMetric <- function(x, uniqueComm, fun, minTaxCount, metricName) {

	############################
	# Apply function to grid cells
	
	# initialize with all NA
	resVal <- rep(NA, length(uniqueComm))
	
	# only apply function to cells with more than minTaxCount taxa
	ind <- which((lengths(uniqueComm) >= minTaxCount & !sapply(uniqueComm, anyNA)) == TRUE)
	
	resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) fun(y))
	
	# for debugging
	# for (i in 1:length(ind)) {
		# f(uniqueComm[[ind[i]]])
	# }
	
	#######################
	# Assemble results
	
	resVal[is.nan(resVal)] <- NA
	cellVec <- numeric(length = length(x[['cellCommInd']]))
	for (i in 1:length(resVal)) {
		cellVec[x[['cellCommInd']] == i] <- resVal[i]
	}
	if (inherits(x[[1]], 'sf')) {
		x[[1]][metricName] <- cellVec
	} else {
		tmp <- terra::rast(x[[1]][[1]])
		tmp[] <- cellVec
		names(tmp) <- metricName
		x[[1]] <- c(x[[1]], tmp)
	}
	attributes(x)$metric <- metricName
	
	# update species richness, as some species may have been dropped for this metric calculation.
	if (inherits(x[[1]], 'sf')) {
		for (i in 1:length(x[['speciesList']])) {
			x[['grid']][which(x[['cellCommInd']] == i), 'spRichness'] <- length(x[['speciesList']][[i]])
		}
	} else {
		tmpVec <- terra::values(x[['grid']][['spRichness']])
		for (i in 1:length(x[['speciesList']])) {
			tmpVec[which(x[['cellCommInd']] == i)] <- length(x[['speciesList']][[i]])
		}
		terra::values(x[['grid']][['spRichness']]) <- tmpVec
		rm(tmpVec)
	}
		
	# update geog species
	x[['geogSpecies']] <- sort(unique(unlist(x[['geogSpecies']])))
	
	return(x)	
}


