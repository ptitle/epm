##'@title Grid Metrics
##'
##'@description Calculate various morphological and phylogenetic community
##'  metrics for every cell in a \code{epmGrid} object. To implement other 
##'  metrics not available here, see \code{\link{customGridMetric}}.
##'
##'@param x object of class \code{epmGrid}
##'
##'@param metric name of metric to use, see Details.
##'
##'@param column If a univariate morphological metric is specified, and the data
##'  in \code{x} are multivariate, which trait should be used? This can also
##'  specify which subset of columns a multivariate metric should be applied to.
##'
##'
##'@param verbose print various messages to the console. Default is TRUE.
##'
##'@return object of class \code{epmGrid} where the grid represents calculations
##'  of the metric at every cell. The species identities per grid cell are those
##'  that had data for the calculation of the metric. If taxa were dropped from
##'  the initial epmGrid object, then they have been removed from this epmGrid.
##'	 If a set of trees was involved, then returns a list of \code{epmGrid} objects.
##'
##'@details Univariate trait metrics 
##'\itemize{ 
##'    \item{mean} 
##'    \item{median}
##'    \item{range} 
##'    \item{variance}
##'    \item{mean_NN_dist:} {mean nearest neighbor distance}
##'    \item{min_NN_dist:} {minimum nearest neighbor distance}
##'	   \item{evenness:} {variance of nearest neighbor distances, 
##' 			larger values imply decreasing evenness}
##'    \item{arithmeticWeightedMean} (see below) 
##'    \item{geometricWeightedMean} (see below) 
##'}
##'Multivariate trait metrics 
##'\itemize{ 
##'    \item{disparity}
##'    \item{partialDisparity:} {contribution of species in each gridcell to 
##'    overall disparity, returned as the ratio of summed partial disparities
##'    to total disparity.} 
##'    \item{range} 
##'    \item{mean_NN_dist:} {mean nearest neighbor distance} 
##'    \item{min_NN_dist:} {minimum nearest neighbor distance}
##'	   \item{evenness:} {variance of nearest neighbor distances, 
##' 		larger values imply decreasing evenness}
##' } 
##' Phylogenetic metrics 
##' \itemize{ 
##'     \item{pd:} {Faith's phylogenetic diversity, including the root} 
##'     \item{meanPatristic}
##'     \item{meanPatristicNN:} {mean nearest neighbor in patristic distance}
##'     \item{minPatristicNN:} {minimum nearest neighbor in patristic distance}
##'	    \item{phyloEvenness:} {variance of nearest neighbor patristic distances,
##' 		larger values imply decreasing evenness}
##'     \item{phyloDisparity:} {sum of squared deviations in patristic distance}
##'     \item{PSV:} {Phylogenetic Species Variability} 
##'     \item{PSR:} {Phylogenetic Species Richness} 
##'     \item{DR:} {non-parametric estimate of speciation rates} 
##' }
##'Range-weighted metrics 
##'\itemize{ 
##'     \item{weightedEndemism:} {Species richness inversely weighted by range size.}
##'     \item{correctedWeightedEndemism:} {Weighted endemism standardized by 
##'     species richness} 
##'     \item{phyloWeightedEndemism:} {Phylogenetic diversity inversely weighted 
##' 			by range size associated with each phylogenetic branch.} 
##' }
##'
##'If data slot contains a pairwise matrix, \code{column} is ignored. Weighted
##'mean options are available where, for each cell, a weighting scheme (inverse
##'of species range sizes) is applied such that small-ranged species are
##'up-weighted, and broadly distributed species are down-weighted. This can be a
##'useful way to lessen the influence of broadly distributed species in the
##'geographic mapping of trait data.
##'
##'It may be desirable to have metrics calculated for a dataset where only taxa
##'shared across geography, traits and phylogeny are included. The function
##'\code{\link{reduceToCommonTaxa}} does exactly that.
##'
##'If a set of trees are associated with the input epmGrid object \code{x}, 
##' then the metric is calculated for each tree, and a list of epmGrid objects 
##' is returned. This resulting list can be summarized with the function
##' \code{\link{summarizeEpmGridList}}. For instance the mean and variance can 
##' be calculated, to show the central tendency of the metric across grid cells, 
##' and to quantify where across geography variability in phylogenetic topography
##' manifests itself.  
##'
##'To implement other metrics not available here, see
##'\code{\link{customGridMetric}}.
##'
##'
##'@references 
##'partial disparity\cr Foote, M. (1993). Contributions of individual taxa to
##'overall morphological disparity. Paleobiology, 19(4), 403–419.
##'https://doi.org/10.1017/s0094837300014056
##'
##'PSV, RSV\cr Helmus, M. R., Bland, T. J., Williams, C. K., & Ives, A. R.
##'(2007). Phylogenetic Measures of Biodiversity. The American Naturalist,
##'169(3), E68–E83. https://doi.org/10.1086/511334
##'
##'DR\cr Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O.
##'(2012). The global diversity of birds in space and time. Nature, 491(7424),
##'444–448. https://doi.org/10.1038/nature11631
##'
##'weighted endemism\cr Crisp, M. D., Laffan, S., Linder, H. P., & Monro, A.
##'(2001). Endemism in the Australian flora. Journal of Biogeography, 28(2),
##'183–198. https://doi.org/10.1046/j.1365-2699.2001.00524.x
##'
##'phylo weighted endemism\cr Rosauer, D., Laffan, S. W., Crisp, M. D.,
##'Donnellan, S. C., & Cook, L. G. (2009). Phylogenetic endemism: a new approach
##'for identifying geographical concentrations of evolutionary history.
##'Molecular Ecology, 18(19), 4061–4072.
##'https://doi.org/10.1111/j.1365-294x.2009.04311.x
##'
##' @examples
##' tamiasEPM <- addPhylo(tamiasEPM, tamiasTree)
##' tamiasEPM <- addTraits(tamiasEPM, tamiasTraits)
##'
##' # univariate morphological example
##' x <- gridMetrics(tamiasEPM, metric='mean', column='V2')
##' plot(x)
##' \donttest{
##' # multivariate morphological
##' x <- gridMetrics(tamiasEPM, metric='disparity')
##' plot(x)
##'
##' # phylogenetic metrics
##' x <- gridMetrics(tamiasEPM, metric='meanPatristic')
##' plot(x)
##' }
##'@export


gridMetrics <- function(x, metric, column = NULL, verbose = TRUE) {
	
	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	if (!identical(names(x), c('grid', 'speciesList', 'cellCommInd', 'geogSpecies', 'cellCount', 'data', 'phylo'))) {
		stop('enmGrid object is missing some components.')
	}

	if (length(metric) > 1) {
		stop('You can only specify one metric.')
	}
	
	if (inherits(x[['phylo']], 'phylo')) {
		x[['phylo']] <- list(x[['phylo']])
		class(x[['phylo']]) <- 'multiPhylo'
	}
	
	metric <- match.arg(metric, choices = c('mean', 'median', 'range', 'variance', 'evenness', 'arithmeticWeightedMean', 'geometricWeightedMean', 'rangePCA', 'disparity', 'mean_NN_dist', 'min_NN_dist', 'pd', 'meanPatristic', 'meanPatristicNN', 'minPatristicNN', 'phyloEvenness', 'phyloDisparity', 'PSV', 'PSR', 'DR', 'weightedEndemism', 'correctedWeightedEndemism', 'phyloWeightedEndemism', 'partialDisparity'))
	
	pairwise <- FALSE
	phyloNeeded <- FALSE
	
	# if data is pairwise matrix, then set flag appropriately
	if (inherits(x[['data']], c('matrix', 'data.frame'))) {
		if (identical(rownames(x[['data']]), colnames(x[['data']]))) {
			if (verbose) message('\t...detected pairwise distance matrix...\n') 
			column <- NULL
			pairwise <- TRUE
			# make the diagonal NA
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
	
	if (inherits(x[['data']], c('matrix', 'data.frame')) | is.vector(x[['data']])) {
		if (inherits(x[['data']], c('matrix', 'data.frame'))) {
			metricType <- 'multiVar'
			if (!is.null(column) & length(column) == 1) {
				metricType <- 'uniVar'
			}
		} else {
			metricType <- 'uniVar'
		}
	} else {
		metricType <- 'none'
	}
	
	# if a subset of data columns are requested, subset the data table
	if (!is.null(column) & inherits(x[['data']], c('matrix', 'data.frame'))) {
		if (length(column) > 1) {
			x[['data']] <- x[['data']][, column]
		} else {
			x[['data']] <- setNames(x[['data']][, column], rownames(x[['data']]))
		}
		column <- NULL
	}
	
	if (metric %in% c('mean', 'median', 'variance', 'arithmeticWeightedMean', 'geometricWeightedMean') & metricType == 'multiVar' & !pairwise) {
		stop('If a univariate metric is requested from a multivariate dataset, a column name must be provided as column.')
	}
	
	if (metric %in% c('rangePCA', 'disparity', 'partialDisparity') & metricType == 'uniVar' & !pairwise) {
		stop('A multivariate metric cannot be applied to a univariate dataset.')
	}
	
	if (metric == 'geometricWeightedMean') {
		if (any(x[['data']] < 0)) {
			stop('Negative trait values cannot be supplied to geometricWeightedMean because the log of a negative number is NA.')
		}
	}
	
	# For metrics that require both trait and phylogenetic data, prune the data and tree to the shared taxon set
	if (metric == 'phylosignal') {
		phyloNeeded <- TRUE
		if (is.vector(x[['data']])) {
			sharedTaxa <- intersect(names(x[['data']]), x[['phylo']][[1]]$tip.label)
			x[['data']] <- x[['data']][sharedTaxa]
		} else {
			sharedTaxa <- intersect(rownames(x[['data']]), x[['phylo']][[1]]$tip.label)
			x[['data']] <- x[['data']][sharedTaxa,]
		}
		x[['phylo']] <- lapply(x[['phylo']], ape::keep.tip, sharedTaxa)
		class(x[['phylo']]) <- 'multiPhylo'
	}
	
	# Prune species list according to metric of interest
	if (metric %in% c('mean', 'median', 'disparity', 'partialDisparity', 'range', 'variance', 'evenness', 'arithmeticWeightedMean', 'geometricWeightedMean', 'rangePCA', 'mean_NN_dist', 'min_NN_dist', 'phylosignal')) {
		
		if (verbose) message('\t...dropping species that are not in trait data...\n')
		# check that there is data in epmGrid object
		if (is.null(x[['data']])) {
			stop('epmGrid object does not contain trait data!')
		}
		
		# prune epmGrid object down to species shared with data
		if (is.vector(x[['data']])) {
			x[['speciesList']] <- intersectList(x[['speciesList']], names(x[['data']]))
		} else {
			x[['speciesList']] <- intersectList(x[['speciesList']], rownames(x[['data']]))
		}
		
	 } else if (metric %in% c('pd', 'meanPatristic', 'meanPatristicNN', 'minPatristicNN', 'phyloEvenness', 'phyloDisparity', 'phyloWeightedEndemism', 'PSV', 'PSR', 'DR')) {
	 	phyloNeeded <- TRUE
	 	
	 	# check that there is a phylogeny in epmGrid object
		if (is.null(x[['phylo']])) {
			stop('epmGrid object does not contain a phylo object!')
		}
		
		if (verbose) message('\t...dropping species that are not in phylo data...\n')
	 	# prune epmGrid object down to species shared with phylogeny
		x[['speciesList']] <- intersectList(x[['speciesList']], x[['phylo']][[1]]$tip.label)
	
	} else if (!metric %in% c('weightedEndemism', 'correctedWeightedEndemism')) {
		stop('Metric not recognized!')
	}
	
	uniqueComm <- x[['speciesList']]
	
	if (phyloNeeded) {
		
		if (length(x[['phylo']]) == 1) {
			# single tree
			## still pass along data because might be needed for some metrics.
			ret <- calcGridMetric(x, uniqueComm, metric, tree = x[['phylo']][[1]], dat = x[['data']], metricType = metricType, pairwise = pairwise, verbose = verbose)
			
		} else {
			# multiPhylo
			# disable progress bars and enable progress bar that will show 
			# progress across trees instead.
			op <- getOption("pboptions")
			opb <- pbapply::pboptions(type = 'none')
			# getOption("pboptions")
			
			if (verbose) message('\t...calculating metric ', metric, ' over ', length(x[['phylo']]), ' trees...\n')
			
			pb <- utils::txtProgressBar(max = length(x[['phylo']]), char = '+', style = 3)
			ret <- vector('list', length(x[['phylo']]))
			for (i in 1:length(x[['phylo']])) {
				utils::setTxtProgressBar(pb, i)
				ret[[i]] <- calcGridMetric(x, uniqueComm, metric, tree = x[['phylo']][[i]], dat = x[['data']], metricType = metricType, pairwise = pairwise, verbose = FALSE)
			}
			close(pb)
			opb <- pbapply::pboptions(op)
		}
		
	
	} else {
		# traits only
		ret <- calcGridMetric(x, uniqueComm, metric, tree = NULL, dat = x[['data']], metricType = metricType, pairwise = pairwise, verbose = verbose)

	}
	
	return(ret)
}
	
	
	
	
	
	
	


calcGridMetric <- function(x, uniqueComm, metric, tree = NULL, dat = NULL, metricType, pairwise, verbose) {

	## ----------------------------------
	## MORPHOLOGY-RELATED METRICS
	
	## UNIVARIATE
	
	if (metric %in% c('mean', 'median', 'variance', 'range', 'mean_NN_dist', 'min_NN_dist', 'evenness') & !pairwise & metricType == 'uniVar') {
		if (verbose) message('\t...calculating univariate metric: ', metric, '...\n')
		resVal <- cellAvg(uniqueComm, trait = dat, stat = metric)
		
		# to be consistent with other metrics, we will have 'variance' return 0 for single-species cells
		if (metric == 'variance') {
			ind <- which(lengths(uniqueComm) == 1 & sapply(uniqueComm, anyNA) == FALSE)
			resVal[ind] <- 0
		}
	}
	
	if (metric %in% c('arithmeticWeightedMean', 'geometricWeightedMean') & metricType == 'uniVar') {
		if (verbose) message('\t...calculating univariate metric: ', metric, '...\n')
			
		# get inverse range sizes (1 / cell counts)
		inverseCellCount <- 1 / x[['cellCount']]
		
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(lengths(uniqueComm) > 0 & sapply(uniqueComm, anyNA) == FALSE)

		if (metric == 'arithmeticWeightedMean') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
				(1 / sum(inverseCellCount[y])) * sum(inverseCellCount[y] * dat[y])
			})
		}
		if (metric == 'geometricWeightedMean') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
				exp(sum(inverseCellCount[y] * log(dat[y])) / sum(inverseCellCount[y]))
			})			
		}	
	}
	
	if (metric %in% c('mean', 'median', 'variance', 'range', 'mean_NN_dist', 'min_NN_dist', 'evenness') & pairwise) {
	
		# if pairwise matrix
		if (verbose) message('\t...calculating pairwise metric: ', metric, '...\n')
		
		resVal <- rep(NA, length(uniqueComm)) # set up with NA
		
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, anyNA) == FALSE)

		if (metric == 'mean') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) mean(dat[y, y], na.rm = TRUE))
		} else if (metric == 'median') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) stats::median(dat[y, y], na.rm = TRUE))
		} else if (metric == 'variance') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) stats::var(as.numeric(dat[y, y]), na.rm = TRUE))
			ind <- which(lengths(uniqueComm) == 1 & sapply(uniqueComm, anyNA) == FALSE)
			resVal[ind] <- 0 # set single-sp cells to zero variance
		} else if (metric == 'mean_NN_dist') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) mean(apply(dat[y, y], 1, min, na.rm = TRUE)))
		} else if (metric == 'min_NN_dist') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) min(apply(dat[y, y], 1, min, na.rm = TRUE)))
		} else if (metric == 'evennness') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) stats::var(apply(dat[y, y], 1, min, na.rm = TRUE)))
		}
	}
	
	## MULTIVARIATE
	
	if (metric == 'disparity' & metricType == 'multiVar') {
		if (verbose) message('\t...calculating multivariate metric: ', metric, '...\n')
		# sum of the diagonal of the covariance matrix
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		# resVal <- rep(NA, length(uniqueComm)) # set up with NA
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(lengths(uniqueComm) > 1)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) sum(diag(cov(dat[y,]))))
	}
	
	if (metric == 'partialDisparity' & metricType == 'multiVar') {
		if (verbose) message('\t...calculating multivariate metric: ', metric, '...\n')
		partialDisp <- getSpPartialDisparities(dat)
		resVal <- pbapply::pbsapply(uniqueComm, function(y) sum(partialDisp[y]) / sum(partialDisp))
	}	
	
	if (metric == 'range' & metricType == 'multiVar') {
		if (verbose) message('\t...calculating multivariate metric: ', metric, '...\n')
		# maximum of the distance matrix (0 if one sp)
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(lengths(uniqueComm) > 1)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) max(dist(dat[y, ])))
	}
	
	if (metric == 'rangePCA' & metricType == 'multiVar') {
		if (verbose) message('\t...calculating multivariate metric: ', metric, '...\n')
		pc <- prcomp(dat)
		# retain 99% of the variation
		keep <- 1:which(cumsum(((pc$sdev)^2) / sum(pc$sdev^2)) >= 0.99)[1]
		if (length(keep) == 1) keep <- 1:ncol(pc$x)
		pc <- pc$x[, keep]
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(lengths(uniqueComm) > 1)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
			sum(apply(pc[y,], 2, function(z) diff(range(z))))
		})	
	}
	
	if (metric == 'mean_NN_dist' & metricType == 'multiVar') {
		if (verbose) message('\t...calculating multivariate metric: ', metric, '...\n')
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(lengths(uniqueComm) > 1)
		morphDists <- as.matrix(dist(dat))
		diag(morphDists) <- NA
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) mean(apply(morphDists[y, y], 1, min, na.rm = TRUE)))
	}

	if (metric == 'min_NN_dist' & metricType == 'multiVar') {
		if (verbose) message('\t...calculating multivariate metric: ', metric, '...\n')
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(lengths(uniqueComm) > 1)
		morphDists <- as.matrix(dist(dat))
		diag(morphDists) <- NA
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) min(apply(morphDists[y, y], 1, min, na.rm = TRUE)))
	}

	if (metric == 'evenness' & metricType == 'multiVar') {
		if (verbose) message('\t...calculating multivariate metric: ', metric, '...\n')
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(lengths(uniqueComm) > 1)
		morphDists <- as.matrix(dist(dat))
		diag(morphDists) <- NA
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) stats::var(apply(morphDists[y, y], 1, min, na.rm = TRUE)))
	}
	
	
	## ----------------------------------
	## PHYLOGENY-RELATED METRICS
	
	if (metric %in% c('pd', 'meanPatristic', 'minPatristicNN', 'meanPatristicNN', 'phyloEvenness', 'phyloDisparity', 'PSV', 'PSR', 'DR')) {
		
		# calculate pairwise patristic distance
		patdist <- cophenetic(tree)
		diag(patdist) <- NA
		
		if (metric == 'pd') {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			
			resVal <- rep(NA, length(uniqueComm)) # set up with NA
			ind <- which(sapply(uniqueComm, anyNA) == FALSE)		

			# this version is overly redundant because it generates nodepath each time, when it is only needed once. 
			# resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) faithPD(tree, y))
			
			spEdges <- ape::nodepath(tree)
			names(spEdges) <- tree$tip.label
				
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
				yy <- spEdges[y]
				yy <- lapply(yy, function(z) setdiff(match(z, tree$edge[,2]), NA))
				sum(tree$edge.length[unique(unlist(yy))])
			})
			
		}
		
		if (metric == 'meanPatristic') {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			# meanPatristic is 0 if 1 species, NA if no species
			patdist[upper.tri(patdist, diag = TRUE)] <- NA
			resVal <- pbapply::pbsapply(uniqueComm, function(y) mean(unlist(patdist[y, y]), na.rm = TRUE))
			resVal[lengths(uniqueComm) == 1] <- 0
			resVal[sapply(uniqueComm, anyNA)] <- NA
		}
		
		if (metric == 'meanPatristicNN') {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			# the mean of the minimum patristic distance for each species present
			resVal <- numeric(length = length(uniqueComm)) # set up with zeros
			resVal[sapply(uniqueComm, anyNA)] <- NA
			ind <- which(lengths(uniqueComm) > 1)
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
				return(mean(apply(patdist[y, y], MARGIN = 1, min, na.rm = TRUE)))
			})
		}

		if (metric == 'minPatristicNN') {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			# the minimum of the minimum patristic distance for each species present
			resVal <- numeric(length = length(uniqueComm)) # set up with zeros
			resVal[sapply(uniqueComm, anyNA)] <- NA
			ind <- which(lengths(uniqueComm) > 1)
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
				return(min(apply(patdist[y, y], MARGIN = 1, min, na.rm = TRUE)))
			})
		}

		if (metric == 'phyloEvenness') {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			# the variance of the minimum patristic distance for each species present
			resVal <- numeric(length = length(uniqueComm)) # set up with zeros
			resVal[sapply(uniqueComm, anyNA)] <- NA
			ind <- which(lengths(uniqueComm) > 1)
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
				return(stats::var(apply(patdist[y, y], MARGIN = 1, min, na.rm = TRUE)))
			})
		}

		if (metric == 'phyloDisparity') {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			# the sum of the squared deviations from the mean
			# value of 0 if 1 species, NA if no species
			
			## note: the mean patristic distance is the mean of pairwise distances, so the sample size is the number of unique pairs, not the number of species. Therefore, in calculating the variance, we will not divide by the number of species, but rather by the number of unique pairs (this is where choose(length(y), 2) comes from). 
			
			resVal <- numeric(length = length(uniqueComm)) # set up with zeros
			resVal[sapply(uniqueComm, anyNA)] <- NA
			ind <- which(lengths(uniqueComm) > 1)

			patdist[upper.tri(patdist, diag = TRUE)] <- NA
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) sum((patdist[y, y] - mean(patdist[y,y], na.rm = TRUE)) ^ 2, na.rm = TRUE) / choose(length(y), 2))
			
			# resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) sum(diag(ape::vcv.phylo(keep.tip(tree, y)))))
						
		}		
		
		if (metric %in% c('PSV', 'PSR')) {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			# Phylogenetic Species Variability is based on the variance-covariance matrix
			# measure of phylogenetic diversity, ranges from 0-1, not confounded by species richness
			
			vmat <- ape::vcv.phylo(tree, corr = TRUE)
			
			resVal <- numeric(length = length(uniqueComm)) # set up with zeros
			resVal[sapply(uniqueComm, anyNA)] <- NA
			ind <- which(lengths(uniqueComm) > 1)
			
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) (length(y) * sum(diag(vmat[y,y])) - sum(vmat[y,y]))/(length(y) * (length(y) - 1)))
		}
		
		if (metric == 'PSR') {
			# Phylogenetic Species Richness is richness * PSV
			# Conveys richness, but where the more phylogenetically distinct the species, the more weight it has. Will be less than or equal to richness.
			
			resVal[ind] <- resVal[ind] * lengths(uniqueComm[ind])
		}
		
		if (metric == 'DR') {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			# DR statistic for speciation rates
			
			tipDR <- DRstat(tree)			
			resVal <- pbapply::pbsapply(uniqueComm, function(y) mean(tipDR[y]))
		}		
	}
	
	## ----------------------------------
	## RANGE-WEIGHTED METRICS
	
	if (metric %in% c('weightedEndemism', 'correctedWeightedEndemism')) {
		if (verbose) message('\t...calculating weighted endemism metric...\n')
		
		# get inverse range sizes (1 / cell counts)
		inverseCellCount <- 1 / x[['cellCount']]

		resVal <- rep(NA, length(uniqueComm)) # set up with NA
		ind <- which(sapply(uniqueComm, anyNA) == FALSE)		
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) sum(inverseCellCount[y]))
	
		# if corrected metric, we will standardize by species richness
		if (metric == 'correctedWeightedEndemism') {
			resVal[ind] <- resVal[ind] / lengths(uniqueComm[ind])
		}
	}
	
	if (metric == 'phyloWeightedEndemism') {
		if (verbose) message('\t...calculating phylogenetic weighted endemism...\n')
		
		# Phylogenetic endemism is total branch length for edges uniting a set of taxa, standardized
		# by total tree branch length
		# divided by
		# sum of geographic area associated with those branches
		
		# list of descendant nodes from root to tip per taxon (includes root and terminal nodes)
		spNodes <- ape::nodepath(tree)
		names(spNodes) <- tree$tip.label
		
		# list of branch indices from root to tip
		spEdges <- lapply(spNodes, function(y) setdiff(match(y, tree$edge[,2]), NA))
		totalPD <- sum(tree$edge.length[unique(unlist(spEdges))])
		rootEdge <- tree$root.edge
		if (is.null(rootEdge)) {
			rootEdge <- 0
		}
		
		# get total geographic area covered by each branch
		allnodes <- 1:ape::Nnode(tree, internal.only = FALSE)
		
		fullSpList <- expandSpeciesCellList(x)
		
		# for each node, what is the geographic area (= number of grid cells)
		cellsPerNode <- integer(length(allnodes))
		for (i in 1:length(allnodes)) {
			nodeTips <- names(which(sapply(spNodes, function(y) allnodes[i] %in% y) == TRUE))
			cellsPerNode[i] <- length(which(sapply(fullSpList, function(y) any(nodeTips %in% y)) == TRUE))
		}
		
		resVal <- rep(NA, length(uniqueComm)) # set up with NA
		ind <- which(sapply(uniqueComm, anyNA) == FALSE)
	
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
			unitingEdges <- unique(unlist(spEdges[y]))
			unitingNodes <- unique(unlist(spNodes[y]))
			# relativePD <- sum(tree$edge.length[unitingEdges]) / totalPD
			
			# which descendant nodes belong to this set of branches?
			sum(c(rootEdge, tree$edge.length[unitingEdges] / 1) / cellsPerNode[unitingNodes])
		})
	
		# compare to:
		# commMat <- generateOccurrenceMatrix(x, sites = 'all')
		# testPE <- phyloregion::phylo_endemism(commMat, tree)
	}
		
	## ----------------------------------
	## REMAP RESULTS TO RELEVANT CELLS AND ASSIGN TO GRID
	
	resVal[is.nan(resVal)] <- NA
	cellVec <- numeric(length = length(x[['cellCommInd']]))
	for (i in 1:length(resVal)) {
		cellVec[x[['cellCommInd']] == i] <- resVal[i]
	}
	if (inherits(x[[1]], 'sf')) {
		x[[1]][metric] <- cellVec
	} else {
		tmp <- terra::rast(x[[1]][[1]])
		tmp[] <- cellVec
		names(tmp) <- metric
		x[[1]] <- c(x[[1]], tmp)
	}
	attributes(x)$metric <- metric
	
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
	
	if (!is.null(tree)) {
		x[['phylo']] <- tree
	}
		
	# update geog species
	x[['geogSpecies']] <- sort(unique(unlist(x[['geogSpecies']])))
	
	return(x)	
}




phyloTips <- function(spNodes, node) {
	names(which(sapply(spNodes, function(y) node %in% y) == TRUE))
}



