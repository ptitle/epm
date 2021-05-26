##' @title Grid Metrics
##' 
##' @description Calculate various morphological and phylogenetic community
##' metrics for every cell in a \code{epmGrid} object. 
##' 
##' @param x object of class \code{epmGrid}
##' 
##' @param metric name of metric to use, see Details. 
##'
##' @param var If a univariate morphological metric is specified, and the 
##' 	data in \code{x} are multivariate, which trait should be used?
##' 	This can also specify which subset of columns a multivariate metric should be applied to.
##' 
##' @param nreps Number of repetitions for Foote metric distribution.
##'
##' @param verbose Intended primarily for debugging, prints progress to the console
##' 
##' @return object of class \code{epmGrid} where the grid represents
##' 	calculations of the metric at every cell. 
##' 
##' @details 
##' 	Univariate trait metrics
##' 	\itemize{
##' 		\item{mean}
##' 		\item{median}
##' 		\item{range}
##'			\item{mean_NN_dist:} {mean nearest neighbor distance}
##'			\item{min_NN_dist:} {minimum nearest neighbor distance}
##' 		\item{variance}
##' 		\item{arithmeticWeightedMean} (see below)
##' 		\item{geometricWeightedMean} (see below)
##'			\item{phylosignal:} {Blomberg's K for phylogenetic signal}
##' 	}
##' 	Multivariate trait metrics
##'		\itemize{
##'			\item{disparity} 
##'			\item{partialDisparity:} {contribution of species in each grid cell to overall disparity}
##' 		\item{range}
##' 		\item{rangePCA}
##' 		\item{mean_NN_dist:} {mean nearest neighbor distance}
##'			\item{min_NN_dist:} {minimum nearest neighbor distance}
##' 		\item{phylosignal:} {Blomberg's K for phylogenetic signal, as implemented in geomorph::physignal}
##' 	}
##' 	Phylogenetic metrics
##' 	\itemize{
##'			\item{pd:} {Faith's phylogenetic diversity, including the root}
##'			\item{meanPatristic}
##'			\item{meanPatristicNN:} {mean nearest neighbor in patristic distance}
##'			\item{minPatristicNN:} {minimum nearest neighbor in patristic distance}
##'			\item{phyloDisparity:} {sum of squared deviations in patristic distance}
##'			\item{PSV:} {Phylogenetic Species Variability}
##'			\item{DR:} {non-parametric estimate of speciation rates}
##' 	}
##' 	Range-weighted metrics
##' 	\itemize{
##'			\item{weightedEndemism}
##'			\item{correctedWeightedEndemism:} {Weighted endemism standardized by species richness}
##'			\item{phyloWeightedEndemism:}
##' 	}
##'
##'		If data slot contains a pairwise matrix, \code{var} is ignored.
##'		Weighted mean options are available where, for each cell, a weighting scheme 
##' 	(inverse of species range sizes) is applied such that small-ranged species are 
##' 	up-weighted, and broadly distributed species are down-weighted. 
##' 	This can be a useful way to lessen the influence of broadly distributed species 
##' 	in the geographic mapping of trait data. 
##'
##' @examples
##' tamiasEPM <- addPhylo(tamiasEPM, tamiasTree)
##' tamiasEPM <- addTraits(tamiasEPM, tamiasTraits)
##'
##' # univariate morphological example
##' x <- gridMetrics(tamiasEPM, metric='mean', var='V2')
##' plot(x) 
##'
##' # multivariate morphological
##' x <- gridMetrics(tamiasEPM, metric='disparity')
##' plot(x)
##' 
##' # phylogenetic metrics
##' x <- gridMetrics(tamiasEPM, metric='meanPatristic')
##' plot(x)
##' 
##' @export


gridMetrics <- function(x, metric, var = NULL, nreps = 20, verbose = FALSE) {
	
	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}
	
	if (!identical(names(x), c('grid', 'speciesList', 'cellCommInd', 'geogSpecies', 'cellCount', 'data', 'phylo'))) {
		stop('enmGrid object is missing some components.')
	}

	if (length(metric) > 1) {
		stop('You can only specify one metric.')
	}
	
	metric <- match.arg(metric, choices = c('mean', 'median', 'range', 'variance', 'arithmeticWeightedMean', 'geometricWeightedMean', 'rangePCA', 'disparity', 'mean_NN_dist', 'min_NN_dist', 'pd', 'meanPatristic', 'meanPatristicNN', 'minPatristicNN', 'phyloDisparity', 'PSV', 'DR', 'weightedEndemism', 'correctedWeightedEndemism', 'phyloWeightedEndemism', 'phylosignal', 'partialDisparity'))
	
	pairwise <- FALSE
	
	# if data is pairwise matrix, then set flag appropriately
	if (inherits(x[['data']], c('matrix', 'data.frame'))) {
		if (identical(rownames(x[['data']]), colnames(x[['data']]))) {
			if (verbose) message('\t...detected pairwise distance matrix...\n') 
			var <- NULL
			pairwise <- TRUE
			# make the diagonal and lower triangle NA
			x[['data']][lower.tri(x[['data']], diag = TRUE)] <- NA
			if (all(is.na(x[['data']][upper.tri(x[['data']])]))) {
				stop('There are no values in the upper triangle of the pairwise matrix.')
			}
		}
	}
	
	# if var is defined but data are univariate, disable var. 
	if (is.vector(x[['data']]) & !is.null(var)) {
		var <- NULL
	}
		
	if (!is.null(var) & inherits(x[['data']], c('matrix', 'data.frame'))) {
		if (is.character(var)) {
			if (any(!var %in% colnames(x[['data']]))) {
				stop('var not a valid column name of the data.')
			}
		} else if (is.numeric(var)) {
			if (any(!all(var %in% 1:ncol(x[['data']])))) {
				stop('Requested data column indices are out of range.')
			}
		}
	}
	
	if (inherits(x[['data']], c('matrix', 'data.frame')) | is.vector(x[['data']])) {
		if (inherits(x[['data']], c('matrix', 'data.frame'))) {
			metricType <- 'multiVar'
			if (!is.null(var) & length(var) == 1) {
				metricType <- 'uniVar'
			}
		} else {
			metricType <- 'uniVar'
		}
	} else {
		metricType <- 'none'
	}
	
	# if a subset of data columns are requested, subset the data table
	if (!is.null(var) & inherits(x[['data']], c('matrix', 'data.frame'))) {
		if (length(var) > 1) {
			x[['data']] <- x[['data']][, var]
		} else {
			x[['data']] <- setNames(x[['data']][, var], rownames(x[['data']]))
		}
		var <- NULL
	}
	
	if (metric %in% c('mean', 'median', 'variance', 'arithmeticWeightedMean', 'geometricWeightedMean') & metricType == 'multiVar' & !pairwise) {
		stop('If a univariate metric is requested from a multivariate dataset, a column name must be provided as var.')
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
		if (is.vector(x[['data']])) {
			sharedTaxa <- intersect(names(x[['data']]), x[['phylo']]$tip.label)
			x[['data']] <- x[['data']][sharedTaxa]
		} else {
			sharedTaxa <- intersect(rownames(x[['data']]), x[['phylo']]$tip.label)
			x[['data']] <- x[['data']][sharedTaxa,]
		}
		x[['phylo']] <- ape::keep.tip(x[['phylo']], sharedTaxa)
	}
	
	# Prune species list according to metric of interest
	if (metric %in% c('mean', 'median', 'disparity', 'partialDisparity', 'range', 'variance', 'arithmeticWeightedMean', 'geometricWeightedMean', 'rangePCA', 'mean_NN_dist', 'min_NN_dist', 'phylosignal')) {
		
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
		
	 } else if (metric %in% c('pd', 'meanPatristic', 'meanPatristicNN', 'minPatristicNN', 'phyloDisparity', 'phyloWeightedEndemism', 'PSV', 'DR')) {
	 	
	 	# check that there is a phylogeny in epmGrid object
		if (is.null(x[['phylo']])) {
			stop('epmGrid object does not contain a phylo object!')
		}
		
		if (verbose) message('\t...dropping species that are not in phylo data...\n')
	 	# prune epmGrid object down to species shared with phylogeny
		x[['speciesList']] <- intersectList(x[['speciesList']], x[['phylo']]$tip.label)
	
	} else if (!metric %in% c('weightedEndemism', 'correctedWeightedEndemism')) {
		stop('Metric not recognized!')
	} else {
		# set empty entries to NA
		emptyEntries <- which(lengths(x[[2]]) == 0)
		x[['speciesList']][emptyEntries] <- rep(list(NA), length(emptyEntries))
	}
	
	uniqueComm <- x[['speciesList']]

	## ----------------------------------
	## MORPHOLOGY-RELATED METRICS
	
	## UNIVARIATE
	
	if (metric %in% c('mean', 'median', 'variance', 'range', 'mean_NN_dist', 'min_NN_dist') & !pairwise & metricType == 'uniVar') {
		if (verbose) message('\t...calculating univariate metric: ', metric, '...\n')
		trait <- x[['data']]
		resVal <- cellAvg(uniqueComm, trait = trait, stat = metric)
		
		# to be consistent with other metrics, we will have 'variance' return 0 for single-species cells
		if (metric == 'variance') {
			ind <- which(lengths(uniqueComm) == 1 & sapply(uniqueComm, anyNA) == FALSE)
			resVal[ind] <- 0
		}
	}
	
	if (metric %in% c('arithmeticWeightedMean', 'geometricWeightedMean') & metricType == 'uniVar') {
		if (verbose) message('\t...calculating univariate metric: ', metric, '...\n')
		trait <- x[['data']]
			
		# get inverse range sizes (1 / cell counts)
		inverseCellCount <- 1 / x[['cellCount']]
		
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(lengths(uniqueComm) > 0 & sapply(uniqueComm, anyNA) == FALSE)

		if (metric == 'arithmeticWeightedMean') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
				(1 / sum(inverseCellCount[y])) * sum(inverseCellCount[y] * trait[y])
			})
		}
		if (metric == 'geometricWeightedMean') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
				exp(sum(inverseCellCount[y] * log(trait[y])) / sum(inverseCellCount[y]))
			})			
		}	
	}
	
	if (metric %in% c('mean', 'median', 'variance') & pairwise) {
	
		# if pairwise matrix
		if (verbose) message('\t...calculating pairwise metric: ', metric, '...\n')
		
		resVal <- rep(NA, length(uniqueComm)) # set up with NA
		
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, anyNA) == FALSE)

		if (metric == 'mean') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) mean(unlist(x[['data']][y, y]), na.rm = TRUE))
		} else if (metric == 'median') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) stats::median(unlist(x[['data']][y, y]), na.rm = TRUE))
		} else if (metric == 'variance') {
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) stats::var(unlist(x[['data']][y, y]), na.rm = TRUE))
			ind <- which(lengths(uniqueComm) == 1 & sapply(uniqueComm, anyNA) == FALSE)
			resVal[ind] <- 0
		}
	}
	
	if (metric == 'phylosignal' & metricType == 'uniVar') {
		if (verbose & metricType == 'uniVar') message('\t...calculating univariate metric: ', metric, '...\n')
		resVal <- rep(NA, length(uniqueComm)) # set up with NA
		ind <- which(lengths(uniqueComm) > 2)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) geomorph::physignal(as.matrix(x[['data']][y]), phy = ape::keep.tip(x$phylo, y), iter = 999, print.progress = FALSE)$phy.signal)	
	}

	## MULTIVARIATE
	
	if (metric == 'disparity' & metricType == 'multiVar') {
		if (verbose) message('\t...calculating multivariate metric: ', metric, '...\n')
		# sum of the diagonal of the covariance matrix
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		# resVal <- rep(NA, length(uniqueComm)) # set up with NA
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, length) > 1)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) sum(diag(cov(x[['data']][y,]))))
	}
	
	if (metric == 'partialDisparity' & metricType == 'multiVar') {
		if (verbose) message('\t...calculating multivariate metric: ', metric, '...\n')
		partialDisp <- getSpPartialDisparities(x[['data']])
		# sum of the diagonal of the covariance matrix
		resVal <- pbapply::pbsapply(uniqueComm, function(y) sum(partialDisp[y]))
	}	
	
	if (metric == 'range' & metricType == 'multiVar') {
		if (verbose) message('\t...calculating multivariate metric: ', metric, '...\n')
		# maximum of the distance matrix (0 if one sp)
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, length) > 1)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) max(dist(x[['data']][y, ])))
	}
	
	if (metric == 'rangePCA' & metricType == 'multiVar') {
		if (verbose) message('\t...calculating multivariate metric: ', metric, '...\n')
		pc <- prcomp(x[['data']])
		# retain 99% of the variation
		keep <- 1:which(cumsum(((pc$sdev)^2) / sum(pc$sdev^2)) >= 0.99)[1]
		if (length(keep) == 1) { keep <- 1:ncol(pc$x) }
		pc <- pc$x[,keep]
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, length) > 1)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
			sum(apply(pc[y,], 2, function(z) diff(range(z))))
		})	
	}
	
	if (metric == 'mean_NN_dist' & metricType == 'multiVar') {
		if (verbose) message('\t...calculating multivariate metric: ', metric, '...\n')
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, length) > 1)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) mean(nnDist(x[['data']][y, ], Nrep = nreps)$mean_dist))
	}

	if (metric == 'min_NN_dist' & metricType == 'multiVar') {
		if (verbose) message('\t...calculating multivariate metric: ', metric, '...\n')
		resVal <- numeric(length = length(uniqueComm)) # set up with zeros
		resVal[sapply(uniqueComm, anyNA)] <- NA
		ind <- which(sapply(uniqueComm, length) > 1)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) min(nnDist(x[['data']][y, ], Nrep = nreps)$mean_dist))
	}
	
	if (metric == 'phylosignal' & metricType == 'multiVar') {
		if (verbose & metricType == 'multiVar') message('\t...calculating multivariate metric: ', metric, '...\n')
		resVal <- rep(NA, length(uniqueComm)) # set up with NA
		ind <- which(lengths(uniqueComm) > 2)
		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) geomorph::physignal(as.matrix(x[['data']][y, ]), phy = ape::keep.tip(x$phylo, y), iter = 999, print.progress = FALSE)$phy.signal)	
	}
	
	## ----------------------------------
	## PHYLOGENY-RELATED METRICS
	
	if (metric %in% c('pd', 'meanPatristic', 'minPatristicNN', 'meanPatristicNN', 'phyloDisparity', 'PSV', 'DR')) {
		
		# calculate pairwise patristic distance
		patdist <- cophenetic(x[['phylo']])
		diag(patdist) <- NA
		
		if (metric == 'pd') {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			mat <- matrix(0, nrow = length(uniqueComm), ncol = length(x$phylo$tip.label))
			colnames(mat) <- x$phylo$tip.label
			for (i in 1:length(uniqueComm)) {
				mat[i, colnames(mat) %in% uniqueComm[[i]]] <- 1
			}
			resVal <- picante::pd(mat, x$phylo, include.root = TRUE)$PD
			resVal[sapply(uniqueComm, anyNA)] <- NA
		}
		
		if (metric == 'meanPatristic') {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			# meanPatristic is 0 if 1 species, NA if no species
			patdist[upper.tri(patdist, diag = TRUE)] <- NA
			resVal <- pbapply::pbsapply(uniqueComm, function(y) mean(unlist(patdist[y, y]), na.rm = TRUE))
			resVal[sapply(uniqueComm, length) == 1] <- 0
			resVal[sapply(uniqueComm, anyNA)] <- NA
		}
		
		if (metric == 'meanPatristicNN') {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			# the mean of the minimum patristic distance for each species present
			resVal <- numeric(length = length(uniqueComm)) # set up with zeros
			resVal[sapply(uniqueComm, anyNA)] <- NA
			ind <- which(sapply(uniqueComm, length) > 1)
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
				return(mean(apply(patdist[y, y], MARGIN = 1, min, na.rm = TRUE)))
			})
		}

		if (metric == 'minPatristicNN') {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			# the minimum of the minimum patristic distance for each species present
			resVal <- numeric(length = length(uniqueComm)) # set up with zeros
			resVal[sapply(uniqueComm, anyNA)] <- NA
			ind <- which(sapply(uniqueComm, length) > 1)
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
				return(min(apply(patdist[y, y], MARGIN = 1, min, na.rm = TRUE)))
			})
		}

		if (metric == 'phyloDisparity') {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			# the sum of the squared deviations from the mean
			# value of 0 if 1 species, NA if no species
			
			resVal <- numeric(length = length(uniqueComm)) # set up with zeros
			resVal[sapply(uniqueComm, anyNA)] <- NA
			ind <- which(sapply(uniqueComm, length) > 1)

			patdist[upper.tri(patdist, diag = TRUE)] <- NA
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) sum((patdist[y, y] - mean(patdist[y,y], na.rm = TRUE)) ^ 2, na.rm = TRUE) / choose(length(y), 2))
		}		
		
		if (metric == 'PSV') {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			# Phylogenetic Species Variability is based on the variance-covariance matrix
			# measure of phylogenetic diversity, ranges from 0-1, not confounded by species richness
			
			vmat <- ape::vcv.phylo(x[['phylo']], corr = TRUE)
			
			resVal <- numeric(length = length(uniqueComm)) # set up with zeros
			resVal[sapply(uniqueComm, anyNA)] <- NA
			ind <- which(sapply(uniqueComm, length) > 1)
			
			resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) (length(y) * sum(diag(vmat[y,y])) - sum(vmat[y,y]))/(length(y) * (length(y) - 1)))
		}
		
		if (metric == 'DR') {
			if (verbose) message('\t...calculating phylo metric: ', metric, '...\n')
			# DR statistic for speciation rates
			
			tipDR <- DRstat(x[['phylo']])			
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
		spEdges <- getRootToTipEdges(x[['phylo']])
		if (!'edgeArea' %in% names(x)) {
			if (verbose) message('\t...calculating branch-specific range sizes...\n')

			# phylo branch ranges must be based on full list of cell communities, therefore we need to expand the speciesList
			fullSpList <- expandSpeciesCellList(x)
			x[['edgeArea']] <- do.call(cbind, phyloBranchRanges(x[['phylo']], convertNAtoEmpty(fullSpList), spEdges))
		}
		tipIndVec <- sapply(x[['phylo']]$tip.label, function(z) which(x[['phylo']]$tip.label == z))
		resVal <- rep(NA, length(uniqueComm)) # set up with NA
		ind <- which(sapply(uniqueComm, anyNA) == FALSE)

		resVal[ind] <- pbapply::pbsapply(uniqueComm[ind], function(y) {
			commEdges <- unique(unlist(spEdges[tipIndVec[y]])) + 1
			sub <- x[['edgeArea']][commEdges,]
			if (inherits(sub, 'numeric')) {
				sub <- matrix(sub, nrow = 1)
			}
			sum(sub[,1] / sub[,2])
		})
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
	## Maybe we don't want to do this, as we want the species richness to be our best understanding.
	# if (inherits(x[[1]], 'sf')) {
		# for (i in 1:length(x[['speciesList']])) {
			# x[['grid']][which(x[['cellCommInd']] == i), 'spRichness'] <- length(x[['speciesList']][[i]])
		# }
	# } else {
		# tmpVec <- terra::values(x[['grid']][['spRichness']])
		# for (i in 1:length(x[['speciesList']])) {
			# tmpVec[which(x[['cellCommInd']] == i)] <- length(x[['speciesList']][[i]])
		# }
		# terra::values(x[['grid']][['spRichness']]) <- tmpVec
		# rm(tmpVec)
	# }
		
	# update geog species
	x[['geogSpecies']] <- sort(unique(unlist(x[['geogSpecies']])))
	
	return(x)	
}




# # function to calculate phylo branch specific range areas

# phyloBranchRanges <- function(x) {

	# if (!'speciesRaster' %in% class(x)) {
		# stop('x must be of class speciesRaster.')
	# }
	
	# # check that there is a phylogeny in speciesRaster object
	# if (is.null(x[['phylo']])) {
		# stop('speciesRaster object does not contain a phylo object!')
	# }
	
	# phylo <- x[['phylo']]
	
	# allLeaves <- phangorn::Descendants(phylo, phylo$edge[,2], type = 'tips')
	# allLeaves <- lapply(allLeaves, function(y) phylo$tip.label[y])
	
	# branchTable <- matrix(nrow = length(phylo$edge.length), ncol = 2)
	# colnames(branchTable) <- c('branchLength', 'branchArea')
	# branchTable[,1] <- phylo$edge.length
	
	# for (i in 1:length(phylo$edge.length)) {
		
		# inCell <- sapply(x[['speciesList']], function(y) any(allLeaves[[i]] %in% y))
		# branchTable[i, 2] <- sum(inCell)	
	# }

	# x[['edgeArea']] <- branchTable
	# return(x)	
# }



# getCommEdges2 <- function(phylo, comm, indList) {
	
	# tipInd <- sapply(comm, function(z) which(phylo$tip.label == z))
	# edges <- unique(unlist(indList[tipInd]))
	# edges <- edges + 1
	# # edgeColors <- ifelse(1:length(phylo$edge.length) %in% edges, 'red','black')
	# # edgeWidths <- ifelse(1:length(phylo$edge.length) %in% edges, 2, 1)
	# # plot(phylo, tip.color=ifelse(phylo$tip.label %in% comm, 'blue','gray'), root.edge=TRUE, edge.color = edgeColors, edge.width = edgeWidths)
	# # nodelabels(node=nodes, frame='circle', bg='red', cex=0.1)
	
	# return(edges)
# }




# # library(raster)
# library(maptools)
# library(Rcpp)
# # convert polygon ranges to raster
# ranges <- rasterStackFromPolyList(tamiasPolyList, resolution = 20000)
# x <- createSpeciesRaster(ranges = ranges)
# x <- addPhylo_speciesRaster(x, tamiasTree)


convertNAtoEmpty <- function(spCellList) {
	ind <- which(sapply(spCellList, anyNA) == TRUE)
	spCellList[ind] <- rep(list('empty'), length(ind))
	return(spCellList)
}



# # phylo2 <- x[['phylo']]
# class(phylo2) <- 'list'

# sourceCpp('~/Desktop/testing.cpp')

# tipEdges <- getRootToTipEdges(phylo2)

# q <- phyloBranchRanges(phylo2, spCellList, tipEdges)
# q <- do.call(cbind, q)



