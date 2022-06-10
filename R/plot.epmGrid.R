##' @title Plot epmGrid
##'
##' @description Plot a epmGrid object. This function uses the tmap package 
##' for plotting by default.
##'
##' @param x object of class \code{epmGrid}
##' @param log boolean; should the cell values be logged?
##' @param colorRampRange numeric vector of min and max value for scaling the 
##' color ramp. Automatically inferred if set to \code{NULL}. 
##' This is relevant if multiple plots are desired on the same scale. 
##' See \code{\link{getMultiMapRamp}}. 
##' @param legend boolean; should legend be included?
##' @param col either a vector of color names that will be interpolated, 
##' or a color ramp function that takes an integer 
##' (see for example \code{\link{colorRampPalette}}).
##'	@param basemap if \code{'none'}, then only the grid is plotted. 
##'		If \code{'worldmap'}, then vector map is plotted.
##' 	If \code{'interactive'}, then the plot is sent to the web browser.
##' @param minTaxCount an integer, or 'auto'. Should cells 
##' containing certain numbers of taxa be grayed out? For example, 
##' should single-taxon cells be ignored because the metric only makes sense 
##' for multi-taxon cells? This is predetermined for all metrics in
##'  \code{\link{gridMetrics}} if  minTaxCount = 'auto'.
##' @param ignoredColor color for ignored cells. See details.
##' @param lwd grid cell border width
##' @param borderCol color for grid cell borders
##' @param alpha opacity of all colors and borders, ranging from 0 
##' (fully transparent) to 1 (fully opaque)
##' @param includeFrame boolean; include frame around plot?
##' @param use_tmap boolean; if FALSE, plotting will be done via sf 
##' instead of tmap package
##' @param fastPoints Intended for debugging purposes. For hex grids and \code{use_tmap = F}, plot points instead of polygons. Helpful for sorting out plotting 
##' details without waiting for slow polygon plotting.
##'	@param title text to add to the plot
##' @param add logical, add to existing plot?
##' @param ... additional arguments that can be passed to sf::plot or terra::plot 
##' 		if \code{use_tmap = FALSE}
##'
##'
##' @details If \code{x} is a metric as generated with \code{gridMetrics} 
##' that returns 0 for single-species cells, then those cells 
##' (that have a value of 0) will be plotted in gray (or any color as specified 
##' with \code{ignoredColor}) if \code{minTaxCount = 'auto'}. You can specify
##' other values as well. For instance, if you use the function
##' \code{\link{customGridMetric}} to calculate phylogenetic signal, which is 
##' a metric that only makes sense for cells with 3 or more taxa, then you could 
##' then specify \code{minTaxCount = 3}. Setting \code{minTaxCount = 1} shows all
##' cells with data.
##'
##'	If the tmap package is not installed, then this function will default 
##'	to plotting with \code{sf::plot}.
##' 
##'	If you would like more control over the legend, then plot with 
##'	\code{tmap = FALSE} and \code{legend = FALSE}, and then call the function
##'	\code{\link{addLegend}}. 
##'
##' @return Nothing is returned if plotting with tmap (the default). 
##' If plotting with \code{use_tmap = FALSE}, and if the plot is directed to 
##' a variable, then this variable will contain relevant information to be passed
##' on to the function \code{\link{addLegend}}: 
##'
##' @author Pascal Title
##' 
##' @examples
##' plot(tamiasEPM)
##' 
##' plot(tamiasEPM, legend = FALSE, use_tmap = FALSE, col = viridisLite::inferno)
##' addLegend(tamiasEPM, location = 'top', ramp = viridisLite::inferno)
##' 
##' # Example for how to plot multiple epmGrids on the same color scale
##' # for illustration purposes, we will compare weighted endemism to
##' # phylogenetic weighted endemism
##' library(tmap)
##'
##' tamiasEPM <- addPhylo(tamiasEPM, tamiasTree)
##' epm1 <- gridMetrics(tamiasEPM, metric='weightedEndemism')
##' epm2 <- gridMetrics(tamiasEPM, metric='phyloWeightedEndemism')
##' # get global min and max values
##' minmax <- getMultiMapRamp(epm1, epm2)
##' 
##' map1 <- plot(epm1, colorRampRange = log(minmax), log = TRUE, legend = FALSE)
##' map2 <- plot(epm2, colorRampRange = log(minmax), log = TRUE, legend = FALSE)
##' # tmap_arrange(map1, map2)
##' \donttest{
##' # view your plot in the web-browser as a dynamic plot.
##' plot(tamiasEPM, basemap = 'interactive')
##' }
##'
##' # Adding a custom legend, and passing along arguments via params
##' xx <- plot(tamiasEPM, use_tmap = FALSE, legend = FALSE, 
##' col = viridisLite::magma)
##' addLegend(tamiasEPM, params = xx, location = 'bottom')
##'
##' @rdname plot
##' @aliases plot.epmGrid
##' @export
##' @export plot.epmGrid

plot.epmGrid <- function(x, log = FALSE, legend = TRUE, col, basemap = 'worldmap', colorRampRange = NULL, minTaxCount = 'auto', ignoredColor = gray(0.9), lwd, borderCol = 'black', alpha = 1, includeFrame = FALSE, use_tmap = TRUE, fastPoints = FALSE, title = '', add = FALSE, ...) {
	
	# x = tamiasEPM; log = FALSE; legend = TRUE; basemap = 'worldmap'; colorRampRange = NULL; ignoredColor = gray(0.9); lwd = 0.25; borderCol = 'black'; includeFrame = FALSE; use_tmap = TRUE; alpha = 1; add = FALSE; fastPoints = FALSE; minTaxCount = 'auto'; title = ''
	
	if (!inherits(x, 'epmGrid')) {
		stop('Object must be of class epmGrid')
	}
	
	if (!basemap %in% c('worldmap', 'interactive', 'none')) {
		stop("basemap argument must be 'none', 'worldmap' or 'interactive'.")
	}
	
	if (use_tmap & !requireNamespace('tmap', quietly = TRUE)) {
		message('\ttmap package not installed -- defaulting to sf/terra plot')
		use_tmap <- FALSE
	}
	
	if (!is.numeric(minTaxCount) & !all(minTaxCount == 'auto')) {
		stop('minTaxCount must be either an integer or "auto".')
	}
	if (is.numeric(minTaxCount)) {
		if (minTaxCount == 0) {
			stop('minTaxCount should be 1 or greater.')
		}
	}
		
	if (fastPoints) {
		use_tmap <- FALSE
	}
	
	if (add) {
		use_tmap <- FALSE
		basemap <- 'none'
		legend <- FALSE
	}
	
	plotMetric <- attributes(x)$metric
	
	# if x is a epmGrid that represents a metric that only makes sense for communities with multiple species, 
	# then single species cells have a value of zero, and we will plot those cells as gray.
	# We must also then specify a new minimum bound for the color palette to be the min of multi-sp cells.

	if (minTaxCount == 'auto') {
		if (plotMetric %in% c('range', 'mean_NN_dist', 'min_NN_dist', 'evenness', 'variance', 'disparity', 'rangePCA', 'meanPatristic', 'meanPatristicNN', 'minPatristicNN', 'phyloEvenness', 'phyloDisparity', 'PSV', 'PSR')) {
			minTaxCount <- 2
		} else {
			minTaxCount <- 1
		}
	}

	if (minTaxCount > 1) {
		# determine which cells have too few species
		## we do minTaxCount - 1 because we want to identify the cells to exclude.
		tooFewInd <- spCountIndex(x, 1:(minTaxCount - 1))
		
		# check
		# table(lengths(x[[2]][x[['cellCommInd']][tooFewInd]]))
		
		if (inherits(x[[1]], 'sf')) {
			# add on cells that do not have species
			# singleSpCells <- sort(c(singleSpCells, which(x[['cellCommInd']] %in% which(sapply(x[['speciesList']], anyNA) == TRUE))))
			grid_multiSp <- x[[1]][- tooFewInd,]
			grid_singleSp <- x[[1]][tooFewInd,]
		} else {
			grid_multiSp <- x[[1]][plotMetric]
			grid_singleSp <- terra::rast(x[[1]][plotMetric])
			grid_multiSp[tooFewInd] <- NA
			grid_singleSp[] <- NA
			grid_singleSp[tooFewInd] <- 1
		}	
	} 	
	
	ncolors <- 100
	isInt <- FALSE
	if (inherits(x[[1]], 'sf')) {
		if (all(integercheck(sf::st_drop_geometry(x[[1]])[, plotMetric]))) {
			isInt <- TRUE
			if (max(x[[1]][[plotMetric]]) <= 10) {
				ncolors <- max(x[[1]][[plotMetric]])
			}
		}	
	} else if (inherits(x[[1]], 'SpatRaster')) {
		samp <- sample(as.vector(stats::na.omit(terra::values(x[[1]][plotMetric]))), size = 1000, replace = TRUE)
		if (all(integercheck(samp))) {
			isInt <- TRUE
			if (max(terra::minmax(x[[1]][plotMetric])) <= 10) {
				ncolors <- max(terra::minmax(x[[1]][plotMetric]))
			}
		}
	}

	if (missing(col)) {
		# colramp <- sf::sf.colors
		# default color ramp: viridis turbo, but without the most extreme values
		colramp <- function(n) viridisLite::turbo(n = n, begin = 0.1, end = 0.9)
		# colramp <- grDevices::colorRampPalette(c('blue', 'cyan', 'yellow', 'red'))
	} else {
		if (inherits(col, 'function')) {
			colramp <- col
		} else if (inherits(col, 'character') & length(col) > 1) {
				colramp <- grDevices::colorRampPalette(col)
		} else {
			stop('color not understood.')
		}
	}

	if (missing(lwd)) {
		lwd <- 0.15
	}
	
	metricName <- plotMetric
	
	####################################
	## HEXAGONAL GRIDS
	
	if (inherits(x[[1]], 'sf')) {
	
		if (log) {
			x[[1]]$loggedMetric <- log(x[[1]][[plotMetric]])
			metricName <- paste0('log ', plotMetric)
			isInt <- FALSE
			plotMetric <- 'loggedMetric'
		}
		
		if (is.null(colorRampRange)) {
			breaks <- NULL
		} else {
			breaks <- seq(min(colorRampRange), max(colorRampRange), length.out = 5)
		}
					
		if (use_tmap) {
		    
		    if (basemap != 'interactive' & getOption("tmap.mode") == 'view') {
		        tmap::tmap_mode('plot')
		    }
		
			if (basemap == 'worldmap') {

				# wrld <- sf::st_transform(worldmap, crs = sf::st_crs(x[[1]]))
				map <- tmap::tm_shape(worldmap, is.master = FALSE, projection = sf::st_crs(x[[1]]), bbox = x[[1]]) + tmap::tm_borders(lwd = 0.25)
			}
			
			if (basemap == 'none') {
				map <- NULL
			}
			
			if (basemap == 'interactive') {
				map <- NULL
				tmap::tmap_mode('view')
			}
			
			if (isInt & ncolors <= 10) {
				tmapStyle <- 'pretty'
			} else {
				tmapStyle <- 'cont'
			}
			
			if (minTaxCount <= 1) {
				
				map <- map + tmap::tm_shape(x[[1]]) + tmap::tm_fill(plotMetric, palette = colramp(ncolors), legend.show = legend, title = title, breaks = breaks, style = tmapStyle, midpoint = NA, alpha = alpha, legend.reverse = T) + tmap::tm_borders(col = borderCol, lwd = lwd, alpha = alpha) + tmap::tm_layout(frame = includeFrame, legend.outside = TRUE)
								
			} else {
				
				map <- map + tmap::tm_shape(grid_multiSp, is.master = TRUE, bbox = x[[1]]) + tmap::tm_fill(plotMetric, palette = colramp(ncolors), legend.show = legend, title = title, breaks = breaks, style = tmapStyle, midpoint = NA, alpha = alpha, legend.reverse = T) + tmap::tm_borders(col = borderCol, lwd = lwd, alpha = alpha) + tmap::tm_shape(grid_singleSp) + tmap::tm_fill(plotMetric, palette = ignoredColor, legend.show = FALSE, alpha = alpha) + tmap::tm_borders(col = borderCol, lwd = lwd, alpha = alpha) + tmap::tm_layout(frame = includeFrame, legend.outside = TRUE)
			
			}
		
			map
	
		} else {
			# use sf plotting
			
			# if (legend) {
				# key.pos <- 4
			# } else {
				# key.pos <- NULL
			# }
			
			if (is.null(colorRampRange)) {
				if (minTaxCount <= 1) {
					valRange <- range(x[[1]][[plotMetric]], na.rm = TRUE)
				} else {
					valRange <- range(grid_multiSp[[plotMetric]], na.rm = TRUE)
				}
				breaks <- seq(min(valRange), max(valRange), length.out = ncolors + 1)
			} else {
				breaks <- seq(min(colorRampRange), max(colorRampRange), length.out = ncolors + 1)
			}
			
			colors <- colramp(ncolors)
			colors <- grDevices::adjustcolor(colors, alpha.f = alpha)
			borderColor <- grDevices::adjustcolor('black', alpha.f = alpha)
			ignoredColor <- grDevices::adjustcolor(ignoredColor, alpha.f = alpha)
			
			# helpful for debugging because much faster than plotting polygons
			if (fastPoints) {
				if (minTaxCount <= 1) {
					cols <- colors[findInterval(x[[1]][[plotMetric]], breaks, rightmost.closed = TRUE)]
					plot(sf::st_centroid(sf::st_geometry(x[[1]][plotMetric])), pch = 20, reset = FALSE, col = cols, cex = 0.5, add = add)
				} else {
					cols <- colors[findInterval(grid_multiSp[[plotMetric]], breaks, rightmost.closed = TRUE)]
					plot(sf::st_centroid(sf::st_geometry(grid_singleSp[plotMetric])), pch = 20, reset = FALSE, col = ignoredColor, cex = 0.5, add = add)
					plot(sf::st_centroid(sf::st_geometry(grid_multiSp[plotMetric])), pch = 20, reset = FALSE, col = cols, cex = 0.5, add = TRUE)	
				}
			} else {
	

				if (minTaxCount <= 1) {
					cols <- colors[findInterval(x[[1]][[plotMetric]], breaks, rightmost.closed = TRUE)]
					plot(sf::st_geometry(x[[1]][plotMetric]), axes = includeFrame, main = NULL, key.pos = NULL, lwd = lwd, reset = FALSE, col = cols, border = borderColor, add = add, ...)
				} else {
					cols <- colors[findInterval(grid_multiSp[[plotMetric]], breaks, rightmost.closed = TRUE)]
					plot(sf::st_geometry(grid_singleSp[plotMetric]), axes = includeFrame, main = NULL, key.pos = NULL, col = ignoredColor, border = NA, reset = FALSE, add = add, ...)
					plot(sf::st_geometry(grid_multiSp[plotMetric]), lwd = lwd, col = cols, border = borderCol, add = TRUE, reset = FALSE)
				}
			



# # 				if (!plotSingleCells) {
					# plot(x[[1]][plotMetric], axes = includeFrame, main = NULL, key.pos = NULL, lwd = lwd, reset = FALSE, pal = colors, border = borderColor, breaks = breaks, add = add, ...)
				# } else {
					# plot(grid_singleSp[plotMetric], axes = includeFrame, main = NULL, key.pos = NULL, pal = ignoredColor, border = NA, reset = FALSE, add = add, ...)
					# plot(grid_multiSp[plotMetric], lwd = lwd, pal = colors, border = borderCol, breaks = breaks, add = TRUE, reset = FALSE)
				# }
			}
			
			if (legend) {
	
				if (minTaxCount <= 1) {
					if (isInt & ncolors <= 10) {
						if (inherits(x[[1]], 'sf')) {
							nTicks <- max(x[[1]][[plotMetric]])
						} else if (inherits(x[[1]], 'SpatRaster')) {
							nTicks <- max(terra::minmax(x[[1]][plotMetric]))
						}
					} else {
						nTicks <- 3
					}
					
					minmax <- colorRampRange

					addLegend(x[[1]][plotMetric], location = 'topright', ramp = colramp, isInteger = isInt, ncolors = ncolors, nTicks = nTicks, minmax = minmax)

				} else {

					if (isInt & ncolors <= 10) {
						if (inherits(x[[1]], 'sf')) {
							nTicks <- max(grid_multiSp[plotMetric])
						} else if (inherits(x[[1]], 'SpatRaster')) {
							nTicks <- max(terra::minmax(grid_multiSp[plotMetric]))
						}
					} else {
						nTicks <- 3
					}
					
					if (is.null(colorRampRange)) {
						minmax <- range(sf::st_drop_geometry(grid_multiSp[plotMetric]))
					} else {
						minmax <- colorRampRange
					}

					addLegend(grid_multiSp[plotMetric], location = 'topright', ramp = colramp, isInteger = isInt, ncolors = ncolors, nTicks = nTicks, minmax = minmax)
				}
			}
	
			if (basemap == 'worldmap') {
				# add map for context
				wrld <- sf::st_transform(worldmap, crs = sf::st_crs(x[[1]]))
				wrld <- sf::st_cast(wrld, 'MULTILINESTRING')
				grXY <- graphics::par("usr")
				clip <- sf::st_make_grid(sf::st_as_sf(rbind.data.frame(grXY[c(1,3)], grXY[c(2,4)]), coords = 1:2, crs = sf::st_crs(x[[1]])), n = 1)
				wrld <- sf::st_intersection(wrld, clip)
				wrld <- sf::st_combine(wrld)				
				# graphics::clip(grXY[1], grXY[2], grXY[3], grXY[4]) # this ensures that world map is constrained to plot region
				graphics::plot(wrld, add = TRUE, lwd = lwd)
			}
			
			return(invisible(list(
				log = log,
				minmax = valRange,
				minTaxCount = minTaxCount,
				colramp = colramp,
				alpha = alpha,
				isInteger = isInt
			)))
		}

	} else if (inherits(x[[1]], 'SpatRaster')) {
		
			# get bounding box of grid cells with data. We will do this so that we can leave out of the plot those grid cells that are empty and that may take up quite a bit of geographic space. Only an issue for rasters.
			datBB <- sf::st_bbox(sf::st_as_sf(as.data.frame(terra::xyFromCell(x[[1]], which(!is.na(terra::values(x[[1]][[plotMetric]]))))), coords = 1:2, crs = sf::st_crs(x[[1]])))
		
				
		if (log) {
			metricName <- paste0('log ', plotMetric)
			if (minTaxCount <= 1) {
				metricMap <- log(x[[1]][plotMetric])
			} else {
				metricMap <- log(grid_multiSp[plotMetric])
			}
			isInt <- FALSE
		} else {
			if (minTaxCount <= 1) {
				metricMap <- x[[1]][plotMetric]
			} else {
				metricMap <- grid_multiSp[plotMetric]
			}
		}
		
		if (use_tmap) {
			
			if (is.null(colorRampRange)) {
				breaks <- NULL
			} else {
				breaks <- seq(min(colorRampRange), max(colorRampRange), length.out = 5)
			}
			    
		    if (basemap != 'interactive' & getOption("tmap.mode") == 'view') {
		        tmap::tmap_mode('plot')
		    }
					
			if (basemap == 'worldmap') {

				# wrld <- sf::st_transform(worldmap, crs = sf::st_crs(x[[1]]))
				map <- tmap::tm_shape(worldmap, is.master = FALSE, projection = sf::st_crs(x[[1]]), bbox = datBB) + tmap::tm_borders(lwd = 0.25)
			}
			
			if (basemap == 'none') {
				map <- NULL
			}
			
			if (basemap == 'interactive') {
				map <- NULL
				tmap::tmap_mode('view')
			}
			
			if (isInt & ncolors <= 10) {
				tmapStyle <- 'cat'
			} else {
				tmapStyle <- 'cont'
			}
			
			if (minTaxCount <= 1) {
				map <- map + tmap::tm_shape(metricMap, bbox = datBB) + tmap::tm_raster(palette = colramp(ncolors), legend.show = legend, title = title, breaks = breaks, style = tmapStyle, midpoint = NA, alpha = alpha, legend.reverse = T) + tmap::tm_layout(frame = includeFrame, legend.outside = TRUE)
				
			} else {
				
				map <- map + tmap::tm_shape(metricMap, is.master = TRUE, bbox = datBB) + tmap::tm_raster(palette = colramp(ncolors), legend.show = legend, title = title, breaks = breaks, style = tmapStyle, midpoint = NA, alpha = alpha, legend.reverse = T) + tmap::tm_shape(grid_singleSp) + tmap::tm_raster(palette = ignoredColor, legend.show = FALSE, alpha = alpha) + tmap::tm_layout(frame = includeFrame, legend.outside = TRUE)
			
			}
		
			map
	
		} else {

			# use terra plotting
			
			datBB2 <- terra::ext(x[[1]])
			datBB2[1] <- datBB['xmin']
			datBB2[2] <- datBB['xmax']
			datBB2[3] <- datBB['ymin']
			datBB2[4] <- datBB['ymax']
					

			if (is.null(colorRampRange)) {
				valRange <- range(terra::minmax(metricMap))
			} else {
				valRange <- colorRampRange
			}
			

			if (minTaxCount <= 1) {
				terra::plot(terra::crop(metricMap, datBB2), col = colramp(ncolors), axes = includeFrame, legend = FALSE, plg = list(shrink = 0.7, title = metricName), range = valRange, alpha = alpha, add = add, ...)
				
			} else {
				
				terra::plot(terra::crop(metricMap, datBB2), col = colramp(ncolors), axes = includeFrame, legend = FALSE, plg = list(shrink = 0.7, title = metricName), range = valRange, alpha = alpha, add = add, ...)
				terra::plot(grid_singleSp, col = ignoredColor, axes = includeFrame, legend = FALSE, range = valRange, add = TRUE, alpha = alpha)
				
			}
			
			if (legend) {
	
				if (minTaxCount <= 1) {
					if (isInt & ncolors <= 10) {
						if (inherits(x[[1]], 'sf')) {
							nTicks <- max(x[[1]][[plotMetric]])
						} else if (inherits(x[[1]], 'SpatRaster')) {
							nTicks <- max(terra::minmax(x[[1]][plotMetric]))
						}
					} else {
						nTicks <- 3
					}
					
					minmax <- valRange

					addLegend(x[[1]][plotMetric], location = 'topright', ramp = colramp, isInteger = isInt, ncolors = ncolors, nTicks = nTicks, minmax = minmax)

				} else {

					if (isInt & ncolors <= 10) {
						if (inherits(x[[1]], 'sf')) {
							nTicks <- max(grid_multiSp[plotMetric])
						} else if (inherits(x[[1]], 'SpatRaster')) {
							nTicks <- max(terra::minmax(grid_multiSp[plotMetric]))
						}
					} else {
						nTicks <- 3
					}
					
					if (is.null(colorRampRange)) {
						minmax <- range(terra::minmax(grid_multiSp[plotMetric]))
					} else {
						minmax <- colorRampRange
					}

					addLegend(grid_multiSp[plotMetric], location = 'topright', ramp = colramp, isInteger = isInt, ncolors = ncolors, nTicks = nTicks, minmax = minmax)
				}
			}
			
	
			if (basemap == 'worldmap') {
				# add map for context
				wrld <- sf::st_transform(worldmap, crs = sf::st_crs(x[[1]]))
				wrld <- sf::st_cast(wrld, 'MULTILINESTRING')
				grXY <- graphics::par("usr")
				clip <- sf::st_make_grid(sf::st_as_sf(rbind.data.frame(grXY[c(1,3)], grXY[c(2,4)]), coords = 1:2, crs = sf::st_crs(x[[1]])), n = 1)
				wrld <- sf::st_intersection(wrld, clip)
				wrld <- sf::st_combine(wrld)				
				# graphics::clip(grXY[1], grXY[2], grXY[3], grXY[4]) # this ensures that world map is constrained to plot region
				graphics::plot(wrld, add = TRUE, lwd = lwd)
			}
			
			return(invisible(list(
				log = log,
				minmax = valRange,
				minTaxCount = minTaxCount,
				colramp = colramp,
				alpha = alpha,
				isInteger = isInt
			)))

		}
	}	
}
		
integercheck <- function(x) {
	abs(x - round(x)) < .Machine$double.eps ^ 0.5
}


