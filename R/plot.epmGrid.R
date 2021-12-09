##' @title Plot epmGrid
##'
##' @description Plot a epmGrid object. This function uses the tmap package for plotting by default.
##'
##' @param x object of class \code{epmGrid}
##' @param log boolean; should the cell values be logged?
##' @param colorRampRange numeric vector of min and max value for scaling the color
##' 	ramp. Automatically inferred if set to \code{NULL}. This is relevant if multiple
##' 	plots are desired on the same scale. See \code{\link{getMultiMapRamp}}. 
##' 		Not intended for mapview option.
##' @param legend boolean; should legend be included?
##' @param col either a vector of color names that will be interpolated, or a color ramp
##' 	function that takes an integer (see for example \code{\link{colorRampPalette}}).
##'	@param basemap if \code{'none'}, then only the grid is plotted. 
##'		If \code{'worldmap'}, then vector map is plotted.
##' 	If \code{'interactive'}, then the \code{mapview} package is used.
##' @param singleSpCol color for single-species cells. See details.
##' @param lwd grid cell border width
##' @param borderCol color for grid cell borders
##' @param alpha opacity of all colors and borders, ranging from 0 (fully transparent) 
##' 		to 1 (fully opaque)
##' @param includeFrame boolean; include frame around plot?
##' @param use_tmap boolean; if FALSE, plotting will be done via sf instead of tmap package
##' @param fastPoints Intended for debugging purposes. For hex grids and use_tmap = F, 
##' 	plot points instead of polygons.
##' @param add logical, add to existing plot?
##' @param ... additional arguments that can be passed to sf::plot or terra::plot 
##' 		if \code{use_tmap = FALSE}
##'
##'
##' @details If \code{x} is a metric as generated with \code{gridMetrics} that returns 0 
##' 	for single-species cells, then those cells (that have a value of 0) will be plotted in gray (or any color
##' 	as specified with \code{singleSpCol}).
##'
##'		If the tmap package is not installed, then this function will default to plotting with \code{sf::plot}.
##' 
##' @return Nothing is returned. 
##'
##' @author Pascal Title
##' 
##' @examples
##' plot(tamiasEPM)
##' 
##' plot(tamiasEPM, legend = FALSE)
##' # addLegend(tamiasEPM['spRichness'], location = 'top', ramp=c('blue','yellow','red'))
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
##' minmax <- getMultiMapRamp(list(epm1, epm2))
##' 
##' map1 <- plot(epm1, colorRampRange = log(minmax), log = TRUE, legend = FALSE)
##' map2 <- plot(epm2, colorRampRange = log(minmax), log = TRUE, legend = FALSE)
##' # tmap_arrange(map1, map2)
##' \donttest{
##' # use mapview for plotting
##' plot(tamiasEPM, basemap = 'interactive')
##' }
##'
##' @rdname plot
##' @aliases plot.epmGrid
##' @export

plot.epmGrid <- function(x, log = FALSE, legend = TRUE, col, basemap = 'worldmap', colorRampRange = NULL, singleSpCol = gray(0.9), lwd, borderCol = 'black', alpha = 1, includeFrame = FALSE, use_tmap = TRUE, fastPoints = FALSE, add = FALSE, ...) {
	
	# x = tamiasEPM; log = FALSE; legend = TRUE; basemap = 'worldmap'; colorRampRange = NULL; singleSpCol = gray(0.9); lwd = 0.25; borderCol = 'black'; includeFrame = FALSE; use_tmap = TRUE; alpha = 1
	
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
	if (plotMetric %in% c('range', 'mean_NN_dist', 'min_NN_dist', 'variance', 'disparity', 'rangePCA', 'meanPatristic', 'meanPatristicNN', 'minPatristicNN', 'phyloDisparity', 'PSV', 'PSR')) {
		# determine which cells have just 1 species
		singleSpCells <- singleSpCellIndex(x)
		
		if (inherits(x[[1]], 'sf')) {
			# add on cells that do not have species
			singleSpCells <- sort(c(singleSpCells, which(x[['cellCommInd']] %in% which(sapply(x[['speciesList']], anyNA) == TRUE))))
			grid_multiSp <- x[[1]][- singleSpCells,]
			grid_singleSp <- x[[1]][singleSpCells,]
		} else {
			grid_multiSp <- x[[1]][plotMetric]
			grid_singleSp <- terra::rast(x[[1]][plotMetric])
			grid_multiSp[singleSpCells] <- NA
			grid_singleSp[] <- NA
			grid_singleSp[singleSpCells] <- 1
		}
				
		plotSingleCells <- TRUE
		
	} else {
		plotSingleCells <- FALSE
	}	
	
	
	ncol <- 1000
	isInt <- FALSE
	if (inherits(x[[1]], 'sf')) {
		if (is.integer(sf::st_drop_geometry(x[[1]])[, plotMetric]) & max(x[[1]][[plotMetric]]) <= 10) {
			ncol <- max(x[[1]][[plotMetric]])
			isInt <- TRUE
		}	
	} else if (inherits(x[[1]], 'SpatRaster')) {
		samp <- sample(as.vector(stats::na.omit(terra::values(x[[1]][plotMetric]))), size = 1000, replace = TRUE)
		if (all(as.integer(samp) == samp) & max(terra::minmax(x[[1]][plotMetric])) <= 10) {
			ncol <- max(terra::minmax(x[[1]][plotMetric]))
			isInt <- TRUE
		}
	}

	if (missing(col)) {
		# colramp <- sf::sf.colors
		# default color ramp: viridis turbo, but without the most extreme values
		colramp <- function(n) viridisLite::turbo(n = n, begin = 0.1, end = 0.9)
		# colramp <- grDevices::colorRampPalette(c('blue', 'cyan', 'yellow', 'red'))
	} else {
		if (class(col) == 'function') {
			colramp <- col
		} else if (class(col) == 'character' & length(col) > 1) {
				colramp <- grDevices::colorRampPalette(col)
		} else {
			stop('color not understood.')
		}
	}

	if (missing(lwd)) {
		lwd <- 0.15
	}

					
	if (is.null(colorRampRange)) {
		breaks <- NULL
	} else {
		breaks <- seq(min(colorRampRange), max(colorRampRange), length.out = ncol + 1)
	}
			
		
	metricName <- plotMetric
		
	if (inherits(x[[1]], 'sf')) {
	
	
				
		if (log) {
			x[[1]]$loggedMetric <- log(x[[1]][[plotMetric]])
			metricName <- paste0('log ', plotMetric)
			isInt <- FALSE
			plotMetric <- 'loggedMetric'
		}
		
		# if (!is.null(breaks)) {
			# x[[1]]$color <- colramp(ncol)[findInterval(unlist(sf::st_drop_geometry(x[[1]][, plotMetric])), breaks, all.inside = TRUE)]
		# }
		
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
			
			if (isInt) {
				tmapStyle <- 'pretty'
			} else {
				tmapStyle <- 'cont'
			}
			
			if (!plotSingleCells) {
				map <- map + tmap::tm_shape(x[[1]]) + tmap::tm_fill(plotMetric, palette = colramp(ncol), legend.show = legend, title = metricName, breaks = breaks, style = tmapStyle, midpoint = NA, alpha = alpha) + tmap::tm_borders(col = borderCol, lwd = lwd, alpha = alpha) + tmap::tm_layout(frame = includeFrame, legend.outside = TRUE)
				
			} else {
				
				map <- map + tmap::tm_shape(grid_multiSp, is.master = TRUE, bbox = x[[1]]) + tmap::tm_fill(plotMetric, palette = colramp(ncol), legend.show = legend, title = metricName, breaks = breaks, style = tmapStyle, midpoint = NA, alpha = alpha) + tmap::tm_borders(col = borderCol, lwd = lwd, alpha = alpha) + tmap::tm_shape(grid_singleSp) + tmap::tm_fill(plotMetric, palette = singleSpCol, legend.show = FALSE, alpha = alpha) + tmap::tm_borders(col = borderCol, lwd = lwd, alpha = alpha) + tmap::tm_layout(frame = includeFrame, legend.outside = TRUE)
			
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
				if (!plotSingleCells) {
					valRange <- range(x[[1]][[plotMetric]], na.rm = TRUE)
				} else {
					valRange <- range(grid_multiSp[[plotMetric]], na.rm = TRUE)
				}
				breaks <- seq(min(valRange), max(valRange), length.out = ncol + 1)
			} else {
				breaks <- seq(min(colorRampRange), max(colorRampRange), length.out = ncol + 1)
			}
			
			colors <- colramp(ncol)
			colors <- grDevices::adjustcolor(colors, alpha.f = alpha)
			borderColor <- grDevices::adjustcolor('black', alpha.f = alpha)
			singleSpCol <- grDevices::adjustcolor(singleSpCol, alpha.f = alpha)
			
			# helpful for debugging because much faster than plotting polygons
			if (fastPoints) {
				if (!plotSingleCells) {
					cols <- colors[findInterval(x[[1]][[plotMetric]], breaks, rightmost.closed = TRUE)]
					plot(sf::st_centroid(sf::st_geometry(x[[1]][plotMetric])), pch = 20, reset = FALSE, col = cols, cex = 0.5)
				} else {
					cols <- colors[findInterval(grid_multiSp[[plotMetric]], breaks, rightmost.closed = TRUE)]
					plot(sf::st_centroid(sf::st_geometry(grid_singleSp[plotMetric])), pch = 20, reset = FALSE, col = singleSpCol, cex = 0.5)
					plot(sf::st_centroid(sf::st_geometry(grid_multiSp[plotMetric])), pch = 20, reset = FALSE, col = cols, cex = 0.5, add = TRUE)	
				}
			} else {
	

				if (!plotSingleCells) {
					cols <- colors[findInterval(x[[1]][[plotMetric]], breaks, rightmost.closed = TRUE)]
					plot(sf::st_geometry(x[[1]][plotMetric]), axes = includeFrame, main = NULL, key.pos = NULL, lwd = lwd, reset = FALSE, col = cols, border = borderColor, add = add, ...)
				} else {
					cols <- colors[findInterval(grid_multiSp[[plotMetric]], breaks, rightmost.closed = TRUE)]
					plot(sf::st_geometry(grid_singleSp[plotMetric]), axes = includeFrame, main = NULL, key.pos = NULL, col = singleSpCol, border = NA, reset = FALSE, add = add, ...)
					plot(sf::st_geometry(grid_multiSp[plotMetric]), lwd = lwd, col = cols, border = borderCol, add = TRUE, reset = FALSE)
				}
			



# # 				if (!plotSingleCells) {
					# plot(x[[1]][plotMetric], axes = includeFrame, main = NULL, key.pos = NULL, lwd = lwd, reset = FALSE, pal = colors, border = borderColor, breaks = breaks, add = add, ...)
				# } else {
					# plot(grid_singleSp[plotMetric], axes = includeFrame, main = NULL, key.pos = NULL, pal = singleSpCol, border = NA, reset = FALSE, add = add, ...)
					# plot(grid_multiSp[plotMetric], lwd = lwd, pal = colors, border = borderCol, breaks = breaks, add = TRUE, reset = FALSE)
				# }
			}
			
			if (legend) {
	
				if (!plotSingleCells) {
					if (isInt) {
						if (inherits(x[[1]], 'sf')) {
							nTicks <- max(x[[1]][[plotMetric]])
						} else if (inherits(x[[1]], 'SpatRaster')) {
							nTicks <- max(terra::minmax(x[[1]][plotMetric]))
						}
					} else {
						nTicks <- 3
					}

					addLegend(x[[1]][plotMetric], location = 'topright', ramp = colramp, isInteger = isInt, ncolors = ncol, nTicks = nTicks)

				} else {

					if (isInt) {
						if (inherits(x[[1]], 'sf')) {
							nTicks <- max(grid_multiSp[plotMetric])
						} else if (inherits(x[[1]], 'SpatRaster')) {
							nTicks <- max(terra::minmax(grid_multiSp[plotMetric]))
						}
					} else {
						nTicks <- 3
					}

					addLegend(grid_multiSp[plotMetric], location = 'topright', ramp = colramp, isInteger = isInt, ncolors = ncol, nTicks = nTicks, minmax = range(sf::st_drop_geometry(grid_multiSp[plotMetric])))
				}
			}
	
			if (basemap == 'worldmap') {
				# add map for context
				wrld <- sf::st_transform(worldmap, crs = sf::st_crs(x[[1]]))
				grXY <- graphics::par("usr")
				graphics::clip(grXY[1], grXY[2], grXY[3], grXY[4]) # this ensures that world map is constrained to plot region
				graphics::plot(wrld, add = TRUE, lwd = lwd)
			}
		}

	} else if (inherits(x[[1]], 'SpatRaster')) {
		
			# get bounding box of grid cells with data. We will do this so that we can leave out of the plot those grid cells that are empty and that may take up quite a bit of geographic space. Only an issue for rasters.
			datBB <- sf::st_bbox(sf::st_as_sf(as.data.frame(terra::xyFromCell(x[[1]], which(!is.na(terra::values(x[[1]][[plotMetric]]))))), coords = 1:2, crs = sf::st_crs(x[[1]])))
		
				
		if (log) {
			metricName <- paste0('log ', plotMetric)
			if (!plotSingleCells) {
				metricMap <- log(x[[1]][plotMetric])
			} else {
				metricMap <- log(grid_multiSp[plotMetric])
			}
			isInt <- FALSE
		} else {
			if (!plotSingleCells) {
				metricMap <- x[[1]][plotMetric]
			} else {
				metricMap <- grid_multiSp[plotMetric]
			}
		}
		
		if (use_tmap) {
		    
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
			
#			if (isInt) {
#				tmapStyle <- 'pretty'
#			} else {
				tmapStyle <- 'cont'
#			}
			
			if (!plotSingleCells) {
				map <- map + tmap::tm_shape(metricMap, bbox = datBB) + tmap::tm_raster(palette = colramp(ncol), legend.show = legend, title = metricName, breaks = breaks, style = tmapStyle, midpoint = NA, alpha = alpha) + tmap::tm_layout(frame = includeFrame, legend.outside = TRUE)
				
			} else {
				
				map <- map + tmap::tm_shape(metricMap, is.master = TRUE, bbox = datBB) + tmap::tm_raster(palette = colramp(ncol), legend.show = legend, title = metricName, breaks = breaks, style = tmapStyle, midpoint = NA, alpha = alpha) + tmap::tm_shape(grid_singleSp) + tmap::tm_raster(palette = singleSpCol, legend.show = FALSE, alpha = alpha) + tmap::tm_layout(frame = includeFrame, legend.outside = TRUE)
			
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
			
		
			if (!plotSingleCells) {
				terra::plot(terra::crop(metricMap, datBB2), col = colramp(ncol), axes = includeFrame, legend = FALSE, plg = list(shrink = 0.7, title = metricName), range = valRange, alpha = alpha, add = add, ...)
				
			} else {
				
				terra::plot(terra::crop(metricMap, datBB2), col = colramp(ncol), axes = includeFrame, legend = FALSE, plg = list(shrink = 0.7, title = metricName), range = valRange, alpha = alpha, add = add, ...)
				terra::plot(grid_singleSp, col = singleSpCol, axes = includeFrame, legend = FALSE, range = valRange, add = TRUE, alpha = alpha)
				
			}
			
			if (legend) {
	
				if (!plotSingleCells) {
					if (isInt) {
						if (inherits(x[[1]], 'sf')) {
							nTicks <- max(x[[1]][[plotMetric]])
						} else if (inherits(x[[1]], 'SpatRaster')) {
							nTicks <- max(terra::minmax(x[[1]][plotMetric]))
						}
					} else {
						nTicks <- 3
					}

					addLegend(x[[1]][plotMetric], location = 'topright', ramp = colramp, isInteger = isInt, ncolors = ncol, nTicks = nTicks)

				} else {

					if (isInt) {
						if (inherits(x[[1]], 'sf')) {
							nTicks <- max(grid_multiSp[plotMetric])
						} else if (inherits(x[[1]], 'SpatRaster')) {
							nTicks <- max(terra::minmax(grid_multiSp[plotMetric]))
						}
					} else {
						nTicks <- 3
					}

					addLegend(grid_multiSp[plotMetric], location = 'topright', ramp = colramp, isInteger = isInt, ncolors = ncol, nTicks = nTicks)
				}
			}
			
	
			if (basemap == 'worldmap') {
				# add map for context
				wrld <- sf::st_transform(worldmap, crs = sf::st_crs(x[[1]]))
				grXY <- graphics::par("usr")
				graphics::clip(grXY[1], grXY[2], grXY[3], grXY[4]) # this ensures that world map is constrained to plot region
				graphics::plot(wrld, add = TRUE, lwd = lwd)
			}
		}
	}
}
		
		
		
		
		
		
		
		
		
		
		
		
	
# library(ggplot2)
# gg <- ggplot(spRas[[1]]) + geom_sf(aes(fill = spRichness), color = NA) + scale_fill_viridis_c(option = "plasma")
# gg

		
		
			
		
		
		
		
		
		
		
# plot(grid_singleSp[plotMetric], extent = sf::st_bbox(x[[1]]), axes = includeFrame, main = NULL, pal = singleSpCol, border = NA, reset = FALSE)
# plot(grid_multiSp[plotMetric], key.pos = key.pos, lwd = lwd, reset = FALSE, pal = colramp(ncol), breaks = breaks, add = TRUE)

# use tmap_arrange() for multi-plots

# diy with legend: show in example how to invisibly return and add legend where you want it


# make basemap light gray, add graticules behind it, add arguments for basemap color, and true/false graticules




# tmap::tm_shape(x[[1]]) + tmap::tm_fill(col = 'color', title = metricName, border.col = borderCol, lwd = lwd) + tmap::tm_layout(frame = includeFrame, legend.outside = TRUE)


# tmap::tm_shape(x[[1]]) + tmap::tm_fill(plotMetric, palette = 'viridis', legend.show = legend, title = metricName, breaks = breaks) + tmap::tm_borders(col = borderCol, lwd = lwd) + tmap::tm_layout(frame = includeFrame)

# tmap::tm_shape(x[[1]]) + tmap::tm_fill('loggedMetric', palette = colramp(ncol), legend.show = TRUE, title = metricName, breaks = seq(min(colorRampRange), max(colorRampRange), length.out = 5), style = 'cont', as.count = isInt) + tm_layout(legend.only = TRUE)



