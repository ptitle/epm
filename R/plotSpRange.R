##' @title plot a single species' range
##'
##' @description Plot one species' geographic range, as encoded in the epmGrid object.
##'
##' @param x object of class \code{epmGrid}
##' @param taxon taxon to plot
##' @param taxonColor color for plotting taxon's range
##' @param basemap if \code{'none'}, then only the grid is plotted. 
##'		If \code{'worldmap'}, then vector map is plotted.
##' 	If \code{'interactive'}, then the plotting is done via your web browser.
##' @param lwd grid cell border width
##' @param alpha opacity of all colors and borders, ranging from 0 
##' (fully transparent) to 1 (fully opaque)
##' @param use_tmap if false, use sf or terra packages for plotting
##' @param add logical. If TRUE, adds the gridded taxon range to existing plot. 
##'
##' @return nothing is returned
##'
##' @author Pascal Title
##'
##' @examples
##' tamiasEPM
##'
##' plotSpRange(tamiasEPM, 'Tamias_speciosus')
##' 
##' @export




# plot one species' geographic range, as encoded in the epmGrid object.

plotSpRange <- function(x, taxon, taxonColor = 'orange', basemap = 'worldmap', lwd = 0.5, alpha = 1, use_tmap = TRUE, add = FALSE) {
	
	if (use_tmap & !requireNamespace('tmap', quietly = TRUE)) {
		message('\ttmap package not installed -- defaulting to sf/terra plot')
		use_tmap <- FALSE
	}

	if (!inherits(x, 'epmGrid')) {
		stop('x must be of class epmGrid.')
	}

	if (!taxon %in% x[['geogSpecies']]) {
		stop('taxon not found.')
	}
	
	if (add) {
		use_tmap <- FALSE
		basemap <- 'none'
	}
	
	comm <- which(sapply(x[['speciesList']], function(x) taxon %in% x) == TRUE)
	cellInd <- which(x[['cellCommInd']] %in% comm)
	
	if (use_tmap) {

		if (basemap != 'interactive') {
		
			if (basemap == 'worldmap') {
				map <- tmap::tm_shape(worldmap, is.master = FALSE, projection = sf::st_crs(x[[1]]), bbox = x[[1]]) + tmap::tm_borders(lwd = 0.3)
			} else {
				map <- NULL
			}
			
			if (inherits(x[[1]], 'sf')) {
				map <- map + tmap::tm_shape(x[[1]]) + tmap::tm_polygons(col = gray(0.95), border.col = gray(0.90), legend.show = FALSE, alpha = alpha) + tmap::tm_layout(frame = FALSE) + tmap::tm_shape(x[['grid']][cellInd, ]) + tmap::tm_polygons(col = taxonColor, border.col = 'black', lwd = lwd, alpha = alpha)
	
			} else {
				
				datBB <- sf::st_bbox(sf::st_as_sf(as.data.frame(terra::xyFromCell(x[[1]], which(!is.na(terra::values(x[[1]][['spRichness']]))))), coords = 1:2, crs = sf::st_crs(x[[1]])))
	
				fullGrid <- x[[1]]['spRichness']
				fullGrid[!is.na(fullGrid)] <- 1
				spGrid <- terra::rast(fullGrid)
				spGrid[cellInd] <- 1
				map <- map + tmap::tm_shape(fullGrid, bbox = datBB) + tmap::tm_raster(pal = gray(0.95), legend.show = FALSE, alpha = alpha) + tmap::tm_layout(frame = FALSE) + tmap::tm_shape(spGrid) + tmap::tm_raster(pal = taxonColor, alpha = alpha, legend.show = FALSE)
			}
			map
		
		} else {
		
			tmap::tmap_mode('view')
			
			if (inherits(x[[1]], 'sf')) {
				map <- tmap::tm_shape(x[[1]]) + tmap::tm_polygons(col = gray(0.95), border.col = gray(0.90), alpha = alpha, legend.show = FALSE) + tmap::tm_layout(frame = FALSE) + tmap::tm_shape(x[['grid']][cellInd, ]) + tmap::tm_polygons(col = taxonColor, border.col = 'black', lwd = lwd, alpha = alpha)
	
			} else {
	
				fullGrid <- x[[1]]['spRichness']
				fullGrid[!is.na(fullGrid)] <- 1
				spGrid <- terra::rast(fullGrid)
				spGrid[cellInd] <- 1
				map <- tmap::tm_shape(fullGrid) + tmap::tm_raster(pal = gray(0.95), alpha = alpha, legend.show = FALSE) + tmap::tm_shape(spGrid) + tmap::tm_raster(pal = taxonColor, alpha = alpha, legend.show = FALSE)
			}
	
			map
		}
		
	} else {
		
		# don't use t_map
		
		if (inherits(x[[1]], 'sf')) {
			
			if (!add) {
				plot(sf::st_combine(sf::st_geometry(x[['grid']])), col = gray(0.95, alpha = alpha), lwd = lwd, border = NA)
			}
			
			plot(sf::st_geometry(x[['grid']])[cellInd], add = TRUE, lwd = lwd, col = grDevices::adjustcolor(taxonColor, alpha.f = alpha), border = grDevices::adjustcolor('black', alpha.f = alpha))
				
		} else {
			
			spGrid <- terra::rast(x[['grid']][[1]])
			spGrid[cellInd] <- 1
			if (!add) {
				terra::plot(x[['grid']]['spRichness'], col = gray(0.95, alpha = alpha), axes = FALSE, legend = FALSE)
			}
			terra::plot(spGrid, col = grDevices::adjustcolor(taxonColor, alpha.f = alpha), add = TRUE, axes = FALSE, legend = FALSE)
				
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





