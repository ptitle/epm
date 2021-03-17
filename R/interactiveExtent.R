##' @title Interactively choose extent
##'
##' @description Given a list of polygons, sets up an interactive plot
##' 	to allow the user to draw the desired extent. This can be used to
##' 	define the extent in \link{createEPMgrid}.
##'
##' @param polyList a list of Simple Feature polygons.
##'
##' @return A list with a polygon, and its WKT string 
##'
##' @author Pascal Title
##' 
##' @examples
##' \dontrun{
##' ex <- interactiveExtent(tamiasPolyList)
##'
##' # You can use this as the extent in createEPMgrid
##' grid <- createEPMgrid(tamiasPolyList, resolution = 50000, extent = ex$wkt)
##'
##'}
##' @export




interactiveExtent <- function(polyList) {
	
	#get overall extent
	masterExtent <- sf::st_as_sfc(getExtentOfList(polyList, format = 'sf'))
	
	# coarse template
	# if projected, use 100km, if not, use 20 degrees
	quickRes <- ifelse(sf::st_is_longlat(polyList[[1]]), 20, 100000)

	quickTemplate <- sf::st_make_grid(masterExtent, cellsize = c(quickRes, quickRes), square = FALSE, crs = sf::st_crs(masterExtent))
	quickCentroids <- sf::st_centroid(quickTemplate)
	quickTemplate <- sf::st_sf(quickTemplate, grid_id = 1:length(quickTemplate))
	
	quick <- lapply(polyList, function(x) polyToGridCells(x, method ='centroid', quickTemplate, quickCentroids, subset = 0))

	# flip list from list of cells per species, to list of species per cells
	cellList <- vector('list', nrow(quickTemplate))
	for (i in 1:length(quick)) {
		cellList[quick[[i]]] <- lapply(cellList[quick[[i]]], function(x) append(x, names(quick)[i]))
	}
	cellList <- lapply(cellList, unique)

	# uniqueComm <- unique(cellList)
	# uniqueComm2 <- sapply(uniqueComm, function(y) paste(y, collapse = '|'))
	# cellList2 <- sapply(cellList, function(y) paste(y, collapse = '|'))

	# cellCommVec <- integer(length = length(cellList))
	# for (i in 1:length(uniqueComm)) {
		# cellCommVec[which(cellList2 == uniqueComm2[i])] <- i
	# }
	
	# rm(cellList2, uniqueComm2)

	# add unique community index in as attribute
	# quickTemplate$uniqueComm <- cellCommVec
	quickTemplate$spRichness <- lengths(cellList)
	
	# remove empty gridcells
	quickTemplate <- quickTemplate[which(lengths(cellList) > 0),]

	# add map for context
	if (!sf::st_is_longlat(polyList[[1]])) {
		wrld <- sf::st_transform(worldmap, crs = sf::st_crs(polyList[[1]]))
	} else {
		wrld <- worldmap
	}
	
	message('\n\tAn interactive coarse-grain map has been displayed.\n')
	message('\n\tPlease wait until plot is completed......', appendLF = FALSE)
	
	plot(quickTemplate['spRichness'], border = NA, key.pos = NULL, main = NULL, reset = FALSE)
	graphics::plot(wrld, add = TRUE, lwd = 1, border = gray(0.5))
	
	#ggplot(quickTemplate) + geom_sf(aes(fill = spRichness), lwd = 0.5, show.legend = FALSE) + scale_fill_gradientn(colours = sf::sf.colors()) + geom_sf(data = wrld, fill = NA)
	
	message('done!\n')
	graphics::title(main = 'Define your extent polygon.')

	message('\tClick on the map to create a polygon that will define the extent of the epmGrid.')
	message('\tRight-clicking will close the polygon and terminate the interactive plot.\n\n')
	
	userPoly <- terra::draw(x = 'polygon', col = 'red', xpd = NA)
	userPoly <- terra::as.data.frame(userPoly, geom = "hex")
	userPoly$geometry <- structure(as.list(userPoly$geometry), class = "WKB")
	userPoly <- sf::st_as_sf(userPoly)
	sf::st_crs(userPoly) <- sf::st_crs(polyList[[1]])
	
	# display call so user can use this extent in the future
	wkt <- sf::st_as_text(sf::st_geometry(userPoly))
	
	grDevices::dev.off()
	
	return(list(poly = userPoly, wkt = wkt))
}