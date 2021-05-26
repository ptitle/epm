##' @title Create epmGrid object
##' 
##' @description Takes a list of polygons or point occurrences and creates a epmGrid object.
##'
##'
##' @param spDat a list of polygon objects (sf or sp), named with taxon names. Alternatively, 
##' a set of occurrence records. See details for more information.    
##'
##' @param resolution vertical and horizontal spacing of grid cells, in units 
##'		of the polygons' or points' projection.
##'
##' @param method approach used for gridding. Either \code{centroid} or \code{areaCutoff}. See details below.
##'
##' @param cellType either \code{hexagon} or \code{square}. See details below.
##' 
##' @param coverCutoff the percent that a species range must cover a grid cell to be considered present.
##'
##' @param retainSmallRanges boolean; should small ranged species be dropped or preserved.
##'		See details.
##'
##' @param extent if 'auto', then the maximal extent of the polygons will be used. 
##' 	If not 'auto', can be a SpatialPolygon, sf object, or raster, in which case the resulting epmGrid
##'		will be cropped and masked with respect to the polygon, or a spatial coordinates object, 
##' 	from which an extent object will be generated, or a numeric vector of length 4 
##' 	with minLong, maxLong, minLat, maxLat. If 'interactive', then an interactive plot
##'		will appear in which the user can draw the desired polygon extent. That extent will then be returned.
##'
##' @param percentWithin The percentage of a species range that must be within the defined extent in order
##' 	for that species to be included. This filter can be used to exclude species whose range barely enters
##' 	the area of interest. The default value of 0 will disable this filter. If \code{extent == 'auto'},
##'		then this filter will also have no effect, as the extent is defined by the species' ranges.
##'
##' @param checkValidity if \code{TRUE}, then check polygon validity and repair if needed, 
##' 	using sf::st_make_valid and the lwgeom package. 
##'
##' @param crs if supplying occurrence records in a non-spatial format, then you must specify the crs.
##' For unprojected long/lat data, you can simply provide \code{crs = 4326}.
##'
##' @param nThreads if > 1, then employ parallel computing. 
##'
##' @param template an object of class \code{SpatRaster} or \code{RasterLayer} that can be used to 
##' 	get extent and resolution. If \code{cellType = 'square'}, then the template will be used as the 
##' 	reference grid.
##' 
##' @param nGroups break up the grid into this many groups for processing. This will alleviate memory 
##' 	usage for datasets that are very high resolution (not relevant for large numbers of species).
##'
##' @param verbose if TRUE, list out all species that are dropped/excluded, rather than counts.
##'
##'
##' @details 
##'		If \code{cellType = 'hexagon'}, then the grid is made of polygons via the sf package.
##' 	If \code{cellType = 'square'}, then the grid is a raster generated via the terra package.
##' 	Hexagonal cells have several advantages, including being able to be of different sizes (if the grid is in 
##' 	unprojected long/lat), and may be able to more naturally follow coastlines and non-linear features.
##'		However, the raster-based square cells will be much less memory intensive for high resolution datasets. 
##' 	Choice of grid type matters more for spatial resolution (total number of cells), than for number of species.
##'
##' 
##' 	In the polygon-to-grid conversion process, two approaches are implemented. 
##'		For \code{method = 'centroid'}, a range polygon registers in a cell if the polygon 
##' 	overlaps with the cell centroid. 
##' 	For \code{method = 'areaCutoff'}, a range polygon registers in a cell if it covers that cell by 
##' 	at least \code{coverCutoff} fraction of the cell.
##' 	If \code{retainSmallRanges = FALSE}, then species whose ranges are so small that no 
##' 	cell registers as present will be dropped. If \code{retainSmallRanges = TRUE}, then the 
##' 	cell that contains the majority of the the small polygon will be considered as present, even if it's a 
##' 	small percent of the cell.
##' 
##'		If \code{retainSmallRanges = TRUE}, and an extent is provided, then species may still be dropped if they
##' 	fall outside of that extent.
##'
##'		In interactive mode for defining the extent, the user can draw a bounding polygon on a 
##'		map. The drawn polygon will then be printed to the console so that the user can provide 
##'		that bounding polygon in future calls as the extent.
##'		
##'		If input data consist of occurrence records rather than polygons, then a couple of formats are possible:
##'     \enumerate{
##'		  \item You can provide a list of species-specific spatial point objects. 
##'	 	  \item You can provide a single spatial object, where points have a taxon attribute.
##'		  \item You can provide a list of non-spatial species-specific dataframes.
##'		  \item You can provide a single non-spatial dataframe.
##'		}
##'		
##'		For options (1) and (3), the taxon names must be provided as the list names.
##'		For options (3) and (4), the columns must be 'taxon', 'x' and 'y' (or 'long', 'lat').
##'		For options (3) and (4), as these are non-spatial, you must provide a crs object to the
##'		\code{crs} argument, so that the function knows what projection to use.
##'
##'		Any SpatialPolygon or SpatialPoints objects are converted to objects of class \code{sf}.
##'	
##'     For \code{extent = 'interactive'}, you can additionally specify some bounding coordinates.
##'     This can be helpful if the interactive map is too broad in extent, making it difficult to draw
##'     the extent that you want. Instead, specify \code{extent = list('interactive', c(xmin, xmax, ymin, ymax))}.
##' 	In interactive mode, the basemap is from \url{www.naturalearthdata.com}. 
##' 
##' @return an object of class \code{epmGrid}. If \code{extent = 'interactive'}, then a polygon is returned.
##' 
##' @author Pascal Title
##'
##' @examples
##' library(sf)
##' # example dataset: a list of 24 chipmunk distributions as polygons
##' head(tamiasPolyList)
##' 
##' # hexagonal grid
##' tamiasEPM <- createEPMgrid(tamiasPolyList, resolution = 50000, 
##' 	cellType = 'hexagon', method = 'centroid')
##' tamiasEPM
##'
##' # square grid
##' tamiasEPM2 <- createEPMgrid(tamiasPolyList, resolution = 50000, 
##' 	cellType = 'square', method = 'centroid')
##' tamiasEPM2
##' 
##' #######
##' # With point occurrences
##' ## demonstrating all possible input formats
##' 
##' # list of sf spatial objects
##' spOccList <- lapply(tamiasPolyList, function(x) st_sample(x, size = 10, type= 'random'))
##' tamiasEPM <- createEPMgrid(spOccList, resolution = 100000, cellType = 'hexagon')
##'
##' # list of coordinate tables
##' spOccList2 <- lapply(spOccList, function(x) st_coordinates(x))
##' tamiasEPM <- createEPMgrid(spOccList, resolution = 100000, cellType = 'square')
##' 
##' # single table of coordinates
##' spOccList3 <- spOccList2
##' for (i in 1:length(spOccList3)) {
##' 	spOccList3[[i]] <- cbind.data.frame(taxon = names(spOccList3)[i], spOccList3[[i]])
##' 	colnames(spOccList3[[i]]) <- c('taxon', 'X', 'Y')
##' }
##' spOccList3 <- do.call(rbind, spOccList3)
##' rownames(spOccList3) <- NULL
##' spOccList3[, "taxon"] <- as.character(spOccList3[, "taxon"])
##' tamiasEPM <- createEPMgrid(spOccList, resolution = 100000, cellType = 'square')
##'
##' # a single labeled spatial object
##' spOccList4 <- st_as_sf(spOccList3[, c("taxon", "X", "Y")], coords = c("X","Y"), 
##' crs = st_crs(spOccList[[1]]))
##' tamiasEPM <- createEPMgrid(spOccList, resolution = 100000, cellType = 'square')
##' 
##' 
##' 
##' @export

createEPMgrid <- function(spDat, resolution = 50000, method = 'centroid', cellType = 'hexagon', coverCutoff = 0.1, retainSmallRanges = TRUE, extent = 'auto', percentWithin = 0, checkValidity = FALSE, crs = NULL, nThreads = 1, template = NULL, nGroups = 1, verbose = FALSE) {
	
	# spDat <- tamiasPolyList; resolution = 50000; method = 'centroid'; cellType = 'hexagon'; coverCutoff = 0.1; retainSmallRanges = TRUE; extent = 'auto'; percentWithin = 0; checkValidity = FALSE; nThreads = 1; template = NULL; nGroups = 1; verbose = FALSE
	
	# test with occurrences
	# spOccList <- lapply(tamiasPolyList, function(x) sf::st_sample(x, size = 10, type= 'random'))
	# spDat <- spOccList; resolution = 50000; method = 'centroid'; cellType = 'hexagon'; coverCutoff = 0.1; retainSmallRanges = TRUE; extent = 'auto'; checkValidity = FALSE; nThreads = 1; template = NULL; nGroups = 1; verbose = FALSE
	
    if (inherits(spDat, 'Spatial')) {
        # if class SpatialPolygons, convert to sf
        spDat <- sf::st_as_sf(spDat)
    }
    
    if (!inherits(spDat, 'Spatial') & inherits(spDat[[1]], 'Spatial')) {
        # if class SpatialPolygons, convert to sf
        for (i in 1:length(spDat)) {
            spDat[[i]] <- sf::st_as_sf(spDat[[i]])
        }
    }

	# determine whether input is occurrences or polygons
    # and if occurrences, get crs
    
    # list of sf objects (only acceptable format for polygons)
    if (!inherits(spDat, c('sf', 'sfc')) & inherits(spDat[[1]], c('sf', 'sfc'))) {
        if (unique(as.character(sf::st_geometry_type(spDat[[1]]))) %in% c('MULTIPOLYGON', 'POLYGON', 'GEOMETRY')) {
            datType <- 'polygons'
        } else if (unique(as.character(sf::st_geometry_type(spDat[[1]]))) == 'POINT') {
            datType <- 'points'
            proj <- sf::st_crs(spDat[[1]])
        } else {
        	stop('Geometry type not recognized.')
        }
        
    # single sf object
    } else if (inherits(spDat, c('sf', 'sfc')) & !inherits(spDat[[1]], c('sf', 'sfc'))) {
        datType <- 'points'
        proj <- sf::st_crs(spDat)
    
    # single non-spatial data table
    } else if (!inherits(spDat, c('sf', 'sfc')) & inherits(spDat, c('data.frame', 'matrix'))) {
        datType <- 'points'
        if (!is.null(crs)) {
            proj <- crs
        } else {
            stop('Input occurrence data need a crs but none provided.')
        }
    
    # list of non-spatial tables
    } else if (!inherits(spDat, c('data.frame', 'matrix')) & inherits(spDat[[1]], c('data.frame', 'matrix')) & !inherits(spDat[[1]], c('sf', 'sfc'))) {
        datType <- 'points'
        if (!is.null(crs)) {
            proj <- crs
        } else {
            stop('Input occurrence data need a crs but none provided.')
        }
        
    } else {
        stop('Format of spDat not recognized.')
    }
        
    if (datType == 'points') {
        percentWithin <- 0
        method <- 'centroid'
        retainSmallRanges <- FALSE
        nGroups <- 1
       
        # reformat occurrence data if necessary
        ## all input formats are returned as a list of species-specific tables
        occList <- occurrenceFormatting(spDat)
      
        # create list of sf objects
        spDat <- lapply(occList, function(x) sf::st_as_sf(x, coords = 2:3, crs = proj))
       
        message('\t', 'Detected ', length(occList), ' taxa with point data.')
    } 
    
    
	cellType <- match.arg(cellType, c('hexagon', 'square'))
	
	if (!cellType %in% c('hexagon', 'square')) {
		stop('cellType must be either hexagon or square.')
	}

	if (is.null(names(spDat))) {
		stop('List must be named with species names.')
	}
	
	if (anyDuplicated(names(spDat))) {
		stop('Some taxon names are duplicated.')
	}
	
	if (coverCutoff > 1) {
		stop('coverCutoff is a fraction.')
	}
	
	if (nThreads > parallel::detectCores()) {
		stop('nThreads is greater than what is available.')
	}
	
	if (inherits(extent, c('SpatRaster', 'RasterLayer'))) {
		if (inherits(extent, 'RasterLayer')) {
			extent <- as(extent, 'SpatRaster')
		}
	    if (sum(c(terra::is.lonlat(extent, perhaps = TRUE), sf::st_is_longlat(spDat[[1]]))) == 1) {
	        stop('Raster provided as extent has a different projection from input data.')
	    }
		extent <- as.vector(terra::ext(extent))
	}
	
	if (!is.null(template)) {
		if (!inherits(template, c('SpatRaster', 'RasterLayer'))) {
			stop('template must be a RasterLayer or SpatRaster.')
		}
		if (inherits(template, 'RasterLayer')) {
			template <- as(template, 'SpatRaster')
		}
	    if (sum(c(terra::is.lonlat(template, perhaps = TRUE), sf::st_is_longlat(spDat[[1]]))) == 1) {
	        stop('Raster provided as template has a different projection from input data.')
	    }
	    resolution <- terra::res(template)
		extent <- as.vector(terra::ext(template))
        if (cellType == 'hexagon') {
            stop('Use of the template argument is intended for square-cell grids only.')
        }
	}
	
	# if lwgeom package available, and checking is requested, check sf polygon validity. 
	if (!requireNamespace('lwgeom', quietly = TRUE) & checkValidity) {
		stop('For polygon validity checking, the lwgeom package is required.')
	}
	
	if (requireNamespace('lwgeom', quietly = TRUE) & checkValidity) {
		for (i in 1:length(spDat)) {
			if (inherits(spDat[[i]], 'sf')) {
				if (!any(sf::st_is_valid(spDat[[i]]))) {
					message('\trepairing polygon ', i)
					spDat[[i]] <- sf::st_make_valid(spDat[[i]])
				}
			}
		}
	}

	# test that all have same CRS
	if (!any(sapply(spDat, function(x) identical(sf::st_crs(spDat[[1]]), sf::st_crs(x))))) {
		stop('proj4string or EPSG of all polygons must match.')
	}

	if (sf::st_is_longlat(spDat[[1]])) {
	    proj <- sf::st_crs(4326)
	} else {
	    proj <- sf::st_crs(spDat[[1]])
	}
	
	# if WKT string, then convert to sf polygon
	if (inherits(extent, 'character')) {
		if (grepl('POLYGON \\(', extent)) {
			extent <- sf::st_as_sfc(extent, crs = sf::st_crs(spDat[[1]]))
		}
	}

	if (inherits(extent, 'character') | 'interactive' %in% unlist(extent)) {
		
		if (!any(extent %in% c('auto', 'interactive'))) {
			stop("If extent is a character vector, it can only be 'auto' or 'interactive'.")
		}
		
		if (all(extent == 'auto')) {
		
			#get overall extent
			masterExtent <- getExtentOfList(spDat, format = 'sf')
			percentWithin <- 0
		}
	} else if (is.numeric(extent) & length(extent) == 4) {
		# use user-specified bounds
		masterExtent <- sf::st_bbox(spDat[[1]])
		masterExtent[[1]] <- extent[1]
		masterExtent[[2]] <- extent[3]
		masterExtent[[3]] <- extent[2]
		masterExtent[[4]] <- extent[4]
		
	} else if (inherits(extent, c('SpatialPolygons', 'SpatialPolygonsDataFrame', 'SpatialPoints', 'SpatialPointsDataFrame', 'sf', 'sfc'))) {
		
		if (inherits(extent, c('SpatialPolygons', 'SpatialPolygonsDataFrame', 'SpatialPoints', 'SpatialPointsDataFrame'))) {
			extent <- sf::st_as_sf(extent)
		}
			
		if (!is.na(sf::st_crs(extent))) {
			if (!identical(sf::st_crs(extent), proj)) {
				extent <- sf::st_transform(extent, crs = proj)
			}
		} else {
			sf::st_crs(extent) <- proj
		}
	  
		if (inherits(extent, 'sfc')) {
			extent <- sf::st_sf(extent)
		}
		
		# get extent from spatial object
		# masterExtent <- sf::st_bbox(extent)
		masterExtent <- extent

	} else if (inherits(extent, 'Extent')) {
		
		masterExtent <- sf::st_bbox(spDat[[1]])
		masterExtent[[1]] <- extent@xmin
		masterExtent[[2]] <- extent@ymin
		masterExtent[[3]] <- extent@xmax
		masterExtent[[4]] <- extent@ymax
		
	} else {
		stop("extent must be 'auto', a spatial object or a vector with minLong, maxLong, minLat, maxLat.")
	}
	
	# interactive extent: if this option is selected, a coarse richness grid
	# will be plotted, so that the user can designate an extent polygon
	wkt <- NULL
	if (inherits(extent, 'character') | 'interactive' %in% unlist(extent)) {
        if (!all(extent == 'auto')) {
	    	if (all(extent == 'interactive')) {
    			interactive <- interactiveExtent(spDat, nThreads = nThreads)
    		} else if ('interactive' %in% unlist(extent)) {
    		    bb <- extent[[which(sapply(extent, function(x) !all(x == 'interactive')) == TRUE)]]
    		    if (length(bb) != 4 | !inherits(bb, 'numeric')) {
    		        stop("If you are trying to provide bounding coordinates for the interactive extent, then something is wrong.\n You should specify extent = list('interactive', c(xmin, xmax, ymin, ymax))")
    		    }
    		    interactive <- interactiveExtent(spDat, bb = bb, nThreads = nThreads)
    		}
    		extent <- interactive$poly
    		wkt <- interactive$wkt
    			
    		return(wkt)
        }
    }
	
	# percentWithin: Implement a filter based on the percentage that each species' range overlaps the extent.
	## This allows us to provide some threshold for species inclusion. For instance, 5% would imply that a 
	## species must have at least 5% of its range within the extent to be considered. 
	## Default will be 0, indicating that there is no filtering.
	if (percentWithin > 0) {
		
		if (inherits(masterExtent, 'sf')) {
			extentPoly <- masterExtent
		} else {
			extentPoly <- sf::st_as_sfc(sf::st_bbox(sf::st_sf(geom = sf::st_sfc(sf::st_point(masterExtent[c('xmin', 'ymin')]), sf::st_point(masterExtent[c('xmax', 'ymax')])), crs = proj)))
		}
		
		includeSp <- logical(length(spDat))
		for (i in 1:length(spDat)) {
			overlap <- sf::st_intersection(sf::st_geometry(spDat[[i]]), extentPoly)
			if (as.numeric(sum(sf::st_area(overlap)) / sum(sf::st_area(sf::st_geometry(spDat[[i]])))) >= percentWithin) {
				includeSp[i] <- TRUE
			} else {
				includeSp[i] <- FALSE
			}
		}
		
		if (any(includeSp == FALSE)) {
		    if (verbose) {
		        msg <- paste0('The following species were excluded due to the percentWithin filter:\n\t', paste(names(spDat)[which(includeSp == FALSE)], collapse = '\n\t'))
		    } else {
		        msg <- paste0('\n\t', sum(includeSp == FALSE), ' species ', ifelse(sum(includeSp == FALSE) == 1, 'was', 'were'), ' excluded due to the percentWithin filter.')
		    }
		    message(msg)
		}
		
		# exclude those species that did not satisfy the filter
		spDat <- spDat[includeSp]
	}
	
	
	# here, we take one of two routes:
	## if hexagonal grid, then use sf, if square grid, use terra
	
	uniqueSp <- sort(names(spDat))
	spDat <- spDat[uniqueSp]

	
    if (cellType == 'hexagon') {
		
    	# poly = spDat; method = method; coverCutoff = coverCutoff; extentVec = masterExtent; resolution = resolution; crs = proj; nGroups = nGroups; retainSmallRanges = retainSmallRanges; nThreads = nThreads
    	spGridList <- polyToHex(poly = spDat, method = method, coverCutoff = coverCutoff, extentVec = masterExtent, resolution = resolution, crs = proj, nGroups = nGroups, retainSmallRanges = retainSmallRanges, nThreads = nThreads)
		
    } else if (cellType == 'square') {
		
    	spGridList <- polyToTerra(poly = spDat, method = method, coverCutoff = coverCutoff, extentVec = masterExtent, resolution = resolution, crs = proj, retainSmallRanges = retainSmallRanges, template = template, nThreads = nThreads)
		
    } else {
    	stop('Cell type not supported.')
    }
	
    gridTemplate <- spGridList[[1]]
    spGridList <- spGridList[[2]]
  
	# update
	smallSp <- which(lengths(spGridList) == 0)
	
	# if small ranged species are not to be preserved, drop them
	# OR if some species are outside of proposed extent, drop them
	if (length(smallSp) > 0) {
		spGridList <- spGridList[ - smallSp]
        if (verbose) {
    		msg <- paste0('\nThe following species are being dropped:\n\t', paste(names(smallSp), collapse = '\n\t'))
        } else {
            msg <- paste0('\n', length(smallSp), ' species ', ifelse(length(smallSp) == 1, 'is', 'are'), ' being dropped.\n')
        }
    	message(msg)
	}
	
	# flip list from list of cells per species, to list of species per cells
	if (inherits(gridTemplate, 'SpatRaster')) {
		cellList <- vector('list', terra::ncell(gridTemplate))
	} else if (inherits(gridTemplate, 'sf')) {
		cellList <- vector('list', nrow(gridTemplate))
	} else {
		stop('gridTemplate class mismatch.')
	}
	for (i in 1:length(spGridList)) {
		cellList[spGridList[[i]]] <- lapply(cellList[spGridList[[i]]], function(x) append(x, names(spGridList)[i]))
	}
	cellList <- lapply(cellList, unique)
	
	uniqueComm <- unique(cellList)
	uniqueComm2 <- sapply(uniqueComm, function(y) paste(y, collapse = '|'))
	cellList2 <- sapply(cellList, function(y) paste(y, collapse = '|'))

	cellCommVec <- integer(length = length(cellList))
	for (i in 1:length(uniqueComm)) {
		cellCommVec[which(cellList2 == uniqueComm2[i])] <- i
	}
	
	rm(cellList2, uniqueComm2)
	
	# set empty cells to NA, rather than NULL
	uniqueComm[sapply(uniqueComm, is.null)] <- NA

	# add unique community index and species richness in as attributes
	if (inherits(gridTemplate, 'SpatRaster')) {
		gridTemplate <- c(gridTemplate, terra::rast(gridTemplate))
		terra::values(gridTemplate) <- cbind(cellCommVec, lengths(cellList))
		names(gridTemplate) <- c('uniqueComm', 'spRichness')
		gridTemplate$spRichness[gridTemplate$spRichness == 0] <- NA
	} else if (inherits(gridTemplate, 'sf')) {
		gridTemplate$uniqueComm <- cellCommVec
		gridTemplate$spRichness <- lengths(cellList)
	} else {
		stop('gridTemplate class mismatch.')
	}
	
	# remove empty gridcells
	if (inherits(gridTemplate, 'sf')) {
		gridTemplate <- gridTemplate[which(lengths(cellList) > 0),]
		cellCommVec <- cellCommVec[which(lengths(cellList) > 0)]
	}

	# prepare output object
	obj <- vector('list', length = 7)
	names(obj) <- c('grid', 'speciesList', 'cellCommInd', 'geogSpecies', 'cellCount', 'data', 'phylo')

	obj[['grid']] <- gridTemplate
	obj[['speciesList']] <- uniqueComm
	obj[['cellCommInd']] <- cellCommVec
	obj[['geogSpecies']] <- sort(unique(names(spGridList)))
	obj[['cellCount']] <- lengths(spGridList)
	attr(obj, 'resolution') <- resolution
	if (inherits(gridTemplate, 'sf')) {
		attr(obj, 'crs') <- sf::st_crs(gridTemplate)$input
	} else {
		attr(obj, 'crs') <- terra::crs(gridTemplate, proj4 = TRUE)
	}
	if (inherits(gridTemplate, 'sf')) {
		attr(obj, 'projected') <- !sf::st_is_longlat(gridTemplate)
	} else {
		attr(obj, 'projected') <- !sf::st_is_longlat(proj)
	}
	attr(obj, 'gridType') <- cellType
	attr(obj, 'metric') <- 'spRichness'
	
	class(obj) <- 'epmGrid'
	return(obj)	
}


# Function to take a polygon and return occupied grid cells
# subset allows us to operate on a subset of grid cells at a time, if necessary
polyToGridCells <- function(poly, method, gridTemplate, gridCentroids, coverCutoff, subset) {

	if (subset == 0) {
		clusterInd <- 1:nrow(gridTemplate)
	} else {
		clusterInd <- which(gridTemplate$clusters == subset)
	}
	
	if (method == 'centroid') {
		tmp <- clusterInd[unlist(sf::st_intersects(poly, gridCentroids[clusterInd, ]))]
	}

	if (method == 'areaCutoff') {
		tmp <- unlist(sf::st_intersects(poly, gridTemplate[clusterInd, ]))
		tmpWithin <- clusterInd[tmp][lengths(sf::st_within(gridTemplate[clusterInd[tmp],], poly)) > 0]
		tmpCrosses <- setdiff(clusterInd[tmp], tmpWithin)
		# tmpWithin <- tmp[lengths(sf::st_within(gridTemplate[tmp,], poly)) > 0]
		# tmpCrosses <- setdiff(tmp, tmpWithin)
				
		if (length(tmpWithin) > 0 & length(tmpCrosses) > 0) {
			tmpCrosses <- tmpCrosses[as.vector(sf::st_area(sf::st_intersection(sf::st_union(poly), gridTemplate[tmpCrosses,])) / sf::st_area(gridTemplate[tmpCrosses,])) >= coverCutoff]
			tmp <- union(tmpWithin, tmpCrosses)
		}
	}
	
	return(tmp)
}










## if retainSmallSpecies == TRUE, then the cell that has the greatest occupancy for a species (regardless of the amount) is marked as present. 

polyToHex <- function(poly, method, coverCutoff, extentVec, resolution, crs, nGroups, retainSmallRanges, nThreads) {
	
	# Generate template
	if (!inherits(extentVec, 'sf')) {
		masterExtent <-  sf::st_as_sfc(extentVec)
	} else {
		masterExtent <- extentVec
	}
	gridTemplate <- sf::st_make_grid(masterExtent, cellsize = c(resolution, resolution), square = FALSE, crs = sf::st_crs(masterExtent))
	
	# if extent was polygon, then mask the grid template
	if (inherits(masterExtent, c('sf', 'sfc'))) {
		gridTemplate <- gridTemplate[lengths(sf::st_intersects(gridTemplate, masterExtent)) > 0,]
	}	
	
	gridTemplate <- sf::st_sf(gridTemplate, grid_id = 1:length(gridTemplate))
	
	gridCentroids <- suppressWarnings(sf::st_centroid(gridTemplate))
	
	# do kmeans clustering on centroid coords.
	if (nGroups > 1) {
		clusters <- stats::kmeans(sf::st_coordinates(gridCentroids), centers = nGroups)
		gridTemplate$clusters <- clusters$cluster
	}
	
	# prep for parallel computing
	if (nThreads == 1) {
		cl <- NULL
	} else if (nThreads > 1) {
		if (.Platform$OS.type != 'windows') {
			cl <- nThreads
		} else {
			cl <- parallel::makeCluster(nThreads)
			parallel::clusterExport(cl, c('polyToGridCells', 'method', 'gridTemplate', 'gridCentroids', 'coverCutoff'))
		}
	}
	
	if (nGroups > 1) {
		subsetsList <- list()
		for (i in 1:nGroups) {
			
			message('\tProcessing subgroup ', i, ' of ', nGroups, '...')
			
			# get grid cells per species
			subsetsList[[i]] <- pbapply::pblapply(poly, function(x) polyToGridCells(x, method = method, gridTemplate, gridCentroids, coverCutoff, subset = i), cl = cl)
		}
		
		spGridList <- lapply(1:length(subsetsList[[1]]), function(i) {
			unlist(lapply(subsetsList, function(x) x[[i]]))
		})
		names(spGridList) <- names(subsetsList[[1]])
		rm(subsetsList, clusters)	
	
	
	} else {
		
	    if (unique(as.character(sf::st_geometry_type(poly[[1]]))) == 'MULTIPOLYGON') {
		    spGridList <- pbapply::pblapply(poly, function(x) polyToGridCells(x, method = method, gridTemplate, gridCentroids, coverCutoff, subset = 0), cl = cl)
		    # for (i in 1:length(poly)) {
		    #     message(i)
		    #     tmp <- epm:::polyToGridCells(poly[[i]], method = method, gridTemplate, gridCentroids, coverCutoff, subset = 0)
		    # }
	    } else {
	        # for points
	        spGridList <- pbapply::pblapply(poly, function(x) unlist(sf::st_intersects(x, gridTemplate)), cl = cl)
	    }
	}

	if (nThreads > 1 & .Platform$OS.type == 'windows') parallel::stopCluster(cl)
		
	# if we are retaining small ranged species that would otherwise be dropped,
	# then we will search for species that did not register in any grid cell
	smallSp <- which(lengths(spGridList) == 0)
	
	if (retainSmallRanges) {
		if (length(smallSp) > 0) {
			for (i in 1:length(smallSp)) {
			    # message('\t', i)

				tmp <- unlist(sf::st_intersects(poly[[smallSp[i]]], gridTemplate))
				
				if (length(tmp) > 0) {
				
  				    # of grid cells intersected by small-ranged species, keep the cell most occupied.
  				    # this is to avoid going from a species that would not appear in any cell, to a species occuring in multiple cells with 1% coverage.
  				    tmp <- tmp[which.max(as.vector(sf::st_area(sf::st_intersection(sf::st_union(poly[[smallSp[i]]]), gridTemplate[tmp,])) / sf::st_area(gridTemplate[tmp,])))]
				
  				    spGridList[[smallSp[i]]] <- tmp

  			    }
			}
			rescued <- setdiff(names(smallSp), names(which(lengths(spGridList) == 0)))
			if (length(rescued) > 0) {
				msg <- paste0(length(rescued), ' small-ranged species ', ifelse(length(rescued) > 1, 'were', 'was'), ' preserved:\n\t', paste(rescued, collapse='\n\t'))
				message(msg)
			}
		}
	}
	
	return(list(gridTemplate, spGridList))
}	
	








# square gridcells via the terra package
 
### if method == 'centroid', then species belongs to cell if its range intersects a cell's centroid.
### if method == 'areaCutoff', then species belongs to cell if it covers some specified percent of cell area.

## if retainSmallSpecies == TRUE, then the cell that has the greatest occupancy for a species (regardless of the amount) is marked as present. 


polyToTerra <- function(poly, method, coverCutoff, extentVec, resolution, crs, retainSmallRanges, template = NULL, nThreads) {

	# Generate template
	if (is.null(template)) {
		if (inherits(extentVec, 'sf')) {
			bb <- getExtentOfList(extentVec, format = 'sf')
		} else {
			bb <- extentVec
		}
		gridTemplate <- terra::rast(xmin = bb['xmin'], xmax = bb['xmax'], ymin = bb['ymin'], ymax = bb['ymax'], resolution = c(resolution, resolution), crs = crs$input)
	} else {
		gridTemplate <- terra::rast(template)
		bb <- extentVec
	}
	
	# if extent was polygon, then mask cells that fall outside
	if (inherits(extentVec, 'sf')) {
	
		xx <- terra::cells(gridTemplate, y = terra::vect(extentVec), weights = TRUE)
		if (terra::ncell(gridTemplate) != length(unique(xx[, 'cell']))) {
			goodCells <- unique(xx[, 'cell'])
			cellsToExclude <- TRUE
		} else {
			cellsToExclude <- FALSE
		}
	} else {
		cellsToExclude <- FALSE
	}
	
	# prep for parallel computing
	if (nThreads == 1) {
		cl <- NULL
	} else if (nThreads > 1) {
		if (.Platform$OS.type != 'windows') {
			cl <- nThreads
		} else {
			cl <- parallel::makeCluster(nThreads)
			parallel::clusterExport(cl, c('method', 'gridTemplate', 'coverCutoff'))
		}
	}


	spGridList <- pbapply::pblapply(poly, function(x) {
	# for (i in 1:length(poly)) {
		# x <- poly[[i]]
		# message('\t', i)
	    
	    if (all(unique(as.character(sf::st_geometry_type(x))) == 'POINT')) {
	        
	        presenceCells <- terra::cellFromXY(gridTemplate, sf::st_coordinates(x))
	        if (cellsToExclude) {
	            presenceCells <- intersect(presenceCells, goodCells)
	        }
	        
        } else {
		
    		xx <- terra::cells(gridTemplate, y = terra::vect(x), weights = TRUE)
    		
    		if (nrow(xx) > 0) {
    			# if polygon extends beyond grid template, cell may be NaN
    			if (anyNA(xx[, 'cell'])) {
    				xx <- xx[!is.na(xx[, 'cell']),]
    			}
    					
    			# cells are returned regardless of how much they are covered by polygon
    			
    			# exclude some cells if needed
    			if (cellsToExclude) {
    				xx <- xx[xx[, 'cell'] %in% goodCells, ]
    			}
    			
    			if (!is.matrix(xx)) {
    				# if xx was a single row, the subsetting might cause it to become a vector
    				xx <- matrix(xx, ncol = 3, dimnames = list(NULL, names(xx)))
    			}
    			
                if (nrow(xx) > 0) {
    			    
                    if (method == 'centroid') {
                        centroids <- terra::xyFromCell(gridTemplate, cell = xx[, 'cell'])
                        centroids <- sf::st_as_sf(as.data.frame(centroids), coords = 1:2, crs = sf::st_crs(x))	
                        presenceCells <- xx[, 'cell'][unlist(sf::st_intersects(x, centroids))]
        		
                    } else if (method == 'areaCutoff') {
        		
                        presenceCells <- xx[xx[, 'weights'] >= coverCutoff, 'cell']
                        # terra::xyFromCell(gridTemplate, cell = presenceCells) # for testing
        					
                    }
                } else {
                    presenceCells <- numeric(0)
                }
    		} else {
    		    presenceCells <- numeric(0)
    		}
        }
	    
	    presenceCells
	    
	}, cl = cl)
	
	if (nThreads > 1 & .Platform$OS.type == 'windows') parallel::stopCluster(cl)
	
	## for testing
	# plot(st_geometry(poly[[i]]))
	# testGrid <- rast(gridTemplate)
	# testGrid[xx[, 'cell']] <- 1
	# plot(as.polygons(testGrid, dissolve = F), col = 'orange')	
	# plot(st_geometry(poly[[i]]), add=T)


	# if we are retaining small ranged species that would otherwise be dropped,
	# then we will search for species that did not register in any grid cell
	smallSp <- which(lengths(spGridList) == 0)
	
	if (retainSmallRanges) {
		if (length(smallSp) > 0) {
			
		    gridTemplateHighRes <- terra::rast(terra::disaggregate(gridTemplate, fact = 10))
			
			for (i in 1:length(smallSp)) {
				
				xx <- terra::cells(gridTemplateHighRes, y = terra::vect(poly[[smallSp[i]]]), weights = TRUE)

				if (nrow(xx) > 0) {
									
					smallCells <- xx[which.max(xx[, 'weights']), 'cell']
					smallCells <- terra::cellFromXY(gridTemplate, terra::xyFromCell(gridTemplateHighRes, cell = smallCells))
					if (cellsToExclude) {
						smallCells <- intersect(smallCells, goodCells)
					}
					
					spGridList[[smallSp[i]]] <- smallCells
				}

				# ## alternative approach using point sampling
				# xx <- terra::cells(gridTemplateHighRes, y = terra::vect(poly[[smallSp[i]]]), weights = TRUE)
				# if (nrow(xx) > 0) {
					# pts <- terra::xyFromCell(gridTemplateHighRes, cell = xx[, 'cell'])
				# } else {
					# # if range is too small to register, sample points within the polygon
					# pts <- sf::st_sample(poly[[smallSp[i]]], size = 100)
					# pts <- sf::st_coordinates(pts)
				# }
				# smallCells <- unique(terra::cellFromXY(gridTemplate, pts))
				# spGridList[[smallSp[i]]] <- smallCells		
			}
			
			rescued <- setdiff(names(smallSp), names(which(lengths(spGridList) == 0)))
			if (length(rescued) > 0) {
				msg <- paste0(length(rescued), ' small-ranged species ', ifelse(length(rescued) > 1, 'were', 'was'), ' preserved:\n\t', paste(rescued, collapse='\n\t'))
				message(msg)
			}
		}
	}
	
	# # if extent was polygon, then mask cells that fall outside
	# if (inherits(extentVec, 'sf')) {
	
		# xx <- terra::cells(gridTemplate, y = terra::vect(extentVec), weights = TRUE)
		# if (terra::ncell(gridTemplate) != length(unique(xx[, 'cell']))) {
			# spGridList <- lapply(spGridList, function(x) intersect(x, unique(xx[, 'cell'])))
		# }
	# }
	
	return(list(gridTemplate, spGridList))
}

				
				
				
# internal function for preparing occurrence input for createEPMgrid().
##
## @param occ a table of species coordinates, a list of species-specific tables of coordinates, 
##		a spatial coordinates object, or a species-specific list of spatial coordinate objects (sf or sp).
##		If in table form, each coordinate pair must have an associated species name. If in list form, 
## 	each element of the list must be named with the name of the species.

## @examples
## library(raster)
## library(sf)
## # example dataset: a list of 24 chipmunk distributions as polygons
## # To illustrate usage, we will randomly sample some coordinates from each species polygon
## # and generate a couple of alternative input formats
## 
## # list of sf spatial objects
## spOccList <- lapply(tamiasPolyList, function(x) st_sample(x, size = 10, type= 'random'))
## spStack <- rasterStackFromOccurrences(spOccList, resolution = 50000, crs = st_crs(tamiasPolyList[[1]]))
##
## # list of coordinate tables
## spOccList2 <- lapply(spOccList, function(x) st_coordinates(x))
## spStack <- rasterStackFromOccurrences(spOccList2, resolution = 50000, crs = st_crs(tamiasPolyList[[1]]))
## 
## # single table of coordinates
## spOccList3 <- spOccList2
## for (i in 1:length(spOccList3)) {
## 	spOccList3[[i]] <- cbind.data.frame(taxon = names(spOccList3)[i], spOccList3[[i]])
## 	colnames(spOccList3[[i]]) <- c('taxon', 'X', 'Y')
## }
## spOccList3 <- do.call(rbind, spOccList3)
## rownames(spOccList3) <- NULL
## spOccList3[, 'taxon'] <- as.character(spOccList3[, 'taxon'])
## spStack <- rasterStackFromOccurrences(spOccList3, resolution = 50000, 
## 	coordHeaders = c('X', 'Y'), crs = st_crs(spOccList[[1]]))
##
## # a single labeled spatial object
## spOccList4 <- st_as_sf(spOccList3[, c('taxon', 'X', 'Y')], coords = c('X','Y'), 
## 	crs = st_crs(spOccList[[1]])$proj4string)
## spStack <- rasterStackFromOccurrences(spOccList4, resolution = 50000)

###########################################

## # list of sf spatial objects
## spOccList <- lapply(tamiasPolyList, function(x) st_sample(x, size = 10, type= 'random'))
## occurrenceFormatting(spOccList)
##
## # list of coordinate tables
## spOccList2 <- lapply(spOccList, function(x) st_coordinates(x))
## occurrenceFormatting(spOccList2)
## 
## # single table of coordinates
## spOccList3 <- spOccList2
## for (i in 1:length(spOccList3)) {
## 	spOccList3[[i]] <- cbind.data.frame(taxon = names(spOccList3)[i], spOccList3[[i]])
## 	colnames(spOccList3[[i]]) <- c('taxon', 'X', 'Y')
## }
## spOccList3 <- do.call(rbind, spOccList3)
## rownames(spOccList3) <- NULL
## spOccList3[, 'taxon'] <- as.character(spOccList3[, 'taxon'])
## occurrenceFormatting(spOccList3)
##
## # a single labeled spatial object
## spOccList4 <- st_as_sf(spOccList3[, c('taxon', 'X', 'Y')], coords = c('X','Y'), 
## 	crs = st_crs(spOccList[[1]])$proj4string)
## occurrenceFormatting(spOccList4)
##
## # unprojected list
## spOccList5 <- lapply(spOccList, function(x) st_transform(x, crs = 4326))
##
## # unprojected list of tables
## spOccList6 <- lapply(spOccList5, function(x) st_coordinates(x))


occurrenceFormatting <- function(occ) {
	
	# detect format and convert. Target is a list of separate tables of coordinates, one per species.
	if (inherits(occ, 'list')) {
		
		# do all elements in the list have the same class?
		if (length(unique(sapply(occ, class, simplify = FALSE))) != 1) {
			stop('Not all elements in occ have the same class.')
		}
		
		if (inherits(occ[[1]], c('SpatialPoints', 'SpatialPointsDataFrame', 'sf', 'sfc'))) {

			if (inherits(occ[[1]], c('SpatialPoints', 'SpatialPointsDataFrame'))) {
				occ <- lapply(occ, sf::st_as_sf)
			}
			
			crs <- sf::st_crs(occ[[1]])$proj4string
			
			if (length(unique(sapply(occ, function(x) sf::st_crs(x)$proj4string, simplify = FALSE))) != 1) {
				stop('Not all elements in occ have the same projection.')
			}
			
			occ <- lapply(occ, sf::st_geometry)

			if (inherits(occ[[1]], 'sfc_POINT')) {
				occ <- lapply(occ, sf::st_coordinates)
			}
	
			if (inherits(occ[[1]], 'sfc_MULTIPOINT')) {
				occ <- lapply(occ, function(x) sf::st_coordinates(x)[, 1:2])
			}	
		}
		
		for (i in 1:length(occ)) {
			occ[[i]] <- cbind.data.frame(taxon = names(occ)[i], occ[[i]])
			colnames(occ[[i]]) <- c('taxon', 'long', 'lat')
		}
	
		occ <- do.call(rbind, occ)
		rownames(occ) <- NULL
		occ[, 'taxon'] <- as.character(occ[, 'taxon'])			
			
	} else {
		
		if (inherits(occ, c('SpatialPoints', 'SpatialPointsDataFrame', 'sf', 'sfc'))) {
		
			if (inherits(occ, c('SpatialPoints', 'SpatialPointsDataFrame'))) {
				occ <- sf::st_as_sf(occ)
			}
			
			crs <- sf::st_crs(occ)$proj4string
			
			headers <- colnames(occ)
			headers <- setdiff(headers, c('long', 'lat'))
			occ <- as.data.frame(occ)		
			occ <- cbind.data.frame(occ[, setdiff(headers, 'geometry')], sf::st_coordinates(occ[, 'geometry']))
			colnames(occ) <- c(setdiff(headers, 'geometry'), c('long', 'lat'))
		} else if (!inherits(occ, c('matrix', 'data.frame'))) {
			stop('Input format of occ not recognized.')
		}
	}
	

	# Now, regardless of input format, occ should be a single table of coordinates and associated taxon names
	
	if (!all(c('long', 'lat') %in% colnames(occ))) {
	    if (all(c('x', 'y') %in% tolower(colnames(occ)))) {
	        colnames(occ)[which(tolower(colnames(occ)) == 'x')] <- 'long'
	        colnames(occ)[which(tolower(colnames(occ)) == 'y')] <- 'lat'
	    } else {
    		stop('Coordinate headers not found in occ.')
	    }
	}
	if (!'taxon' %in% colnames(occ)) {
		stop('taxon header not found in occ.')
	}
	
	occ <- as.data.frame(occ, stringsAsFactors = FALSE)
	occ[, 'long'] <- as.numeric(occ[, 'long'])
	occ[, 'lat'] <- as.numeric(occ[, 'lat'])

	occ <- occ[, c('taxon', 'long', 'lat')]
	occ[,'taxon'] <- gsub('^\\s+|\\s+$', '', occ[, 'taxon'])
	occ[,'taxon'] <- gsub('\\s+', '_', occ[, 'taxon'])
	
	# are any records missing critical information?
	if (any(occ[, 'taxon'] == '' | is.na(occ[, 'taxon']))) {
		ind <- which(occ[, 'taxon'] == '' | is.na(occ[, 'taxon']))
		message('\t', length(ind), ' records were dropped because taxon was not specified.')
		occ <- occ[ - ind, ]
	}
	
	if (any(occ[, 'long'] == '' | is.na(occ[, 'long']) | occ[, 'lat'] == '' | is.na(occ[, 'lat']))) {
		ind <- which(occ[, 'long'] == '' | is.na(occ[, 'long']) | occ[, 'lat'] == '' | is.na(occ[, 'lat']))
		message('\t', length(ind), ' records were dropped because they lacked coordinates.')
		occ <- occ[ - ind, ]
	}
	
	# split into list of tables
	occList <- split(occ, occ[, 'taxon'])
	
	return(occList)

}






