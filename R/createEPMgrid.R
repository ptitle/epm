##'@title Create epmGrid object
##'
##'@description Creates an epmGrid object from a range of species-specific
##'  inputs.
##'
##'
##'@param spDat a number of possible input formats are possible. See details
##'  below.
##'
##'@param resolution vertical and horizontal spacing of grid cells, in units of
##'  the polygons' or points' projection.
##'
##'@param method approach used for gridding. Either \code{centroid} or
##'  \code{percentOverlap}. See details below.
##'
##'@param cellType either \code{hexagon} or \code{square}. See details below.
##'
##'@param percentThreshold the percent that a species range must cover a grid
##'  cell to be considered present. Specified as a proportion.
##'
##'@param retainSmallRanges boolean; should small ranged species be dropped or
##'  preserved. See details.
##'
##'@param extent if 'auto', then the maximal extent of the polygons will be
##'  used. If not 'auto', can be a SpatialPolygon, sf object, or raster, in
##'  which case the resulting epmGrid will be cropped and masked with respect to
##'  the polygon; or a spatial coordinates object, from which an extent object
##'  will be generated; or a numeric vector of length 4 with minLong, maxLong,
##'  minLat, maxLat. See \code{\link{interactiveExtent}} to draw your own
##'  extent.
##'
##'@param percentWithin The percentage of a species range that must be within
##'  the defined extent in order for that species to be included. This filter
##'  can be used to exclude species whose range barely enters the area of
##'  interest. The default value of 0 will disable this filter. If \code{extent
##'  == 'auto'}, then this filter will also have no effect, as the extent is
##'  defined by the species' ranges.
##'
##'@param checkValidity if \code{TRUE}, then check polygon validity and repair
##'  if needed, using sf::st_make_valid.
##'
##'@param crs if supplying occurrence records in a non-spatial format, then you
##'  must specify the crs. For unprojected long/lat data, you can simply provide
##'  \code{crs = 4326}.
##'
##'@param nThreads if > 1, then employ parallel computing. This won't
##'  necessarily improve runtime.
##'
##'@param template an object of class \code{SpatRaster} or \code{RasterLayer}
##'  that can be used to get extent and resolution. If \code{cellType =
##'  'square'}, then the template will be used as the reference grid.
##'
##'@param verbose if TRUE, list out all species that are dropped/excluded,
##'  rather than counts.
##'
##'@param use.data.table if \code{'auto'}, this is determined by the size of the
##'  dataset. Primarily intended for debugging.
##'
##'
##'@details Types of accepted inputs for argument \code{spDat}: 
##'\enumerate{
##'\item a list of polygon objects (sf or sp), named with taxon names. 
##'\item a list of SpatRaster or RasterLayer grids, named with taxon names. \item a
##'multi-layer RasterStack or multi-layer SpatRaster. 
##'\item a set of occurrence records, multiple accepted formats, see below. 
##'\item a site-by-taxon presence/absence matrix. 
##'}
##'
##'If input data consist of \strong{occurrence records} rather than polygons,
##'then a couple of formats are possible: 
##'\enumerate{ 
##'\item You can provide a list of species-specific spatial point objects. 
##'\item You can provide a single spatial object, where points have a taxon attribute. 
##'\item You can provide a list of non-spatial species-specific dataframes. 
##'\item You can provide a single non-spatial dataframe. 
##'}
##'
##'For options (1) and (3), the taxon names must be provided as the list names.
##'For options (2) and (4), the columns must be 'taxon', 'x' and 'y' (or 'long',
##''lat'). For options (3) and (4), as these are non-spatial, you must provide a
##'crs object to the \code{crs} argument, so that the function knows what
##'projection to use.
##'
##'It is also possible to supply a \strong{matrix with sites as rows and taxa as
##'columns}. The contents of this matrix must be either 0 or 1. If this is the
##'case, then a raster grid must be supplied under the template argument. This
##'will be the grid system used for converting this presence/absence matrix to
##'an epmGrid object. It is expected that the index order of the grid is the
##'same as the row order of the matrix.
##'
##'If input is a set of \strong{species-specific grids}, then it is expected
##'that all grids belong to the same overall grid system, i.e. that the cells
##'align and that all grids have the same resolution. Grids do not need to have
##'the same extent.
##'
##'Any SpatialPolygon or SpatialPoints objects are converted to objects of class
##'\code{sf}.
##'
##'
##'If \code{cellType = 'hexagon'}, then the grid is made of polygons via the sf
##'package. If \code{cellType = 'square'}, then the grid is a raster generated
##'via the terra package. Hexagonal cells have several advantages, including
##'being able to be of different sizes (if the grid is in unprojected long/lat),
##'and may be able to more naturally follow coastlines and non-linear features.
##'However, the raster-based square cells will be much less memory intensive for
##'high resolution datasets. Choice of grid type matters more for spatial
##'resolution (total number of cells), than for number of species.
##'
##'
##'In the polygon-to-grid conversion process, two approaches are implemented.
##'For \code{method = 'centroid'}, a range polygon registers in a cell if the
##'polygon overlaps with the cell centroid. For \code{method =
##''percentOverlap'}, a range polygon registers in a cell if it covers that cell
##'by at least \code{percentThreshold} fraction of the cell.
##'
##'If \code{retainSmallRanges = FALSE}, then species whose ranges are so small
##'that no cell registers as present will be dropped. If \code{retainSmallRanges
##'= TRUE}, then the cell that contains the majority of the the small polygon
##'will be considered as present, even if it's a small percent of the cell.
##'
##'If \code{retainSmallRanges = TRUE}, and an extent is provided, then species
##'may still be dropped if they fall outside of that extent.
##'
##'You may see the message \code{Failed to compute min/max, no valid pixels found in 
##'sampling. (GDAL error 1)} . This just means that a species did not register in any 
##'grid cells. If you specified \code{retainSmallRanges = TRUE}, then those species will
##'be included in a subsequent step. Therefore, this message can be ignored.
##'
##'For very large datasets, this function will make a determination as to
##'whether or not there is sufficient memory. If there is not, an alternative
##'approach that uses the data.table package will be employed. Please install
##'this R package to take advantage of this feature.
##'
##'@return an object of class \code{epmGrid}.
##'
##'@author Pascal Title
##'
##' @examples
##' library(sf)
##' # example dataset: a list of 24 chipmunk distributions as polygons
##' head(tamiasPolyList)
##'
##' # hexagonal grid
##' \donttest{
##' tamiasEPM <- createEPMgrid(tamiasPolyList, resolution = 50000,
##' 	cellType = 'hexagon', method = 'centroid')
##' tamiasEPM
##' }
##' # square grid
##' tamiasEPM2 <- createEPMgrid(tamiasPolyList, resolution = 50000,
##' 	cellType = 'square', method = 'centroid')
##' tamiasEPM2
##'
##' #######
##' \donttest{
##' # demonstration of site-by-species matrix as input.
##' tamiasEPM2 <- createEPMgrid(tamiasPolyList, resolution = 50000,
##' 	cellType = 'square', method = 'centroid')
##'
##' ## first we will use the function generateOccurrenceMatrix() to get
##' ## a presence/absence matrix
##' pamat <- generateOccurrenceMatrix(tamiasEPM2, sites = 'all')
##'
##' # here, our grid template will be tamiasEPM2[[1]]
##' tamiasEPM2[[1]]
##' xx <- createEPMgrid(pamat, template = tamiasEPM2[[1]])
##'
##'
##' #######
##' # demonstration with grids as inputs
##' ## We will first generate grids from the range polygons
##' ## (you normally would not do this -- you would have grids from some other source)
##'
##' # define the extent that contains all range polygons
##' fullExtent <- terra::ext(terra::vect(tamiasPolyList[[1]]))
##' for (i in 2:length(tamiasPolyList)) {
##' 	fullExtent <- terra::union(fullExtent, terra::ext(terra::vect(tamiasPolyList[[i]])))
##' }
##'
##' # create raster template
##' fullGrid <- terra::rast(fullExtent, res = 50000, crs = terra::crs(terra::vect(tamiasPolyList[[1]])))
##'
##' # now we can convert polygons to a common grid system
##' spGrids <- list()
##' for (i in 1:length(tamiasPolyList)) {
##' 	spGrids[[i]] <- terra::trim(terra::rasterize(terra::vect(tamiasPolyList[[i]]), fullGrid))
##' }
##' names(spGrids) <- names(tamiasPolyList)
##'
##' createEPMgrid(spGrids)
##'
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
##'}
##'
##'
##'@export

createEPMgrid <- function(spDat, resolution = 50000, method = 'centroid', cellType = 'hexagon', percentThreshold = 0.25, retainSmallRanges = TRUE, extent = 'auto', percentWithin = 0, checkValidity = FALSE, crs = NULL, nThreads = 1, template = NULL, verbose = FALSE, use.data.table = 'auto') {
	
	# spDat <- tamiasPolyList; resolution = 50000; method = 'centroid'; cellType = 'hexagon'; percentThreshold = 0.1; retainSmallRanges = TRUE; extent = 'auto'; percentWithin = 0; checkValidity = FALSE; nThreads = 1; template = NULL; verbose = TRUE; use.data.table = 'auto';
	
	# test with occurrences
	# spOccList <- lapply(tamiasPolyList, function(x) sf::st_sample(x, size = 10, type= 'random'))
	# spDat <- spOccList; resolution = 50000; method = 'centroid'; cellType = 'hexagon'; percentThreshold = 0.1; retainSmallRanges = TRUE; extent = 'auto'; checkValidity = FALSE; nThreads = 1; template = NULL; verbose = TRUE
	
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
		if (all(c('taxon', 'x', 'y') %in% colnames(spDat)) | all(c('taxon', 'long', 'lat') %in% colnames(spDat))) {
			datType <- 'points'
			if (!is.null(crs)) {
				proj <- crs
			} else {
				stop('Input occurrence data need a crs but none provided.')
			}
		} else {
			# input must be a site by taxon table
			datType <- 'siteMat'
			if (is.null(template)) {
				stop("Input recognized as site x taxon table, but no grid template provided.")
			}
			
			if (inherits(template, 'RasterLayer')) {
				template <- as(template, 'SpatRaster')
			}
			
			if (!inherits(template, 'SpatRaster')) {
				stop('At this time, only a rectangular grid can be provided for site x taxon matrices.')
			}
		
			if (terra::ncell(template) != nrow(spDat)) {
				stop("Grid template has different number of cells than there are rows in matrix.")
			}
			
			proj <- terra::crs(template)
			cellType <- 'square'
			retainSmallRanges <- FALSE
			percentWithin <- 0
			checkValidity <- FALSE					
		}
	
	# list of non-spatial tables
	} else if (!inherits(spDat, c('data.frame', 'matrix')) & inherits(spDat[[1]], c('data.frame', 'matrix')) & !inherits(spDat[[1]], c('sf', 'sfc'))) {
		datType <- 'points'
		if (!is.null(crs)) {
			proj <- crs
		} else {
			stop('Input occurrence data need a crs but none provided.')
		}
	
	} else if (inherits(spDat, c('SpatRaster', 'RasterStack')) | inherits(spDat[[1]], c('SpatRaster', 'RasterLayer'))) {
		
	# input is list or stack or collection of raster grids
	
		# if spDat is a stack, then we already have rasters with the same extent and resolution
		if (inherits(spDat, 'SpatRaster')) {
			spDatStack <- spDat
		} else if (is.list(spDat)) {
			for (i in 1:length(spDat)) {
				if (inherits(spDat[[i]], 'RasterLayer')) {
					spDat[[i]] <- as(spDat[[i]], 'SpatRaster')
				}
			}
			if (!all(sapply(spDat, class) == 'SpatRaster')) {
				stop('Some of the items in spDat are not recognized as raster grids.')
			}
		
			if (!all(sapply(spDat, function(x) terra::compareGeom(spDat[[1]], x, ext = FALSE, rowcol = FALSE)))) {
				stop('All raster grids are expected to be part of the same grid system.')
			}
			
			fullExtent <- terra::ext(spDat[[1]])
			for (i in 2:length(spDat)) {
				fullExtent <- terra::union(fullExtent, terra::ext(spDat[[i]]))
			}
			
			for (i in 1:length(spDat)) {
				spDat[[i]] <- terra::extend(spDat[[i]], fullExtent)
			}
		
			spDatStack <- terra::rast(spDat)
		}
		
		# this can be converted to a cell by taxon table with accompanying grid template
		spDat <- as.matrix(spDatStack)
		spDat[is.nan(spDat)] <- 0
		template <- terra::rast(spDatStack)
		datType <- 'siteMat'
		cellType <- 'square'
		retainSmallRanges <- FALSE
		percentWithin <- 0
		checkValidity <- FALSE	
		proj <- terra::crs(template)				
	
	} else {
		stop('Format of spDat not recognized.')
	}
		
	if (datType == 'points') {
		percentWithin <- 0
		method <- 'centroid'
		# extent <- 'auto'
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
	method <- match.arg(method, c('centroid', 'percentOverlap'))
	
	if (!cellType %in% c('hexagon', 'square')) {
		stop('cellType must be either hexagon or square.')
	}

	if (!method %in% c('centroid', 'percentOverlap')) {
		stop('method must be either centroid or percentOverlap.')
	}

	if (is.list(spDat) & is.null(names(spDat))) {
		stop('List must be named with species names.')
	}
	
	if (is.list(spDat) & anyDuplicated(names(spDat))) {
		stop('Some taxon names are duplicated.')
	}
	
	if (percentThreshold > 1) {
		stop('percentThreshold is a fraction.')
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
	
	# disable terra progress bar, as it interferes with epm progress bar
	terra::terraOptions(progress = 0)
	
	if (!is.null(template)) {
		if (!inherits(template, c('SpatRaster', 'RasterLayer'))) {
			stop('template must be a RasterLayer or SpatRaster.')
		}
		if (inherits(template, 'RasterLayer')) {
			template <- as(template, 'SpatRaster')
		}
		if (datType != 'siteMat') {
			if ((sum(c(terra::is.lonlat(template, perhaps = TRUE), sf::st_is_longlat(spDat[[1]]))) == 1)) {
				stop('Raster provided as template has a different projection from input data.')
			}
		}
		resolution <- terra::res(template)[1]
		extent <- as.vector(terra::ext(template))
		if (cellType == 'hexagon') {
			stop('Use of the template argument is intended for square-cell grids only.')
		}
	}
		
	if (checkValidity) {
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
	if (datType != 'siteMat') {
		if (!any(sapply(spDat, function(x) identical(sf::st_crs(spDat[[1]]), sf::st_crs(x))))) {
			stop('proj4string or EPSG of all polygons must match.')
		}

		if (sf::st_is_longlat(spDat[[1]])) {
			proj <- sf::st_crs(4326)
		} else {
			proj <- sf::st_crs(spDat[[1]])
		}
	}

	# if data are unprojected, then catch likely mistake where resolution is specified in meters.
	## (data would be in decimal degrees, and resolution would be in the 100s or 1000s)
	# also, if data are projected, and resolution is < 90, then likely a mistake as well.
	if (sf::st_is_longlat(proj) & resolution > 90) {
		stop("Input data are in longlat, therefore resolution must be in decimal degrees.")
	}
	if (!sf::st_is_longlat(proj) & resolution < 90) {
		stop('Input data are projected, but resolution does not seem to be in the same units.')
	}
	
	# if WKT string, then convert to sf polygon
	if (inherits(extent, 'character')) {
		if (grepl('POLYGON \\(', extent)) {
			extent <- sf::st_as_sfc(extent, crs = sf::st_crs(spDat[[1]]))
		}
	}

	if (inherits(extent, 'character')) {
		
		if (all(extent != 'auto')) {
			stop("If extent is a character vector, it must be 'auto'.")
		}
		
		if (all(extent == 'auto')) {
		
			#get overall extent
			masterExtent <- getExtentOfList(spDat)
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
	
	
	# if input is a site-by-taxon table, then process this now and exit
	if (datType == 'siteMat') {
		return(processSiteBySpeciesMatrix(spDat, template))
	}
	
	# here, we take one of two routes:
	## if hexagonal grid, then use sf, if square grid, use terra
	
	uniqueSp <- sort(names(spDat))
	spDat <- spDat[uniqueSp]
	nGroups <- 1 # placeholder
	
	if (cellType == 'hexagon') {
		
		# poly = spDat; method = method; percentThreshold = percentThreshold; extentVec = masterExtent; resolution = resolution; crs = proj; nGroups = nGroups; retainSmallRanges = retainSmallRanges; nThreads = nThreads; verbose = TRUE
		spGridList <- polyToHex(poly = spDat, method = method, percentThreshold = percentThreshold, extentVec = masterExtent, resolution = resolution, crs = proj, nGroups = nGroups, retainSmallRanges = retainSmallRanges, nThreads = nThreads, verbose = verbose)
		
	} else if (cellType == 'square') {
		
			spGridList <- polyToTerra(poly = spDat, method = method, percentThreshold = percentThreshold, extentVec = masterExtent, resolution = resolution, crs = proj, retainSmallRanges = retainSmallRanges, template = template, nThreads = nThreads, verbose = verbose)
		
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
	
	if (inherits(gridTemplate, 'SpatRaster')) {
		nGridCells <- terra::ncell(gridTemplate)
	} else if (inherits(gridTemplate, 'sf')) {
		nGridCells <- nrow(gridTemplate)
	} else {
		stop('gridTemplate class mismatch.')
	}
	
	if (inherits(try(matrix(0, nrow = nGridCells, ncol = length(spGridList))), 'try-error') | use.data.table == TRUE) {
		# too much memory required. Use data.table approach
		
		if (!requireNamespace('data.table', quietly = TRUE)) {
			stop('Not enough memory available. Please install the data.table package to switch to an alternative approach.')
		}

		mat <- data.table::data.table('cell' = integer(), 'sp' = character(), key = 'cell')
		for (i in 1:length(spGridList)) {
			mat <- rbind(mat, data.table::data.table('cell' = spGridList[[i]], 'sp' = names(spGridList)[i]))
		}
	
		# if square grid system, we can't drop empty cells, so add them
		if (inherits(gridTemplate, 'SpatRaster')) {
			mat <- rbind(mat, data.table::data.table('cell' = setdiff(1:terra::ncell(gridTemplate), unique(mat$cell)), 'sp' = NA))
		}
			
		data.table::setkey(mat, 'cell')
		
		# get set of taxa for each cell
		sp <- NULL; cell <- NULL
		cellComms <- mat[, toString(sort(sp)), by = cell]
			
		# get unique communities across cells
		uniqueComm <- unique(cellComms$V1)
		
		# create vector of indices that map unique communities to all communities
		cellCommVec <- match(cellComms$V1, uniqueComm)
		
		# convert uniqueComm to a list of species
		uniqueComm <- strsplit(uniqueComm, ',\\s+')
		uniqueComm <- lapply(uniqueComm, sort)
				
		# remove empty cells if hexagonal grid
		if (inherits(gridTemplate, 'sf')) {
			gridTemplate <- gridTemplate[sort(unique(mat$cell)),]
		}
	
	} else {
		
		# there is enough available memory to just create a cell by sp matrix
				
		# create site by species matrix
		mat <- matrix(0, nrow = nGridCells, ncol = length(spGridList))
		colnames(mat) <- names(spGridList)
		for (i in 1:length(spGridList)) {
			mat[spGridList[[i]], names(spGridList)[i]] <- 1
		}
		
		# create condensed version that encodes species at each site
		# This creates a character vector, concatenating each row in mat with -
		if (requireNamespace('data.table', quietly = TRUE) & (use.data.table == TRUE | use.data.table == 'auto')) {
			matCondensed <- fpaste(data.table::as.data.table(mat), "-")$V1
		} else {
			matCondensed <- do.call(paste, c(as.data.frame(mat), sep="-"))
		}
		uniqueComm <- unique(matCondensed)
		
		# create vector of indices that map unique communities to all communities
		cellCommVec <- match(matCondensed, uniqueComm)
		
		# convert unique community codes to species names
		uniqueComm <- strsplit(uniqueComm, split = '-')
		uniqueComm <- lapply(uniqueComm, as.integer)
		uniqueComm <- lapply(uniqueComm, function(x) sort(names(spGridList)[as.logical(x)]))

		# remove empty cells if hexagonal grid
		if (inherits(gridTemplate, 'sf')) {
			gridTemplate <- gridTemplate[which(lengths(uniqueComm[cellCommVec]) > 0),]
			cellCommVec <- cellCommVec[which(lengths(uniqueComm[cellCommVec]) > 0)]
		}
		
	}

	# add unique community index and species richness in as attributes
	if (inherits(gridTemplate, 'SpatRaster')) {
		gridTemplate <- terra::rast(gridTemplate, nlyrs = 2)
		terra::values(gridTemplate) <- cbind(cellCommVec, lengths(uniqueComm[cellCommVec]))
		names(gridTemplate) <- c('uniqueComm', 'spRichness')
		gridTemplate$spRichness[gridTemplate$spRichness == 0] <- NA
	} else if (inherits(gridTemplate, 'sf')) {
		gridTemplate$uniqueComm <- cellCommVec
		gridTemplate$spRichness <- lengths(uniqueComm[cellCommVec])
	} else {
		stop('gridTemplate class mismatch.')
	}
	
	uniqueComm[lengths(uniqueComm) == 0] <- NA
	
	# set back to default value
	terra::terraOptions(progress = 3)

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
		attr(obj, 'crs') <- terra::crs(gridTemplate, proj = TRUE)
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

.datatable.aware <- TRUE

# courtesy of https://stackoverflow.com/questions/14568662/paste-multiple-columns-together
fpaste <- function(dt, sep = ",") {
	x <- tempfile()
	data.table::fwrite(dt, file = x, sep = sep, col.names = FALSE, showProgress = FALSE)
	data.table::fread(x, sep = "\n", header = FALSE, showProgress = FALSE)
}

## alternate version
polyToHex <- function(poly, method, percentThreshold, extentVec, resolution, crs, nGroups, retainSmallRanges, nThreads, verbose) {
	
	# Combine species polygons, keep only geometry, and return single sf object
	taxonNames <- names(poly)
	poly <- lapply(poly, function(x) sf::st_combine(sf::st_geometry(x)))
	poly <- do.call(c, poly)
	
	# Generate template
	if (!inherits(extentVec, 'sf')) {
		masterExtent <-  sf::st_as_sfc(extentVec)
	} else {
		masterExtent <- extentVec
	}
	gridTemplate <- sf::st_make_grid(masterExtent, cellsize = c(resolution, resolution), square = FALSE, crs = sf::st_crs(masterExtent))
	
	# if extent was polygon, then mask the grid template
	if (inherits(masterExtent, c('sf', 'sfc'))) {
		gridTemplate <- gridTemplate[lengths(sf::st_intersects(gridTemplate, masterExtent)) > 0]
	}	
	
	gridTemplate <- sf::st_sf(gridTemplate, grid_id = 1:length(gridTemplate))
	
	gridCentroids <- sf::st_centroid(sf::st_geometry(gridTemplate))
			
	if (unique(as.character(sf::st_geometry_type(poly[[1]]))) == 'MULTIPOLYGON') {
	
		# keep polygons if there are geometry collections
		if (!all(as.character(sf::st_geometry_type(poly)) %in% c('MULTIPOLYGON', 'POLYGON'))) {
			poly <- sf::st_collection_extract(poly, type = 'POLYGON')
		}
			
		if (method == 'centroid') {
			if (verbose) message('\t\t Converting ranges to grid...')
				# spGridList <- sf::st_intersects(poly, gridCentroids, sparse = FALSE)
				# spGridList <- apply(spGridList, 1, function(x) which(x == TRUE))
				spGridList <- sf::st_intersects(poly, gridCentroids, sparse = TRUE)
		}		
		
		if (method == 'percentOverlap') {
			if (verbose) message('\t\t Converting ranges to grid...')
			tmp <- sf::st_intersects(poly, gridTemplate, sparse = TRUE)
			
			tmpCrosses <- sf::st_intersects(sf::st_cast(poly, 'MULTILINESTRING'), gridTemplate, sparse = TRUE)
		
			# plot(poly[i])
			# plot(gridTemplate[tmpWithin, ], add = T, col = 'blue')
			# plot(gridTemplate[tmpCrosses[[i]], ], add = T, col = 'orange')
			# plot(poly[i], add = TRUE)
			
			# check for and try to resolve geometry issues
			for (i in 1:length(poly)) {
				if (!any(sf::st_is_valid(poly[i]))) {
					poly[i] <- sf::st_make_valid(poly[i])
				}
			}

			spGridList <- vector('list', length(tmp))
			pb <- txtProgressBar(max = length(tmp), style = 3)
			for (i in 1:length(tmp)) {
				# message('\t', i)
				setTxtProgressBar(pb, i)
				if (length(tmp[[i]]) > 0) {
					tmpWithin <- setdiff(tmp[[i]], tmpCrosses[[i]])
					xx <- tmpCrosses[[i]][as.vector(sf::st_area(sf::st_intersection(sf::st_union(poly[i]), gridTemplate[tmpCrosses[[i]],])) / sf::st_area(gridTemplate[tmpCrosses[[i]],])) >= percentThreshold]
					## experimental
					# xx <- tmpCrosses[[i]][geos::geos_area(geos::geos_intersection(geos::as_geos_geometry(sf::st_union(poly[i])), geos::as_geos_geometry(gridTemplate[tmpCrosses[[i]],]))) / geos::geos_area(geos::as_geos_geometry(gridTemplate[tmpCrosses[[i]],])) >= percentThreshold]
					spGridList[[i]] <- union(tmpWithin, xx)
				}
			}
			if (verbose) cat('\n'); close(pb)
		}
	} else {
		# for points
		spGridList <- sf::st_intersects(poly, gridTemplate, sparse = TRUE)
	}

	names(spGridList) <- taxonNames

	# if we are retaining small ranged species that would otherwise be dropped,
	# then we will search for species that did not register in any grid cell
	smallSp <- which(lengths(spGridList) == 0)
	
	if (retainSmallRanges) {
		if (length(smallSp) > 0) {
			if (verbose) message('\t\t Recovering small-ranged taxa...')
			pb <- txtProgressBar(max = length(smallSp), style = 3)
			for (i in 1:length(smallSp)) {
				# if (verbose) message('\t\ttaxon ', i)
				setTxtProgressBar(pb, i)

				spGridList[[smallSp[i]]] <- findTopCells5(poly[smallSp[i]], gridTemplate, resolution)

			}
			if (verbose) cat('\n'); close(pb)
			rescued <- setdiff(names(smallSp), names(which(lengths(spGridList) == 0)))
			if (length(rescued) > 0) {
				msg <- paste0(length(rescued), ' small-ranged species ', ifelse(length(rescued) > 1, 'were', 'was'), ' preserved:\n\t', paste(rescued, collapse='\n\t'))
				message(msg)
			}
		}
	}
	
	return(list(gridTemplate, spGridList))
		
}

findTopCells <- function(x, grid) {
	
	tmp <- unlist(sf::st_intersects(x, grid))
			
	if (length(tmp) > 0) {
			
		# of grid cells intersected by small-ranged species, keep the cell most occupied.
		# this is to avoid going from a species that would not appear in any cell, to a species occuring in multiple cells with 1% coverage.
	  					
		check <- sf::st_intersects(sf::st_cast(x, 'POLYGON'), sf::st_cast(x, 'POLYGON'))
		if (any(lengths(check)) > 1) {
			mergedPoly <- sf::st_union(x)
			mergedPoly <-sf::st_cast(mergedPoly, 'POLYGON')
		} else {
			mergedPoly <- sf::st_cast(x)
		}			
		# now mergedPoly should be a set of separate polygons that do not share borders
		topCells <- integer(length(mergedPoly))
		for (j in 1:length(mergedPoly)) {
			areas <- sapply(tmp, function(y) {
				x1 <- as.vector(sf::st_area(sf::st_intersection(mergedPoly[j], grid[y,])))
				if (length(x1) == 0) x1 <- 0
				x2 <- as.vector(sf::st_area(grid[y,]))
				x1 / x2
			})
			topCells[j] <- tmp[which.max(areas)]	
		}  					
		sort(unique(topCells))
	} else {
		numeric(0)
	}
}

findTopCells2 <- function(x, grid) {
	
	tmp <- unlist(sf::st_intersects(x, grid))
			
	if (length(tmp) > 0) {
			
		# of grid cells intersected by small-ranged species, keep the cell most occupied.
		# this is to avoid going from a species that would not appear in any cell, to a species occuring in multiple cells with 1% coverage.
	  					
		check <- sf::st_intersects(sf::st_cast(x, 'POLYGON'), sf::st_cast(x, 'POLYGON'))
		if (any(lengths(check)) > 1) {
			mergedPoly <- sf::st_union(x)
			mergedPoly <-sf::st_cast(mergedPoly, 'POLYGON')
		} else {
			mergedPoly <- sf::st_cast(x)
		}			
		# now mergedPoly should be a set of separate polygons that do not share borders

		# rather than calculate intersection of range with each grid cell to determine percent coverage, 
		## let's install sample 200 regularly spaced points and determine how many are intersected by range.
		samplePts <- lapply(tmp, function(x) sf::st_sample(grid[x,], size = 200, type = 'regular'))
		
		topCells <- integer(length(mergedPoly))
		for (j in 1:length(mergedPoly)) {		
			ints <- lapply(samplePts, function(x) sf::st_intersects(mergedPoly, x))
			percs <- lengths(unlist(ints, recursive = FALSE)) / sapply(samplePts, length)
			topCells[j] <- tmp[which.max(percs)]	
		}		
		sort(unique(topCells))
	} else {
		numeric(0)
	}
}

findTopCells3 <- function(x, grid) {
	
	tmp <- unlist(sf::st_intersects(x, grid))
			
	if (length(tmp) > 0) {
			
		# of grid cells intersected by small-ranged species, keep the cell most occupied.
		# this is to avoid going from a species that would not appear in any cell, to a species occuring in multiple cells with 1% coverage.
	  					
		check <- sf::st_intersects(sf::st_cast(x, 'POLYGON'), sf::st_cast(x, 'POLYGON'))
		if (any(lengths(check)) > 1) {
			mergedPoly <- sf::st_union(x)
			mergedPoly <- sf::st_cast(mergedPoly, 'POLYGON')
		} else {
			mergedPoly <- sf::st_cast(x)
		}			
		# now mergedPoly should be a set of separate polygons that do not share borders

		# rather than calculate intersection of range with each grid cell to determine percent coverage, 
		## let's instead sample 200 regularly spaced points and determine how many are intersected by range.
		samplePts <- sf::st_sample(grid[tmp,], size = rep(200, length(tmp)), type = 'random', by_polygon = F)
		ptGrid <- sf::st_intersects(grid[tmp,], samplePts)
		ptTest <- sf::st_intersects(mergedPoly, samplePts)
		topCells <- integer(length(mergedPoly))
		for (j in 1:length(mergedPoly)) {
			topCells[j] <- tmp[which.max(lengths(sapply(ptGrid, function(y) intersect(ptTest[[j]], y))) / lengths(ptGrid))]
		}
		
		sort(unique(topCells))

	} else {
		numeric(0)
	}
}

findTopCells4 <- function(x, grid, resolution) {

	tmp <- unlist(sf::st_intersects(x, grid))
			
	if (length(tmp) > 0) {
			
		# of grid cells intersected by small-ranged species, keep the cell most occupied.
		# this is to avoid going from a species that would not appear in any cell, to a species occuring in multiple cells with 1% coverage.
	  					
		check <- sf::st_intersects(sf::st_cast(x, 'POLYGON'), sf::st_cast(x, 'POLYGON'))
		if (any(lengths(check)) > 1) {
			mergedPoly <- sf::st_union(x)
			mergedPoly <- sf::st_cast(mergedPoly, 'POLYGON')
		} else {
			mergedPoly <- sf::st_cast(x, 'POLYGON')
		}			
		# now mergedPoly should be a set of separate polygons that do not share borders

		# rather than calculate intersection of range with each grid cell to determine percent coverage, 
		## let's instead sample 200 regularly spaced points and determine how many are intersected by range.
		# ras <- terra::rast(extent = terra::ext(terra::vect(grid[tmp,])), resolution = rep(length(tmp) * 200, 2), crs = sf::st_crs(grid)$input)
		# ras <- terra::rasterize(terra::vect(grid[tmp,]), ras)
		# samplePts <- terra::crds(ras, df = TRUE)
		# samplePts <- sf::st_as_sf(samplePts, coords = 1:2, crs = sf::st_crs(grid)$input)


		ras <- terra::rast(extent = terra::ext(terra::vect(grid[tmp,])), resolution = rep(resolution/100, 2), crs = sf::st_crs(grid)$input)
		# ras <- terra::rasterize(terra::vect(grid[tmp,]), ras)
		polyCells <- terra::cells(ras, terra::vect(mergedPoly))
		gridCells <- terra::cells(ras, terra::vect(grid[tmp,]))
		topCells <- integer(length(mergedPoly))
		for (j in 1:length(mergedPoly)) {
			# topCells[j] <- tmp[which.max(tableAllCases(gridCells[gridCells[,2] %in% polyCells[polyCells[,1] == j, 2],1], gridCells[,1])[names(table(gridCells[,1]))] / table(gridCells[,1]))]
			xx <- table(gridCells[gridCells[,2] %in% polyCells[polyCells[,1] == j, 2], 1])
			topCells[j] <- names(xx)[which.max(xx / table(gridCells[,1])[names(xx)])]
		}
		topCells <- as.numeric(sort(unique(topCells)))
		
		# converting between tmp and terra raster ID's
		tmp[topCells]
	# cbind(tmp, terra::cells(ras, terra::vect(sf::st_centroid(sf::st_geometry(grid[tmp,])))))
		
		
		# ptGrid <- sf::st_intersects(grid[tmp,], samplePts)
		# ptTest <- sf::st_intersects(mergedPoly, samplePts)
		# topCells <- integer(length(mergedPoly))
		# for (j in 1:length(mergedPoly)) {
		# 	topCells[j] <- tmp[which.max(lengths(sapply(ptGrid, function(y) intersect(ptTest[[j]], y))) / lengths(ptGrid))]
		# }
		# 
		# sort(unique(topCells))
		# 
	} else {
		numeric(0)
	}
}	



findTopCells5 <- function(x, grid, resolution) {

	tmp <- unlist(sf::st_intersects(x, grid))
			
	if (length(tmp) > 0) {
			
		# of grid cells intersected by small-ranged species, keep the cell most occupied.
		# this is to avoid going from a species that would not appear in any cell, to a species occuring in multiple cells with 1% coverage.
	  					
		ras <- terra::rast(extent = terra::ext(terra::vect(grid[tmp,])), resolution = rep(resolution/100, 2), crs = sf::st_crs(grid)$input)
		gridcellras <- terra::rasterize(terra::vect(grid[tmp,]), ras, field = 'grid_id')
		polyRas <- terra::rasterize(terra::vect(x), ras, cover = TRUE)
		clusters <- terra::patches(polyRas, directions = 8, allowGaps = FALSE)

		# and identify most occupied cell per patch
		topCells <- integer(max(terra::minmax(clusters)))
		for (j in 1:max(terra::minmax(clusters))) {
			yy <- terra::mask(polyRas, clusters, maskvalues = j, inverse = TRUE)
			yy <- terra::mask(gridcellras, yy)
			topCells[j] <- as.integer(names(which.max(table(terra::values(yy)))))
		}
		sort(unique(topCells))
	} else {
		numeric(0)
	}
}			

# microbenchmark::microbenchmark(test1 = findTopCells(poly[smallSp[i]], gridTemplate), test2 = findTopCells2(poly[smallSp[i]], gridTemplate), test3 = findTopCells3(poly[smallSp[i]], gridTemplate), test4 = findTopCells4(poly[smallSp[i]], gridTemplate, resolution), test5 = findTopCells5(poly[smallSp[i]], gridTemplate, resolution), times = 10)


# square gridcells via the terra package
 
### if method == 'centroid', then species belongs to cell if its range intersects a cell's centroid.
### if method == 'percentOverlap', then species belongs to cell if it covers some specified percent of cell area.

## if retainSmallSpecies == TRUE, then the cell that has the greatest occupancy for a species (regardless of the amount) is marked as present. 


polyToTerra <- function(poly, method, percentThreshold, extentVec, resolution, crs, retainSmallRanges, template = NULL, nThreads, verbose = verbose) {

	# Generate template
	if (is.null(template)) {
		if (inherits(extentVec, 'sf')) {
			bb <- getExtentOfList(extentVec)
		} else {
			bb <- extentVec
		}
		# add half a cell's width to each dimension to be sure we are encompassing all points or polygons
		gridTemplate <- terra::rast(xmin = bb['xmin'] - resolution/2, xmax = bb['xmax'] + resolution/2, ymin = bb['ymin'] - resolution/2, ymax = bb['ymax'] + resolution/2, resolution = c(resolution, resolution), crs = crs$input)
	} else {
		gridTemplate <- terra::rast(template)
		bb <- extentVec
	}
	
	# if extent was polygon, then mask cells that fall outside
	if (inherits(extentVec, 'sf')) {
	
		gridext <- terra::vect(extentVec)
		gridmask <- terra::rasterize(terra::vect(extentVec), gridTemplate)
		if (anyNA(terra::values(gridmask))) {
			cellsToExclude <- TRUE
		} else {
			cellsToExclude <- FALSE
		}
	} else {
		gridext <- terra::ext(gridTemplate)
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
			parallel::clusterExport(cl, c('method', 'gridTemplate', 'percentThreshold'))
		}
	}


	spGridList <- pbapply::pblapply(poly, function(x) {
	# for (i in 1:length(poly)) {
	 	# x <- poly[[i]]
	 	# message('\t', i)
		
		if (all(unique(as.character(sf::st_geometry_type(x))) == 'POINT')) {
			
			presenceCells <- terra::cellFromXY(gridTemplate, sf::st_coordinates(x))
			if (cellsToExclude) {
				presenceCells <- presenceCells[!is.na(unlist(gridmask[presenceCells])), ]
				# presenceCells <- intersect(presenceCells, goodCells)
			}
			
		} else {
			
				# do the extents overlap? If not, then skip
				if (terra::relate(gridext, terra::ext(terra::vect(x)), relation = 'intersects')[1,1]) {
			
				if (method == 'centroid') {
					
					# rasterize polygon against grid (cells register if midpoint is within polygon)
					xx <- terra::rasterize(terra::vect(x), gridTemplate)
					
					# exclude some cells if needed
					if (cellsToExclude) {
						xx <- terra::mask(xx, gridmask)
					}
					
					presenceCells <- which(!is.na(terra::values(xx)))
					
				} else if (method == 'percentOverlap') {	
					
					# rasterize polygon against grid (returns cover fraction per cell)
					xx <- terra::rasterize(terra::vect(x), gridTemplate, cover = TRUE)
					
					xx[xx < percentThreshold] <- NaN
					
					# exclude some cells if needed
					if (cellsToExclude) {
						xx <- terra::mask(xx, gridmask)
					}
					
					presenceCells <- which(!is.na(terra::values(xx)))
				}
			} else {
				presenceCells <- integer(0)
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
			
			if (verbose) message('\t\t Recovering small-ranged taxa...')
			
			pb <- txtProgressBar(max = length(smallSp), style = 3)
			for (i in 1:length(smallSp)) {
				
				setTxtProgressBar(pb, i)
				
				xx <- terra::rasterize(terra::vect(poly[[smallSp[i]]]), gridTemplate, cover = TRUE)

				# exclude some cells if needed
				if (cellsToExclude) {
					xx <- terra::mask(xx, gridmask)
				}
				
				if (!all(is.na(terra::minmax(xx)))) {
				
					# to handle potentially discontiguous ranges, identify cell clusters
					clusters <- terra::patches(xx, directions = 8, allowGaps = FALSE)
					
					# and identify most occupied cell per patch
					topCells <- integer(max(terra::minmax(clusters)))
					for (j in 1:max(terra::minmax(clusters))) {
						yy <- terra::mask(xx, clusters, maskvalues = j, inverse = TRUE)
						topCells[j] <- which.max(terra::values(yy))
					}
					
					spGridList[[smallSp[i]]] <- sort(unique(topCells))
				}
			}
			
			if (verbose) cat('\n'); close(pb)
			
			rescued <- setdiff(names(smallSp), names(which(lengths(spGridList) == 0)))
			if (length(rescued) > 0) {
				msg <- paste0(length(rescued), ' small-ranged species ', ifelse(length(rescued) > 1, 'were', 'was'), ' preserved:\n\t', paste(rescued, collapse='\n\t'))
				message(msg)
			}
		}
	}
		
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


# internal function to take a site by species matrix and a grid template, and return a epmGrid object
processSiteBySpeciesMatrix <- function(mat, gridTemplate) {
		
	taxonNames <- colnames(mat)
	
	# force to be 0's or 1's
	mat[is.na(mat)] <- 0
	mat[mat != 0] <- 1
	
	# create condensed version that encodes species at each site
	if (requireNamespace('data.table', quietly = TRUE)) {
		matCondensed <- fpaste(data.table::as.data.table(mat), "-")$V1
	} else {
		matCondensed <- do.call(paste, c(as.data.frame(mat), sep="-"))
	}
	uniqueComm <- unique(matCondensed)
	
	# create vector of indices that map unique communities to all communities
	cellCommVec <- match(matCondensed, uniqueComm)
	
	# convert unique community codes to species names
	uniqueComm <- strsplit(uniqueComm, split = '-')
	uniqueComm <- lapply(uniqueComm, as.integer)
	uniqueComm <- lapply(uniqueComm, function(x) sort(taxonNames[as.logical(x)]))

	# add unique community index and species richness in as attributes
	if (inherits(gridTemplate, 'SpatRaster')) {
		gridTemplate <- terra::rast(gridTemplate, nlyrs = 2)
		terra::values(gridTemplate) <- cbind(cellCommVec, lengths(uniqueComm[cellCommVec]))
		names(gridTemplate) <- c('uniqueComm', 'spRichness')
		gridTemplate$spRichness[gridTemplate$spRichness == 0] <- NA
	} else {
		stop('gridTemplate can currently only be a SpatRaster.')
	}
	
	uniqueComm[lengths(uniqueComm) == 0] <- NA

	# prepare output object
	obj <- vector('list', length = 7)
	names(obj) <- c('grid', 'speciesList', 'cellCommInd', 'geogSpecies', 'cellCount', 'data', 'phylo')

	obj[['grid']] <- gridTemplate
	obj[['speciesList']] <- uniqueComm
	obj[['cellCommInd']] <- cellCommVec
	obj[['geogSpecies']] <- sort(taxonNames)
	obj[['cellCount']] <- colSums(mat)
	attr(obj, 'resolution') <- terra::res(gridTemplate)[1]
	if (inherits(gridTemplate, 'sf')) {
		attr(obj, 'crs') <- sf::st_crs(gridTemplate)$input
	} else {
		attr(obj, 'crs') <- terra::crs(gridTemplate, proj = TRUE)
	}
	if (inherits(gridTemplate, 'sf')) {
		attr(obj, 'projected') <- !sf::st_is_longlat(gridTemplate)
	} else {
		attr(obj, 'projected') <- !terra::is.lonlat(gridTemplate)
	}
	attr(obj, 'gridType') <- 'square'
	attr(obj, 'metric') <- 'spRichness'
	
	class(obj) <- 'epmGrid'
	return(obj)
}


