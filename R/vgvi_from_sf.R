#' Viewshed Greenness Visibility Index (VGVI) from sf
#' @description The VGVI expresses the proportion of visible greenness to the total visible area based on a \code{\link[GVI]{viewshed}}.
#' The estimated VGVI values range between 0 and 1, where 0 = no green cells are visible, and 1 = all of the visible cells are green.
#' A distance decay function is applied, to account for the reducing visual prominence of an object in space with increasing distance from the observer.
#'
#' @param sf_start object of class \code{sf}; See ‘Details’
#' @param dsm_data object of class \code{\link[terra]{SpatRaster}}; \code{\link[terra]{SpatRaster}} of the DSM
#' @param dtm_data object of class \code{\link[terra]{SpatRaster}}; \code{\link[terra]{SpatRaster}} of the DTM
#' @param greenspace object of class \code{\link[terra]{SpatRaster}}; Binary \code{\link[terra]{SpatRaster}} of the Greenspace mask. Values must be 1 for Green and 0 for No-Green
#' @param max_distance numeric; Buffer distance to calculate the viewshed
#' @param observer_height numeric > 0; Height of the observer (e.g. 1.7 meters)
#' @param resolution optional; NULL or numeric >= 1; Resolution that the GVI should be aggregated to
#' @param m numeric; See ‘Details’
#' @param b numeric; See ‘Details’
#' @param mode character; 'logit' or 'exponential'. See ‘Details’
#' @param cores numeric; The number of cores to use, i.e. at most how many child processes will be run simultaneously
#' @param chunk_size numeric; Chunk size for parallelization. See ‘Details’ 
#' @param progress logical; Show progress bar?
#'
#' @details 
#' sf_start needs to be a geometry of type POINT, LINESTRING, MULTILINESTRING, POLYGON or MULTIPOLYGON. If sf_start is a LINESTRING or MULTILINESTRING, 
#' points will be generated along the line(s) every "resolution" meters. If sf_start is a POLYGON or MULTIPOLYGON, a grid with resolution = "resolution" 
#' will be generated, and VGVI will be computed for every point.
#' The CRS (\code{\link[sf]{st_crs}}) needs to have a metric unit!
#' 
#' The type of function, used for calculating the distance decay weights, can be defined with the \code{mode} parameter.
#' The argument 'logit' uses the logistic function, d = 1 / (1 + e^(b * (x - m))) and 'exponential' the exponential function d = 1 / (1 + (b * x^m)).
#' The decay function can be visualized using the \code{\link[GVI]{visualizeWeights}} function.
#' 
#' Higher values of chunk_size increase computation time, but may also be more RAM intensive. Also, if you choose a high number of cores, RAM usage increases. too.
#' It is highly recommended to use a Linux or Mac system, for better parallel performance.
#'
#' @return sf_object containing the weighted VGVI values as POINT features, where 0 = no green cells are visible, and 1 = all of the visible cells are green. 
#' @export
#' 
#' @importFrom magrittr %>%
#' @importFrom sf st_crs
#' @importFrom sf st_geometry_type
#' @importFrom sf st_union
#' @importFrom sf st_cast
#' @importFrom sf st_line_sample
#' @importFrom sf st_as_sf
#' @importFrom sf st_bbox
#' @importFrom sf st_buffer
#' @importFrom sf st_coordinates
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom dplyr n
#' @importFrom terra crs
#' @importFrom terra res
#' @importFrom terra rast
#' @importFrom terra crop
#' @importFrom terra mask
#' @importFrom terra vect
#' @importFrom terra xyFromCell
#' @importFrom terra terraOptions
#' @importFrom terra aggregate
#' @importFrom terra extract
#' @importFrom terra cellFromXY
#' @importFrom raster raster
#' @importFrom raster as.raster
#' @importFrom methods is
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%

vgvi_from_sf <- function(sf_start, dsm_data, dtm_data, greenspace,
                         max_distance = 800, observer_height = 1.7, 
                         resolution = NULL, m = 0.5, b = 8,
                         mode = c("logit", "exponential"), 
                         cores = 1, chunk_size = 999, progress = TRUE) {
  
  #### 1. Check input ####
  # sf_start
  valid_sf_types <- c("POINT", "LINESTRING", "MULTILINESTRING", "POLYGON", "MULTIPOLYGON")
  if (!is(sf_start, "sf")) {
    stop("sf_start must be a sf object")
  } else if (sf::st_crs(sf_start)$units != "m") {
    stop("sf_start CRS unit needs to be metric")
  } else if (!as.character(sf::st_geometry_type(sf_start, by_geometry = FALSE)) %in% valid_sf_types) {
    stop("sf_start has no valid geometry")
  }
  
  # dsm_data
  if (!is(dsm_data, "SpatRaster")) {
    stop("dsm_data needs to be a SpatRaster object!")
  } else if (sf::st_crs(terra::crs(dsm_data))$epsg != sf::st_crs(sf_start)$epsg) {
    stop("dsm_data needs to have the same CRS as sf_start")
  }
  
  # dtm_data
  if (!is(dsm_data, "SpatRaster")) {
    stop("dtm_data needs to be a SpatRaster object!")
  } else if (sf::st_crs(terra::crs(dtm_data))$epsg != sf::st_crs(sf_start)$epsg) {
    stop("dtm_data needs to have the same CRS as sf_start")
  }
  
  # greenspace
  if (!is(dsm_data, "SpatRaster")) {
    stop("greenspace needs to be a SpatRaster object!")
  } else if (sf::st_crs(terra::crs(greenspace))$epsg != sf::st_crs(sf_start)$epsg) {
    stop("greenspace needs to have the same CRS as sf_start")
  }
  
  # resolution
  if (is.null(resolution)) {
    resolution = min(terra::res(dsm_data))
  } else if (resolution < 1 || resolution < min(terra::res(dsm_data))) {
    stop("If the resolution differs from dsm_data, it needs to be > 1 and higher than the dsm_data resolution")
  }
  
  #### 2. Convert sf_start to points ####
  if (as.character(sf::st_geometry_type(sf_start, by_geometry = FALSE)) %in% c("LINESTRING", "MULTILINESTRING")) {
    sf_start <- sf_start %>%
      sf::st_union() %>%
      sf::st_cast("LINESTRING") %>%
      sf::st_line_sample(density = 1/resolution) %>%
      sf::st_cast("POINT") %>%
      sf::st_as_sf() %>% 
      dplyr::rename(geom = x)
  } else if (as.character(sf::st_geometry_type(sf_start, by_geometry = FALSE)) %in% c("POLYGON", "MULTIPOLYGON")) {
    sf_start_bbox <- sf::st_bbox(sf_start)
    sf_start <- terra::rast(xmin = sf_start_bbox[1], xmax = sf_start_bbox[3], 
                            ymin = sf_start_bbox[2], ymax = sf_start_bbox[4], 
                            crs = terra::crs(dsm_data), resolution = resolution, vals = 0) %>% 
      terra::crop(sf_start) %>% 
      terra::mask(terra::vect(sf_start)) %>%
      terra::xyFromCell(which(terra::values(.) == 0)) %>%
      as.data.frame() %>% 
      sf::st_as_sf(coords = c("x","y"), crs = sf::st_crs(sf_start)) %>% 
      dplyr::rename(geom = geometry)
  }
  
  sf_start <- sf_start %>% 
    dplyr::mutate(VGVI = as.numeric(NA)) %>% 
    dplyr::select(VGVI)
  
  #### 3. Prepare data for viewshed analysis ####
  # Max AOI
  max_aoi <- sf_start %>% 
    sf::st_buffer(max_distance) %>% 
    dplyr::mutate(id = 1:dplyr::n())
  
  # Crop DSM to max AOI and change resolution
  dsm_data <- terra::crop(dsm_data, terra::vect(max_aoi))
  
  if(resolution != min(raster::res(dsm_data))) {
    terra::terraOptions(progress = 0)
    dsm_data <- terra::aggregate(dsm_data, fact = resolution/terra::res(dsm_data))
    terra::terraOptions(progress = 3)
  }
  
  # Coordinates of start point
  x0 <- sf::st_coordinates(sf_start)[,1]
  y0 <- sf::st_coordinates(sf_start)[,2]
  
  # Observer heights
  height_0_vec <- unlist(terra::extract(dtm_data, cbind(x0, y0)), use.names = F) + observer_height
  
  # Convert to Raster
  dsm_data <- raster::raster(dsm_data)
  
  #### 4. Calculate viewsheds and VGVI ####
  max_aoi_list <- suppressWarnings(split(max_aoi, seq(1, nrow(max_aoi), chunk_size)))
  
  if (progress) {
    pb = txtProgressBar(min = 0, max = length(max_aoi_list), initial = 0, style = 3)
  }
  if (cores > 1 && Sys.info()[["sysname"]] == "Windows") {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  }
  for (j in seq_along(max_aoi_list)) {
    
    this_aoi <- max_aoi_list[[j]]
    this_ids <- this_aoi$id
    this_x0 <- x0[this_ids]
    this_y0 <- y0[this_ids]
    this_height_0_vec <- height_0_vec[this_ids]
    
    #### Calculate viewshed
    if (cores > 1) {
      if (Sys.info()[["sysname"]] == "Windows") {
        par_fun <-  function(i){
          viewshed_raster(this_aoi = this_aoi[i,], dsm_data = dsm_data,
                          x0 = this_x0[i], y0 = this_y0[i], height0 = this_height_0_vec[i],
                          resolution = resolution, id = this_ids[i])
        }
        
        viewshed_list <- foreach::foreach(i=1:nrow(this_aoi)) %dopar% par_fun(i)
      }
      else {
        viewshed_list <- parallel::mclapply(1:nrow(this_aoi), function(i){
          viewshed_raster(this_aoi = this_aoi[i,], dsm_data = dsm_data,
                          x0 = this_x0[i], y0 = this_y0[i], height0 = this_height_0_vec[i],
                          resolution = resolution, id = this_ids[i])},
          mc.cores = cores, mc.preschedule = TRUE)
      }
    } else {
      viewshed_list <- lapply(1:nrow(this_aoi), function(i){
        viewshed_raster(this_aoi = this_aoi[i,], dsm_data = dsm_data,
                        x0 = this_x0[i], y0 = this_y0[i], height0 = this_height_0_vec[i],
                        resolution = resolution, id = this_ids[i])})
    }
    
    # Get id's from raster object (in case the order has been changed during parallelization)
    this_ids <- lapply(viewshed_list, names) %>% 
      unlist() %>% 
      gsub("X", "", .) %>% 
      as.integer()
    
    # Convert to terra rast
    viewshed_list <- lapply(viewshed_list, terra::rast)
    
    #### Calculate VGVI
    # Calculate XY-visible (XYV) table
    this_XYV_table_list <- lapply(viewshed_list, function(x){
      # Get XY coordinates of visible cells
      xy <- x %>% 
        terra::xyFromCell(which(x[] == 1))
      
      # Calculate distance
      centroid <- colMeans(terra::xyFromCell(x, which(!is.na(x[]))))
      dxy = round(sqrt((centroid[1] - xy[,1])^2 + (centroid[2] - xy[,2])^2))
      dxy[dxy==0] = min(dxy[dxy!=0])
      
      # Intersect XY with greenspace mask
      output <- greenspace[terra::cellFromXY(greenspace, xy)] %>% 
        unlist(use.names = FALSE) %>% 
        cbind(dxy, .)
      colnames(output) <- c("dxy", "visible")
      return(output)
    })
    
    # Compute VGVI
    if (cores > 1) {
      if (Sys.info()[["sysname"]] == "Windows") {
        par_fun <-  function(i){
          vgvi_from_XYV_table(XYV_table = i, m = m, b = b, mode = mode)
        }
        
        this_vgvis <- foreach::foreach(i=this_XYV_table_list) %dopar% par_fun(i)
      }
      else {
        this_vgvis <- parallel::mclapply(this_XYV_table_list, vgvi_from_XYV_table,
                                         m = m, b = b, mode = mode,
                                         mc.cores = cores, mc.preschedule = TRUE)
      }
    } else {
      this_vgvis <- lapply(this_XYV_table_list, vgvi_from_XYV_table,
                           m = m, b = b, mode = mode)
    }
    
    # Update sf_start
    sf_start[this_ids,1] <- unlist(this_vgvis)
    
    # Update ProgressBar
    if (progress) setTxtProgressBar(pb, j)
  }
  if (cores > 1 && Sys.info()[["sysname"]] == "Windows") {
    parallel::stopCluster(cl)
  }
  
  return(sf_start)
}
