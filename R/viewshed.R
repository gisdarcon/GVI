#' @title Viewshed
#' @description Computes the viewshed of a point on a Digital Surface Model map.
#'
#' @param sf_start object of class \code{sf} with one point; Starting location
#' @param max_distance numeric; Buffer distance to calculate the viewshed
#' @param dsm_data object of class \code{\link[terra]{SpatRaster}}; \code{\link[terra]{SpatRaster}} of the DSM
#' @param dtm_data object of class \code{\link[terra]{SpatRaster}}; \code{\link[terra]{SpatRaster}} of the DTM
#' @param observer_height numeric > 0; Height of the observer (e.g. 1.7 meters)
#' @param resolution optional; numeric >= 1; Resolution that the GVI should be aggregated to
#' @param plot optional; Plot DSM and GVI
#'
#' @return object of class \code{\link[terra]{SpatRaster}}
#' @export
#' 
#' @importFrom magrittr %>%
#' @importFrom sf st_buffer
#' @importFrom sf st_coordinates
#' @importFrom sf st_crs
#' @importFrom sf st_geometry_type
#' @importFrom terra crs
#' @importFrom terra extract
#' @importFrom terra res
#' @importFrom terra crop
#' @importFrom terra mask
#' @importFrom terra vect
#' @importFrom terra aggregate
#' @importFrom terra rowFromY
#' @importFrom terra colFromX
#' @importFrom terra values
#' @importFrom terra ncol
#' @importFrom terra expand
#' @importFrom terra boundaries
#' @importFrom terra xyFromCell
#' @importFrom terra plot
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics par
#' @importFrom graphics points
#' @useDynLib GVI, .registration = TRUE

viewshed <- function(sf_start, dsm_data, dtm_data, 
                     max_distance = 800, observer_height = 1.7, 
                     resolution = NULL, plot = FALSE) {
  #### 1. Check input ####
  # sf_start
  valid_sf_types <- c()
  if (!is(sf_start, "sf")) {
    stop("sf_start must be a sf object")
  } else if (sf::st_crs(sf_start)$units != "m") {
    stop("sf_start CRS unit needs to be metric")
  } else if (!as.character(sf::st_geometry_type(sf_start, by_geometry = FALSE)) == "POINT") {
    stop("sf_start has no valid geometry")
  } else if (nrow(sf_start) > 1) {
    stop("The viewshed function can only calculate the viewshed of one point. Please look into the vgvi_from_sf function.")
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
  
  # resolution
  if (is.null(resolution)) {
    resolution = min(terra::res(dsm_data))
  } else if (resolution < 1 || resolution < min(terra::res(dsm_data))) {
    stop("If the resolution differs from dsm_data, it needs to be > 1 and higher than the dsm_data resolution")
  }
  
  #### 2. Prepare Data for viewshed analysis ####
  # AOI
  this_aoi <- sf_start %>% 
    sf::st_buffer(max_distance)
  
  # Coordinates of start point
  x0 <- sf::st_coordinates(sf_start)[1]
  y0 <- sf::st_coordinates(sf_start)[2]
  
  # Observer height
  height0 <- as.numeric(terra::extract(dtm_data, cbind(x0, y0))) + observer_height
  
  # If the resolution parameter differs from the input-DSM resolution,
  # resample the DSM to the lower resolution.
  if (resolution == min(terra::res(dsm_data))) {
    dsm_data_masked <- terra::crop(dsm_data, this_aoi) %>% 
      terra::mask(terra::vect(this_aoi))
    
    output <- terra::setValues(dsm_data_masked, 0) %>% 
      terra::mask(dsm_data_masked)
  } else {
    terra::terraOptions(progress = 0)
    dsm_data_masked <- terra::crop(dsm_data, this_aoi) %>% 
      terra::aggregate(fact = resolution/terra::res(.)) %>% 
      terra::mask(terra::vect(this_aoi))
    terra::terraOptions(progress = 3)
    
    output <- terra::setValues(dsm_data_masked, 0) %>%
      terra::mask(dsm_data_masked)
  }
  
  #### 3. Compute viewshed ####
  # Start row/col
  r0 <- terra::rowFromY(output, y0)
  c0 <- terra::colFromX(output, x0)
  
  # Convert output raster to matrix
  dsm_mat <- terra::values(dsm_data_masked) %>% 
    matrix(ncol = terra::ncol(dsm_data_masked))
  
  # Calculate boundaries of output raster (boundaries are adjacent to NA values)
  output_boundaries <- terra::expand(output, resolution*2) %>% 
    terra::boundaries()
  
  # Get rows and columns of boundaries cells and convert to list
  xy1 <- terra::xyFromCell(output_boundaries, which(terra::values(output_boundaries) == 1))
  
  # Row/Col of boundary cells
  rc1 <- cbind(terra::rowFromY(output, xy1[,2]), 
               terra::colFromX(output, xy1[,1]))
  
  # Apply lineOfSight (C++) function on every point in rc1
  this_LoS <- LoSCpp(rc1 = rc1, r0 = r0, c0 = c0, observerHeight = height0, dsm_mat = dsm_mat)
  
  # Copy result of lineOfSight to the output raster
  output[this_LoS] <- 1
  
  #### 4. Compare DSM with Visibility ####
  if (plot) {
    graphics::par(mfrow=c(1,2))
    terra::plot(dsm_data_masked, legend = F); graphics::points(x0, y0, col = "red", pch = 20, cex = 2)
    terra::plot(output, legend = F); graphics::points(x0, y0, col = "red", pch = 20, cex = 2)
    graphics::par(mfrow=c(1,1))
  }
  return(output)
}