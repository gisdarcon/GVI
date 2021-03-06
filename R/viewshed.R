#' @title Viewshed
#' @description Computes the viewshed of a point on a Digital Surface Model map.
#'
#' @param observer object of class \code{sf} with one point; Starting location
#' @param max_distance numeric; Buffer distance to calculate the viewshed
#' @param dsm_rast object of class \code{\link[terra]{SpatRaster}}; \code{\link[terra]{SpatRaster}} of the DSM
#' @param dtm_rast object of class \code{\link[terra]{SpatRaster}}; \code{\link[terra]{SpatRaster}} of the DTM
#' @param observer_height numeric > 0; Height of the observer (e.g. 1.7 meters)
#' @param resolution optional; NULL or numeric > 0; Resolution that the GVI raster should be aggregated to. Needs to be a multible of the dsm_rast resolution
#' @param plot optional; Plot the intersect of the buffer around the observer location and the DSM (left DSM; right visibility raster)
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
#' @importFrom terra boundaries
#' @importFrom terra xyFromCell
#' @importFrom terra plot
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics par
#' @importFrom graphics points
#' @useDynLib GVI, .registration = TRUE

viewshed <- function(observer, dsm_rast, dtm_rast, 
                     max_distance = 800, observer_height = 1.7, 
                     resolution = NULL, plot = FALSE) {
  #### 1. Check input ####
  # observer
  if (!is(observer, "sf")) {
    stop("observer must be a sf object")
  } else if (sf::st_crs(observer)$units != "m") {
    stop("observer CRS unit needs to be metric")
  } else if (!as.character(sf::st_geometry_type(observer, by_geometry = FALSE)) == "POINT") {
    stop("observer has no valid geometry")
  } else if (nrow(observer) > 1) {
    stop("The viewshed function can only calculate the viewshed of one point. Please look into the vgvi_from_sf function")
  }
  
  # dsm_rast
  if (!is(dsm_rast, "SpatRaster")) {
    stop("dsm_rast needs to be a SpatRaster object")
  } else if (sf::st_crs(terra::crs(dsm_rast))$epsg != sf::st_crs(observer)$epsg) {
    stop("dsm_rast needs to have the same CRS as observer")
  }
  
  # dtm_rast
  if (!is(dtm_rast, "SpatRaster")) {
    stop("dtm_rast needs to be a SpatRaster object")
  } else if (sf::st_crs(terra::crs(dtm_rast))$epsg != sf::st_crs(observer)$epsg) {
    stop("dtm_rast needs to have the same CRS as observer")
  }
  
  # resolution
  dsm_res <- min(terra::res(dsm_rast))
  if (is.null(resolution)) {
    resolution = dsm_res
  } else if (resolution < min(terra::res(dsm_rast))) {
    stop("resolution must be higher than the resolution of dsm_rast")
  } else if ((resolution %% dsm_res) != 0) {
    stop(paste0("resolution must be a multible of the dsm_rast resolution. Try resolution = ", resolution - (resolution %% dsm_res)))
  }
  rm(dsm_res)
  
  #### 2. Prepare Data for viewshed analysis ####
  # AOI
  this_aoi <- observer %>% 
    sf::st_buffer(max_distance)
  
  # Coordinates of start point
  x0 <- sf::st_coordinates(observer)[1]
  y0 <- sf::st_coordinates(observer)[2]
  
  # Observer height
  height0 <- as.numeric(terra::extract(dtm_rast, cbind(x0, y0))) + observer_height
  
  # If the resolution parameter differs from the input-DSM resolution,
  # resample the DSM to the lower resolution.
  if (resolution == min(terra::res(dsm_rast))) {
    dsm_rast_masked <- terra::crop(dsm_rast, terra::vect(this_aoi)) %>% 
      terra::mask(terra::vect(this_aoi))
    
    output <- terra::setValues(dsm_rast_masked, 0) %>% 
      terra::mask(dsm_rast_masked)
  } else {
    terra::terraOptions(progress = 0)
    dsm_rast_masked <- terra::crop(dsm_rast, terra::vect(this_aoi)) %>% 
      terra::aggregate(fact = resolution/terra::res(.)) %>% 
      terra::mask(terra::vect(this_aoi))
    terra::terraOptions(progress = 3)
    
    output <- terra::setValues(dsm_rast_masked, 0) %>%
      terra::mask(dsm_rast_masked)
  }
  
  #### 3. Compute viewshed ####
  # Start row/col
  r0 <- terra::rowFromY(output, y0)
  c0 <- terra::colFromX(output, x0)
  
  # Convert output raster to matrix
  dsm_mat <- terra::values(dsm_rast_masked) %>% 
    matrix(ncol = terra::ncol(dsm_rast_masked))
  
  # Calculate boundaries of output raster (boundaries are adjacent to NA values)
  output_boundaries <- terra::extend(output, resolution*2) %>% 
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
    terra::plot(dsm_rast_masked, legend = F); graphics::points(x0, y0, col = "red", pch = 20, cex = 2)
    terra::plot(output, legend = F); graphics::points(x0, y0, col = "red", pch = 20, cex = 2)
    graphics::par(mfrow=c(1,1))
  }
  return(output)
}