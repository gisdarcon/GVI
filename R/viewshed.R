#' @title Viewshed
#' @description Computes the viewshed of a point on a Digital Surface Model map.
#'
#' @param sf_start object of class \code{sf} with one point; Starting location
#' @param max_distance numeric; Buffer distance to calculate the viewshed 
#' @param dsm_data object of class \code{\link[terra]{SpatRaster}}; \code{\link[terra]{SpatRaster}} of the DSM
#' @param dtm_data object of class \code{\link[terra]{SpatRaster}}; \code{\link[terra]{SpatRaster}} of the DTM
#' @param observer_height numeric > 0; Height of the observer (e.g. 1.7 meters)
#' @param resolution optional; numeric > 1; Resolution that the GVI should be aggregated to. 
#' @param plot optional; Plot DSM and GVI
#'
#' @return object of class \code{\link[terra]{SpatRaster}}
#' @export
#' 
#' @importFrom magrittr %>%
#' @importFrom sf st_buffer
#' @importFrom sf st_coordinates
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
viewshed <- function(sf_start, dsm_data, dtm_data, 
                     max_distance = 800, observer_height = 1.7, 
                     resolution = NULL, plot = FALSE) {
  
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
  # Also, convert dsm_data_masked to "Raster" object, for faster internal calculation.
  if (is.null(resolution)) {
    resolution = min(terra::res(dsm_data))
  }
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
               terra::colFromX(output, xy1[,1])) %>% 
    split(seq(nrow(.)))
  
  # Apply lineOfSight function on every point in rc1
  this_LoS <- lapply(rc1, GVI::LoS, 
                     r0 = r0, c0 = c0, dsm_mat = dsm_mat, observerHeight = height0)
  
  # Bind list
  this_LoS <- unlist(this_LoS)
  
  # Copy result of lapply to the output raster
  #output[this_LoS[,1]] <- this_LoS[,2]
  output[this_LoS] <- 1
  
  # Compare DSM with Visibilty
  if (plot) {
    par(mfrow=c(1,2))
    terra::plot(dsm_data_masked, legend = F); points(x0, y0, col = "red", pch = 20, cex = 2)
    terra::plot(output, legend = F); points(x0, y0, col = "red", pch = 20, cex = 2)
    par(mfrow=c(1,1))
  }
  return(output)
}