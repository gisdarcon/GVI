#' @importFrom magrittr %>%
#' @importFrom raster crop
#' @importFrom raster mask
#' @importFrom raster xyFromCell
#' @importFrom raster values
#' @importFrom raster extent
#' @importFrom raster crs
#' @importFrom raster raster
#' @importFrom raster res
#' @importFrom raster cellFromXY
#' @importFrom raster rowFromY
#' @importFrom raster colFromX
#' @importFrom raster ncol
#' @importFrom raster boundaries

viewshed_raster <- function(this_aoi, dsm_data, x0, y0,
                     height0, resolution, id) {
  
  # Crop DSM to AOI
  dsm_data_masked <- raster::crop(dsm_data, this_aoi) %>% 
    raster::mask(this_aoi)
  
  # Make empty raster with slightly bigger extent (resolution*2) and fill with 0 value
  valid_xy <- raster::xyFromCell(dsm_data_masked, which(!is.na(raster::values(dsm_data_masked))))
  
  out_ext <- raster::extent(dsm_data_masked)
  out_ext <- raster::extent(c(out_ext[1]-resolution*2, out_ext[2]+resolution*2,
                              out_ext[3]-resolution*2, out_ext[4]+resolution*2))
  
  output <- raster::raster(ext = out_ext, crs = raster::crs(dsm_data_masked),
                           resolution = raster::res(dsm_data_masked), vals = NA)
  output[raster::cellFromXY(output, valid_xy)] <- 0
  names(output) <- id
  
  # Start row/col
  r0 <- raster::rowFromY(dsm_data_masked, y0)
  c0 <- raster::colFromX(dsm_data_masked, x0)
  
  # Convert DSM raster to matrix
  dsm_mat <- raster::values(dsm_data_masked) %>% 
    matrix(ncol = raster::ncol(dsm_data_masked))
  
  # Calculate boundaries of output raster (boundaries are adjacent to NA values)
  output_boundaries <- raster::boundaries(output)
  
  # Get rows and columns of boundaries cells and convert to list
  xy1 <- raster::xyFromCell(output_boundaries, which(raster::values(output_boundaries) == 1))
  
  # Row/Col of boundary cells
  rc1 <- cbind(raster::rowFromY(dsm_data_masked, xy1[,2]), 
               raster::colFromX(dsm_data_masked, xy1[,1]))
  
  # Apply lineOfSight (C++) function on every point in rc1
  this_LoS <- LoSCpp(rc1 = rc1, r0 = r0, c0 = c0, observerHeight = height0, dsm_mat = dsm_mat)
  
  # Copy result of lineOfSight to the output raster
  output <- raster::crop(output, dsm_data_masked)
  output[this_LoS] <- 1

  return(output)
}
