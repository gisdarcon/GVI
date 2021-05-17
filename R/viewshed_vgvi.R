#' @importFrom magrittr %>%
#' @importFrom terra rast
#' @importFrom terra vect
#' @importFrom terra crop
#' @importFrom terra mask
#' @importFrom terra xyFromCell
#' @importFrom terra values
#' @importFrom terra ext
#' @importFrom terra crs
#' @importFrom terra res
#' @importFrom terra cellFromXY
#' @importFrom terra rowFromY
#' @importFrom terra colFromX
#' @importFrom terra ncol
#' @importFrom terra boundaries
#' @importFrom sf st_coordinates
#' @importFrom terra xyFromCell
#' @importFrom terra cellFromXY
#' @importFrom stats integrate

viewshed_vgvi <- function(this_aoi, dsm_path, greenspace_path, x0, y0,
                           height0, resolution, m = 0.5, b = 8, mode = c("logit", "exponential")) {
  # Read DSM and Greenspace
  dsm_data <- terra::rast(dsm_path)
  greenspace <- terra::rast(greenspace_path)
  
  #### A: Viewshed ####
  # Crop DSM to AOI
  dsm_data_masked <- terra::crop(dsm_data, terra::vect(this_aoi)) %>% 
    terra::mask(terra::vect(this_aoi))
  
  output <- terra::setValues(dsm_data_masked, 0) %>%
    terra::mask(dsm_data_masked)
  
  # Start row/col
  r0 <- terra::rowFromY(dsm_data_masked, y0)
  c0 <- terra::colFromX(dsm_data_masked, x0)
  
  # Convert output raster to matrix
  dsm_mat <- terra::values(dsm_data_masked) %>% 
    matrix(ncol = terra::ncol(dsm_data_masked))
  
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
  
  viewshed <- output
  rm(output)
  
  #### B: VGVI ####
  
  # Get XY coordinates of visible cells
  xy <- viewshed %>% 
    terra::xyFromCell(this_LoS)
  
  # Calculate distance
  dxy = round(sqrt((x0 - xy[,1])^2 + (y0 - xy[,2])^2))
  dxy[dxy==0] = min(dxy[dxy!=0])
  
  # Intersect XY with greenspace mask
  output <- greenspace[terra::cellFromXY(greenspace, xy)] %>% 
    unlist(use.names = FALSE) %>% 
    cbind(dxy, .)
  colnames(output) <- c("dxy", "visible")
  
  # Get number of green visible cells and total visible cells per distance
  dxyVisibility <- countGroups(xyVisible = output, uniqueXY = unique(output[,1]))
  
  # Proportion of visible green cells
  raw_GVI <- dxyVisibility$visibleGreen / dxyVisibility$visibleTotal
  
  if (length(dxyVisibility$dxy) == 1) {
    return(raw_GVI)
  } else {
    # Calculate weights from logistic function
    if (mode == c("logit", "exponential") || mode == "logit") {
      logfun <- function(x){
        return(1 / (1 + exp(b * (x - m))))
      }
    } else if(mode == "exponential") {
      logfun <- function(x){
        return(1/(1 + (b * x^m)))
      } 
    } else {
      stop("Currently only logit and exponential are supported")
    }
    
    # Normalize distance
    n <- dxyVisibility$dxy / max(dxyVisibility$dxy)
    
    # Calculate weights by taking the proportion of the integral of each step from the integral of the whole area. 
    big_integral <- stats::integrate(logfun, lower = 0, upper = 1)$value
    min_dist <- min(n)
    decayWeights <- lapply(n, function(i){
      stats::integrate(logfun, lower = i-min_dist, upper = i)$value / big_integral
    }) %>% unlist()
    
    # Proportion of visible green
    return(sum(raw_GVI * decayWeights, na.rm = TRUE)) 
  }
}
