#' @title Viewshed Greenness Visibility Index (VGVI)
#' @description The VGVI expresses the proportion of visible greenness to the total visible area based on a \code{\link[GVI]{viewshed}}.
#' The estimated VGVI values range between 0 and 1, where 0 = no green cells are visible, and 1 = all of the visible cells are green.
#' A distance decay function is applied, to account for the reducing visual prominence of an object in space with increasing distance from the observer. 
#'  
#' @param viewshed object of class \code{\link[terra]{SpatRaster}}; \code{\link[terra]{SpatRaster}} of the \code{\link[GVI]{viewshed}}
#' @param greenspace object of class \code{\link[terra]{SpatRaster}}; Binary \code{\link[terra]{SpatRaster}} of the Greenspace mask. Values must be 1 for Green and 0 for No-Green
#' @param m numeric; See ‘Details’
#' @param b numeric; See ‘Details’
#' @param mode character; 'logit' or 'exponential'. See ‘Details’
#'
#' @details The type of function, used for calculating the distance decay weights, can be defined with the \code{mode} parameter.
#' The argument 'logit' uses the logistic function, d = 1 / (1 + e^(b * (x - m))) and 'exponential' the exponential function d = 1 / (1 + (b * x^m)).
#' The decay function can be visualized using the \code{\link[GVI]{visualizeWeights}} function.
#' 
#' @return numeric; Weighted VGVI, where 0 = no green cells are visible, and 1 = all of the visible cells are green.
#' @export
#' 
#' @importFrom magrittr %>%
#' @importFrom sf st_coordinates
#' @importFrom terra xyFromCell
#' @importFrom terra cellFromXY
#' @importFrom stats integrate
vgvi <- function(viewshed, greenspace, m = 0.5, b = 8, mode = c("logit", "exponential")) {
  # Get XY coordinates of visible cells
  xy <- viewshed %>% 
    terra::xyFromCell(which(viewshed[] == 1))
  
  # Calculate distance
  centroid <- colMeans(terra::xyFromCell(viewshed, which(!is.na(viewshed[]))))
  dxy = round(sqrt((centroid[1] - xy[,1])^2 + (centroid[2] - xy[,2])^2))
  dxy[dxy==0] = min(dxy[dxy!=0])
  
  # Intersect XY with greenspace mask
  output <- greenspace[terra::cellFromXY(greenspace, xy)] %>% 
    unlist(use.names = FALSE) %>% 
    cbind(dxy, .) %>% 
    na.omit()
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