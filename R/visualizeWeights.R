#' @title Visualize Weights
#' @description Helper function to visualize the parameters used in the \code{\link[GVI]{vgvi}} function.
#' 
#' @param viewshed object of class \code{\link[terra]{SpatRaster}}; \code{\link[terra]{SpatRaster}} of the \code{\link[GVI]{viewshed}}.
#' @param m numeric; See ‘Details’.
#' @param b numeric; See ‘Details’.
#' @param mode character; 'logit' or 'exponential'. See ‘Details’.
#'
#' @details The type of function, used for calculating the distance decay weights, can be defined with the \code{mode} parameter.
#' The argument 'logit' uses the logistic function, d = 1 / (1 + e^(b * (x - m))) and 'exponential' the exponential function d = 1 / (1 + (b * x^m)).
#' The decay function can be visualized using the \code{\link[GVI]{viewshed}} function.
#' 
#' @export
#' @importFrom magrittr %>%
#' @importFrom sf st_coordinates
#' @importFrom terra xyFromCell
visualizeWeights <- function(viewshed, m = 0.5, b = 8, mode = c("logit", "exponential")) {
  # Get XY coordinates that are visible
  xy <- viewshed %>% 
    terra::xyFromCell(which(viewshed[] == 1))
  
  # Calculate maximum distance
  max_dist = (terra::nrow(viewshed)/2) * res(viewshed)[1]
  
  if (mode == c("logit", "exponential") || mode == "logit") {
    logfun <- function(x){
      return(1 / (1 + exp((b) * (x - m))))
    }
    
    plot_main <- paste0("Mode: logit\nm: ", m, "    b: ", b)
  } else if(mode == "exponential") {
    logfun <- function(x){
      return(1/(1 + ((b) * x^(m))))
    } 
    plot_main <- paste0("Mode: exponential\nm: ", m, "    b: ", b)
  } else {
    stop("Currently only logit and exponential are supported")
  }
  
  plot(logfun(seq(0, 1, length.out = max_dist)), type = "l", 
       ylab = "Decay Weight (d)", xlab = "Distance [m]",
       main = plot_main)
}