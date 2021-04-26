#' @importFrom stats integrate

vgvi_from_XYV_table <- function(XYV_table, m = 0.5, b = 8, mode = c("logit", "exponential")) {
  # Get number of green visible cells and total visible cells per distance
  dxyVisibility <- countGroups(xyVisible = XYV_table, uniqueXY = unique(XYV_table[,1]))
  
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
      stop("Currently only logit and exponential are supported.")
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