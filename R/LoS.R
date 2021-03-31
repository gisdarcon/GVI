#' @title Line of Sight
#' @description Computes a single Line of Sight (LoS) using a Digital Surface Model (DSM) matrix. Uses a C++ implementation of Bresenham's line algorithm.
#'
#' @param rc1 integer vector; Row and column of the targets location in the DSM matrix.
#' @param r0 integer; Specifying the row of the starting location in the DSM matrix.
#' @param c0 integer; Specifying the column of the starting location in the DSM matrix.
#' @param dsm_mat integer matrix; Matrix derived from a DSM.
#' @param observerHeight numeric > 0; Height of the observer (e.g. 1.8)
#'
#' @return Cell numbers of all visible cells of the DSM. 
#' @export

LoS <- function(rc1, r0, c0, dsm_mat, observerHeight) {
  # End points
  r1 <- rc1[1]
  c1 <- rc1[2]
  
  # Get XY coordinates of all points in the matrix along a line from r0,c0 to r1,c1
  pointsZ <- bresenham(x1 = c0, y1 = r0, x2 = c1, y2 = r1)
  pointsZmat <- matrix(pointsZ, ncol = 2, byrow = T)
  
  # Pull DSM height values from matrix
  dsm_profile <- dsm_mat[pointsZmat]
  
  # Is visible? Current tangent must be greater than max. tangent
  tangents <- tangents(x1 = c0, y1 = r0, height0 = observerHeight, xy2 = pointsZ, dsm_profile = dsm_profile)
  visibile <- isVisible(tangents)
  
  # Get visible cell numbers
  return((pointsZmat[visibile, 2]-1)*nrow(dsm_mat) + pointsZmat[visibile, 1])
}