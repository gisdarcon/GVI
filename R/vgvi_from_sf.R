#' Viewshed Greenness Visibility Index (VGVI) from sf
#' @description The VGVI expresses the proportion of visible greenness to the total visible area based on a \code{\link[GVI]{viewshed}}.
#' The estimated VGVI values range between 0 and 1, where 0 = no green cells are visible, and 1 = all of the visible cells are green.
#' A distance decay function is applied, to account for the reducing visual prominence of an object in space with increasing distance from the observer.
#'
#' @param observer object of class \code{sf}; Observer location(s) from where the VGVI should be computed. See ‘Details’ for valid sf geometry types
#' @param dsm_rast character; File path to the DSM
#' @param dtm_rast character; File path to the DTM
#' @param greenspace_rast character; File path to the binary Greenspace mask. Values of the Greenspace mask must be 1 for Green and 0 for No-Green
#' @param max_distance numeric; Buffer distance to calculate the viewshed
#' @param observer_height numeric > 0; Height of the observer (e.g. 1.7 meters)
#' @param raster_res optional; NULL or numeric > 0; Resolution that the GVI raster should be aggregated to. Needs to be a multible of the dsm_rast resolution
#' @param spacing optional; numeric > 0; Only if \code{observer} is a linestring (or polygon), points on the line (or on a grid) will be generated.
#' The \code{spacing} parameter sets the distance in between the points on the line/grid
#' @param m numeric; See ‘Details’
#' @param b numeric; See ‘Details’
#' @param mode character; 'logit' or 'exponential'. See ‘Details’
#' @param cores numeric; The number of cores to use, i.e. at most how many child processes will be run simultaneously
#' @param chunk_size numeric; Chunk size for parallelization. See ‘Details’
#' @param folder_path optional; Folder path to where the output should be saved continuously. Must not inklude a filename extension (e.g. '.shp', '.gpkg').
#' @param progress logical; Show progress bar and computation time?
#'
#' @details 
#' observer needs to be a geometry of type POINT, LINESTRING, MULTILINESTRING, POLYGON or MULTIPOLYGON. If observer is a LINESTRING or MULTILINESTRING, 
#' points will be generated along the line(s) every "resolution" meters. If observer is a POLYGON or MULTIPOLYGON, a grid with resolution = "resolution" 
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
#' @importFrom sf st_as_sfc
#' @importFrom sf st_cast
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
#' @importFrom terra writeRaster
#' @importFrom terra rast
#' @importFrom methods is
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%

vgvi_from_sf <- function(observer, dsm_rast, dtm_rast, greenspace_rast,
                         max_distance = 800, observer_height = 1.7,
                         raster_res = NULL, spacing = raster_res,
                         m = 0.5, b = 8, mode = c("logit", "exponential"),
                         cores = 1, chunk_size = 999,
                         folder_path = NULL , progress = TRUE) {
  
  #### 1. Check input ####
  # observer
  valid_sf_types <- c("POINT", "MULTIPOINT", "LINESTRING", "MULTILINESTRING", "POLYGON", "MULTIPOLYGON")
  if (!is(observer, "sf")) {
    stop("observer must be a sf object")
  } else if (is.null(sf::st_crs(observer)$units) || sf::st_crs(observer)$units != "m") {
    stop("observer CRS unit needs to be metric")
  } else if (!as.character(sf::st_geometry_type(observer, by_geometry = FALSE)) %in% valid_sf_types) {
    stop("observer has no valid geometry")
  } else if (as.character(sf::st_geometry_type(observer, by_geometry = FALSE)) == "MULTIPOINT") {
    observer <- sf::st_cast(observer, "POINT")
  }
  
  # dsm_rast
  if (!is(dsm_rast, "SpatRaster")) {
    stop("dsm_rast needs to be a SpatRaster object!")
  } else if (sf::st_crs(terra::crs(dsm_rast))$epsg != sf::st_crs(observer)$epsg) {
    stop("dsm_rast needs to have the same CRS as observer")
  }
  
  # dtm_rast
  if (!is(dtm_rast, "SpatRaster")) {
    stop("dtm_rast needs to be a SpatRaster object!")
  } else if (sf::st_crs(terra::crs(dtm_rast))$epsg != sf::st_crs(observer)$epsg) {
    stop("dtm_rast needs to have the same CRS as observer")
  }
  
  # greenspace_rast
  if (!is(greenspace_rast, "SpatRaster")) {
    stop("greenspace_rast needs to be a SpatRaster object!")
  } else if (sf::st_crs(terra::crs(greenspace_rast))$epsg != sf::st_crs(observer)$epsg) {
    stop("greenspace_rast needs to have the same CRS as observer")
  }
  
  # raster_res
  dsm_res <- min(terra::res(dsm_rast))
  if (is.null(raster_res)) {
    raster_res = dsm_res
  } else if (raster_res < min(terra::res(dsm_rast))) {
    stop("raster_res must be higher than the resolution of dsm_rast")
  } else if ((raster_res %% dsm_res) != 0) {
    stop(paste0("raster_res must be a multible of the dsm_rast resolution. Try raster_res = ", raster_res - (raster_res %% dsm_res)))
  }
  rm(dsm_res)
  
  # spacing
  if (is.null(spacing)) {
    spacing <- raster_res
  }
  
  # folder_path
  if (!is.null(folder_path)) {
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    folder_path <- tempfile(pattern = "VGVI_", tmpdir = file.path(folder_path), fileext = ".gpkg")
  }
  
  
  if(progress) {
    message("Preprocessing:")
    pb = txtProgressBar(min = 0, max = 5, initial = 0, style = 2)
  }
  
  
  #### 2. Convert observer to points ####
  if (as.character(sf::st_geometry_type(observer, by_geometry = FALSE)) %in% c("LINESTRING", "MULTILINESTRING")) {
    observer <- observer %>%
      sf::st_union() %>%
      sf::st_cast("LINESTRING") %>%
      sf::st_line_sample(density = 1/spacing) %>%
      sf::st_cast("POINT") %>%
      sf::st_as_sf() %>% 
      dplyr::rename(geom = x)
  } else if (as.character(sf::st_geometry_type(observer, by_geometry = FALSE)) %in% c("POLYGON", "MULTIPOLYGON")) {
    observer_bbox <- sf::st_bbox(observer)
    observer <- terra::rast(xmin = observer_bbox[1], xmax = observer_bbox[3], 
                            ymin = observer_bbox[2], ymax = observer_bbox[4], 
                            crs = terra::crs(dsm_rast), resolution = spacing, vals = 0) %>% 
      terra::crop(terra::vect(observer)) %>% 
      terra::mask(terra::vect(observer)) %>%
      terra::xyFromCell(which(terra::values(.) == 0)) %>%
      as.data.frame() %>% 
      sf::st_as_sf(coords = c("x","y"), crs = sf::st_crs(observer)) %>% 
      dplyr::rename(geom = geometry)
  }
  if (progress) setTxtProgressBar(pb, 1)
  
  
  #### 3. Prepare data for viewshed analysis ####
  # Max AOI
  max_aoi <- observer %>% 
    sf::st_bbox() %>% 
    sf::st_as_sfc() %>% 
    sf::st_buffer(max_distance)
  
  # Crop DSM to max AOI and change resolution
  if(raster_res != min(raster::res(dsm_rast))) {
    dsm_rast <- terra::crop(dsm_rast, terra::vect(max_aoi))
    
    terra::terraOptions(progress = 0)
    dsm_rast <- terra::aggregate(dsm_rast, fact = raster_res/terra::res(dsm_rast))
    terra::terraOptions(progress = 3)
  }
  
  # Coordinates of start point
  x0 <- sf::st_coordinates(observer)[,1]
  y0 <- sf::st_coordinates(observer)[,2]
  
  # Observer heights
  height_0_vec <- unlist(terra::extract(dtm_rast, cbind(x0, y0)), use.names = F) + observer_height
  
  if (progress) setTxtProgressBar(pb, 2)
  
  
  #### 4. Remove points outside the DSM or DTM ####
  invalid_points <- unique(c(
    which(is.na(terra::extract(dsm_rast, cbind(x0, y0)))), # points outside the DSM
    which(is.na(height_0_vec)) # points outside the DTM
  ))
  
  # Remove invalid points
  if (length(invalid_points) > 0) {
    observer <- observer[-invalid_points, ]
    x0 <- x0[-invalid_points]
    y0 <- y0[-invalid_points]
    height_0_vec <- height_0_vec[-invalid_points]
  }
  
  if (progress) setTxtProgressBar(pb, 3)
  
  
  #### 5. Last steps of PreProcessing ####
  # Prepare observer for output
  observer <- observer %>% 
    dplyr::mutate(VGVI = as.numeric(NA),
                  id = 1:dplyr::n()) %>% 
    dplyr::select(id, VGVI)
  
  # Convert to list
  observer_list <- suppressWarnings(
    1:nrow(observer) %>%
      split(seq(1, length(.), chunk_size))
  )
  
  if (progress) setTxtProgressBar(pb, 4)
  
  
  # 6. Save DSM and Greenspace as temp_file
  dsm_path <- tempfile("temp_dsm", fileext = ".tif")
  terra::writeRaster(dsm_rast, dsm_path)
  
  greenspace_path <- tempfile("temp_GS", fileext = ".tif")
  terra::writeRaster(greenspace_rast, greenspace_path)
  
  
  # 7. Final update of Pre-Processing ProgressBar
  if (progress) {
    setTxtProgressBar(pb, 5)
    cat("\n")
  }
  if (length(invalid_points) == 1) {
    message("1 point has been removed, because it was outside of the DSM or DTM")
  } else if (length(invalid_points) > 1) {
    message(paste(length(invalid_points), "points have been removed, because they were outside of the DSM or DTM"))
  }
  
  
  #### 8. Calculate viewsheds and VGVI ####
  if (progress) {
    message(paste0("Computing VGVI for ", nrow(observer), " points:"))
    pb = txtProgressBar(min = 0, max = length(observer_list), initial = 0, style = 3)
    start_time <- Sys.time()
  }
  
  if (cores > 1 && Sys.info()[["sysname"]] == "Windows") {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  }
  
  for (j in seq_along(observer_list)) {
    
    this_aoi <- observer[observer_list[[j]], ] %>% 
      sf::st_buffer(max_distance)
    this_ids <- this_aoi$id
    this_x0 <- x0[this_ids]
    this_y0 <- y0[this_ids]
    this_height_0_vec <- height_0_vec[this_ids]
    
    #### Calculate viewshed
    if (cores > 1) {
      if (Sys.info()[["sysname"]] == "Windows") {
        # Calculate VGVI in parallel
        par_fun <-  function(i){
          viewshed_vgvi(this_aoi = this_aoi[i,],
                        dsm_path = dsm_path, greenspace_path = greenspace_path,
                        x0 = this_x0[i], y0 = this_y0[i], height0 = this_height_0_vec[i],
                        resolution = raster_res, m = m, b = b, mode = mode)
        }
        
        vgvi_list <- foreach::foreach(i=1:nrow(this_aoi)) %dopar% par_fun(i)
      }
      else {
        vgvi_list <- parallel::mclapply(1:nrow(this_aoi), function(i){
          viewshed_vgvi(this_aoi = this_aoi[i,],
                        dsm_path = dsm_path, greenspace_path = greenspace_path,
                        x0 = this_x0[i], y0 = this_y0[i], height0 = this_height_0_vec[i],
                        resolution = raster_res, m = m, b = b, mode = mode)},
          mc.cores = cores, mc.preschedule = TRUE)
      }
    } else {
      vgvi_list <- lapply(1:nrow(this_aoi), function(i){
        viewshed_vgvi(this_aoi = this_aoi[i,],
                      dsm_path = dsm_path, greenspace_path = greenspace_path,
                      x0 = this_x0[i], y0 = this_y0[i], height0 = this_height_0_vec[i],
                      resolution = raster_res,  m = m, b = b, mode = mode)})
    }
    
    # Update observer
    observer[this_ids,2] <- unlist(vgvi_list, use.names = FALSE)
    
    if (!is.null(folder_path)) {
      sf::st_write(observer[this_ids, ], folder_path, append = TRUE, quiet = T)
    }
    
    # Update ProgressBar
    if (progress) setTxtProgressBar(pb, j)
  }
  
  if (cores > 1 && Sys.info()[["sysname"]] == "Windows") {
    parallel::stopCluster(cl)
  }
  
  if (progress) {
    cat("\n")
    time_dif <- round(cores * (as.numeric(difftime(Sys.time(), start_time, units = "s")) / nrow(observer)), 2)
    message(paste("Total runtime:", round(as.numeric(difftime(Sys.time(), start_time, units = "m")))), " mins")
    message(paste("Average time for a single point:", time_dif, "secs"))
  }
  return(observer)
}