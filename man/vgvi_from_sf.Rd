% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vgvi_from_sf.R
\name{vgvi_from_sf}
\alias{vgvi_from_sf}
\title{Viewshed Greenness Visibility Index (VGVI) from sf}
\usage{
vgvi_from_sf(
  observer,
  dsm_rast,
  dtm_rast,
  greenspace_rast,
  max_distance = 800,
  observer_height = 1.7,
  raster_res = NULL,
  spacing = raster_res,
  m = 0.5,
  b = 8,
  mode = c("logit", "exponential"),
  cores = 1,
  chunk_size = 999,
  folder_path = NULL,
  progress = TRUE
)
}
\arguments{
\item{observer}{object of class \code{sf}; Observer location(s) from where the VGVI should be computed. See ‘Details’ for valid sf geometry types}

\item{dsm_rast}{character; File path to the DSM}

\item{dtm_rast}{character; File path to the DTM}

\item{greenspace_rast}{character; File path to the binary Greenspace mask. Values of the Greenspace mask must be 1 for Green and 0 for No-Green}

\item{max_distance}{numeric; Buffer distance to calculate the viewshed}

\item{observer_height}{numeric > 0; Height of the observer (e.g. 1.7 meters)}

\item{raster_res}{optional; NULL or numeric > 0; Resolution that the GVI raster should be aggregated to. Needs to be a multible of the dsm_rast resolution}

\item{spacing}{optional; numeric > 0; Only if \code{observer} is a linestring (or polygon), points on the line (or on a grid) will be generated.
The \code{spacing} parameter sets the distance in between the points on the line/grid}

\item{m}{numeric; See ‘Details’}

\item{b}{numeric; See ‘Details’}

\item{mode}{character; 'logit' or 'exponential'. See ‘Details’}

\item{cores}{numeric; The number of cores to use, i.e. at most how many child processes will be run simultaneously}

\item{chunk_size}{numeric; Chunk size for parallelization. See ‘Details’}

\item{folder_path}{optional; Folder path to where the output should be saved continuously. Must not inklude a filename extension (e.g. '.shp', '.gpkg').}

\item{progress}{logical; Show progress bar and computation time?}
}
\value{
sf_object containing the weighted VGVI values as POINT features, where 0 = no green cells are visible, and 1 = all of the visible cells are green.
}
\description{
The VGVI expresses the proportion of visible greenness to the total visible area based on a \code{\link[GVI]{viewshed}}.
The estimated VGVI values range between 0 and 1, where 0 = no green cells are visible, and 1 = all of the visible cells are green.
A distance decay function is applied, to account for the reducing visual prominence of an object in space with increasing distance from the observer.
}
\details{
observer needs to be a geometry of type POINT, LINESTRING, MULTILINESTRING, POLYGON or MULTIPOLYGON. If observer is a LINESTRING or MULTILINESTRING, 
points will be generated along the line(s) every "resolution" meters. If observer is a POLYGON or MULTIPOLYGON, a grid with resolution = "resolution" 
will be generated, and VGVI will be computed for every point.
The CRS (\code{\link[sf]{st_crs}}) needs to have a metric unit!

The type of function, used for calculating the distance decay weights, can be defined with the \code{mode} parameter.
The argument 'logit' uses the logistic function, d = 1 / (1 + e^(b * (x - m))) and 'exponential' the exponential function d = 1 / (1 + (b * x^m)).
The decay function can be visualized using the \code{\link[GVI]{visualizeWeights}} function.

Higher values of chunk_size increase computation time, but may also be more RAM intensive. Also, if you choose a high number of cores, RAM usage increases. too.
It is highly recommended to use a Linux or Mac system, for better parallel performance.
}
