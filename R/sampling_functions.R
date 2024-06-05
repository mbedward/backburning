#' Generate lines of sample points across back-burning line features
#'
#' Given a set of input line features representing back-burning lines, this
#' function places sample lines at regular intervals along each input feature,
#' perpendicular to the local feature angle, and then generates uniformly spaced
#' points along each sample line. Local angles are determined from a smoothed
#' version of the back-burning line to minimize the influence of any local kinks
#' and turns. Sample lines are located at, and either side of, the mid-point of
#' each input feature. At each sampling line, a sample point is located where
#' the line intersects the back-burning line, then further points are placed
#' along the line at uniform spacing. In cases where a back-burning line is
#' curved or convoluted it is possible for sampling lines from one section to
#' approach or cross other sections. Setting the argument
#' \code{increasing_distance} to \code{TRUE} (the default) will test for such
#' cases and prune sampling lines as required. This can result in different
#' numbers of sample points per line being returned.
#'
#' This function calls the non-exported helper function
#' \code{make_sample_lines()} to generate the sample lines along which the
#' sample points are placed.
#'
#' @param bb_lines An \code{sf} spatial data frame containing one or more line
#'   features representing back-burning lines. The data must have a coordinate
#'   reference system defined with metres as map units (e.g. NSW Lambert/GDA94
#'   EPSG:3308).
#'
#' @param bb_id (character; default "OID") The name of a column in the
#'   \code{bb_lines} data frame that uniquely identifies back-burning line
#'   features.
#'
#' @param line_spacing A single numeric value for the distance in metres
#'   between sampling lines along each back-burning line feature.
#'
#' @param line_half_length A single numeric value specifying the maximum
#'   distance in metres that a sampling line will extend on either side of the
#'   back-burning line.
#'
#' @param point_spacing A single numeric value specifying the uniform distance
#'   between sample points along each sampling line. The first point is always
#'   positioned at the intersection of the sampling line and the reference
#'   back-burning line.
#'
#' @param increasing_distance (logical) If \code{TRUE} (default), there must
#'   always be an increasing distance to the reference back-burning line along
#'   the length of each sampling line. A sampling lines will be truncated if
#'   necessary to avoid approaching or crossing another section of the
#'   back-burning line.
#'
#' @param smoothing_bw A single numeric value for the bandwidth (metres) of the
#'   Gaussian kernel filter used to smooth each back-burning feature. The
#'   default value of 1000m seems to give good results.
#'
#' @return An \code{sf} spatial data frame of sample point features with the
#'   following columns:
#'
#'   \code{refid} Identifier of the input back-burning line feature that a point
#'   is sampling, with values taken from the input data column specified by the
#'   \code{bb_id} argument.
#'
#'   \code{index} Integer index (from 1) of the perpendicular sample line on
#'   which the point was positioned.
#'
#'   \code{segment} One of \code{'X'} (point at intersection with back-burning
#'   line); \code{'L'} (point on left of line); \code{'R'} (point on right of
#'   line). Note that left and right are relative to the digitizing direction
#'   (i.e. order of vertices) of the back-burning line.
#'
#'   \code{linedist} Distance of the point along the sample line segment on
#'   which it was positioned, relative to the reference back-burning line.
#'
#'   \code{refdist} Shortest distance from the point to the reference
#'   back-burning line. This can be less than \code{linedist} when the reference
#'   line is wiggly.
#'
#'   \code{geom} Point geometry, projected in the same coordinate reference
#'   system as the input line features.
#'
#' @examples
#' \dontrun{
#' libary(backburning)
#' library(sf)
#'
#' # Load back-burning line features from a GeoPackage layer, shapefile etc.
#' dat_bb <- st_read(...)
#'
#' # Generate points on sampling lines placed at 1km intervals along each
#' # back-burning line. Sampling lines extend 5km either side of the back-burning
#' # line and points are placed every 500m.
#' #
#' dat_sample_points <- make_sample_points(bb,
#'                                         line_spacing = 1000,
#'                                         line_half_length = 5000,
#'                                         point_spacing = 500)
#' }
#'
#' @export
#
make_sample_points <- function(bb_lines,
                               bb_id = "OID",
                               line_spacing,
                               line_half_length,
                               point_spacing,
                               increasing_distance = TRUE,
                               smoothing_bw = 1000) {

  # Most argument checking will be done by the helper line function.
  # Here we just validate the arguments specific to points.
  #
  checkmate::assert_number(point_spacing, lower = 1, finite = TRUE)
  checkmate::assert_flag(increasing_distance)

  # Point spacing must not be more than line segment length
  if (point_spacing > line_half_length) {
    msg <- glue::glue("The value of line_half_length ({line_half_length}) is not large \\
                       enough for the requested point spacing ({point_spacing})")
    stop(msg)
  }

  # Generate sample lines (this step will also check other argument values)
  #
  dat_sample_lines <- make_sample_lines(bb_lines,
                                        bb_id = bb_id,
                                        line_spacing = line_spacing,
                                        line_half_length = line_half_length,
                                        smoothing_bw = smoothing_bw)

  # Generate points on each sample line by placing a point at the
  # intersection with the back-burning line (the start vertex
  # of each sample line segment) then subsequent points along each
  # segment.
  #
  # Point positions expressed as fractions of line segment length
  dpos <- seq(0.0, 1.0, point_spacing / line_half_length)

  # This will return a geometry list of MULTIPOINT objects: one per line segment
  gpoints <- sf::st_line_sample(dat_sample_lines, sample = dpos)

  # Attribute points with back-burning feature ID and sample line values
  dat_points <- dat_sample_lines
  sf::st_geometry(dat_points) <- gpoints

  # Convert from multi-point to single point features
  dat_points <- suppressWarnings({
    sf::st_cast(dat_points, "POINT")
  })

  # For the point where each pair of sample line segments intersect
  # the back-burning line feature, drop the duplicate record and
  # label the remaining record as 'X'. Also add 'linedist' values.
  #
  dat_points <- dat_points %>%
    dplyr::group_by(across(c(refid, index, segment))) %>% # Truly horrible syntax but it seems to work...

    dplyr::mutate(point_index = dplyr::row_number(),
                  linedist = (point_index-1) * point_spacing) %>%

    dplyr::ungroup() %>%
    dplyr::select(-point_index)

  # Calculate the shortest distance of each point from its reference back-burning line
  IDs <- unique(bb_lines[[bb_id]])
  dat_points_dist <- lapply(IDs, function(id) {
    ibb <- which(bb_lines[[bb_id]] == id)
    ip <- which(dat_points$refid == id)
    d <- as.numeric( sf::st_distance(dat_points[ip, ], bb_lines[ibb, ]) )
    d <- zapsmall(d)
    cbind(dat_points[ip, ], refdist = d)
  })

  dat_points <- do.call(rbind, dat_points_dist)

  # Filter points on each sample line segment to enforce increasing distance
  # from the reference back-burning line
  if (increasing_distance) {
    # Function to check point distances
    fn_flag <- function(d) {
      n <- length(d)
      ok <- rep(TRUE, n)

      if (n > 1) {
        dd <- diff(d)
        ok <- c(TRUE, dd > 0)
        x <- which(!ok)
        if (length(x) > 0) ok[min(x):n] <- FALSE
      }

      ok
    }

    dat_points <- dat_points %>%
      dplyr::group_by(across(c(refid, index, segment))) %>%

      dplyr::mutate(keep = fn_flag(refdist)) %>%

      dplyr::ungroup() %>%

      dplyr::filter(keep) %>%
      dplyr::select(-keep)
  }

  # Label intersection points (linedist == 0) and discard duplicate
  # within each line index
  ii <- with(dat_points, which(segment == 'R' & linedist == 0))
  dat_points$segment[ii] <- 'X'

  ii <- with(dat_points, which(segment == 'L' & linedist == 0))
  dat_points <- dat_points[-ii, ]

  # Return points
  dat_points
}


#' Private helper function to generate sample lines along back-burning line features
#'
#' This is a helper function called by \code{make_sample_points()}. Given a set
#' of input line features representing back-burning lines, this function places
#' sample lines at regular intervals along each input feature. Each sampling
#' line is placed perpendicular to the local angle of the back-burning line.
#' Local angles are determined from a smoothed version of the back-burning line
#' to minimize the influence of any local kinks and turns.
#'
#' @param bb_lines An \code{sf} spatial data frame containing one or more line
#'   features representing back-burning lines. The data must have a coordinate
#'   reference system defined with metres as map units (e.g. NSW Lambert/GDA94
#'   EPSG:3308).
#'
#' @param bb_id (character; default "OID") The name of a column in the
#'   \code{bb_lines} data frame that uniquely identifies back-burning line
#'   features.
#'
#' @param line_spacing A single numeric value for the distance in metres
#'   between sample lines along each back-burning line feature.
#'
#' @param line_half_length A single numeric value specifying the maximum
#'   distance in metres that a sample line will extend on either side of the
#'   back-burning line.
#'
#' @param smoothing_bw A single numeric value for the bandwidth (metres) of the
#'   Gaussian kernel filter used to smooth each back-burning feature. The
#'   default value of 1000m seems to give good results.
#'
#' @return An \code{sf} spatial data frame of sample line features with the
#'   following columns:
#'
#'   \code{refid} Identifier of the input back-burning line feature on which the
#'   sample line is positioned, with values taken from the input data column
#'   specified by the \code{bb_id} argument.
#'
#'   \code{index} Integer index (from 1) of the perpendicular sample line on
#'   which the point was positioned.
#'
#'   \code{segment} One of \code{'L'} (sample line extending left from the
#'   reference line); \code{'R'} (sample line extending right from the reference
#'   line). Note that left and right are relative to the digitizing direction
#'   (i.e. order of vertices) of the back-burning line.
#'
#'   \code{geom} Line geometry, projected in the same coordinate reference
#'   system as the input line features.
#'
#' @seealso [make_sample_points()]
#'
#' @noRd
#
make_sample_lines <- function(bb_lines,
                              bb_id = "OID",
                              line_spacing = 500,
                              line_half_length = 5000,
                              smoothing_bw = 1000) {

  checkmate::assert_class(bb_lines, "sf")

  CRS <- sf::st_crs(bb_lines)
  if (is.na(CRS)) stop("A cooordinate reference system must be set for the back-burning line features")

  x <- CRS$units
  if (is.null(x) || x != "m") stop("Map units for back-burning line features must be metres")

  checkmate::assert_string(bb_id, min.chars = 1)
  if (!bb_id %in% colnames(bb_lines)) {
    msg <- glue::glue("The value of argument bb_id {bb_id} is not a column in the input line features data frame")
    stop(msg)
  }

  checkmate::assert_number(line_spacing, finite = TRUE, lower = 1)
  checkmate::assert_number(line_half_length, finite = TRUE, lower = 1)

  checkmate::assert_number(smoothing_bw, finite = TRUE, lower = 1)

  res <- lapply(seq_len(nrow(bb_lines)), function(index) {
    # Get the identifier for this line feature
    FeatureID <- bb_lines[[bb_id]][index]

    # Get the line feature and densify its vertices
    g <- sf::st_geometry(bb_lines[index, ])

    # Prepare a smoothed version of the line feature using Gaussian kernel smoothing
    gsmooth <- smoothr::smooth(g, method = "ksmooth", bandwidth = smoothing_bw)

    # Densify the vertices of the original and smoothed features
    VertexDistance <- 1.0
    g <- smoothr::smooth(g, method = "densify", max_distance = VertexDistance)
    gsmooth <- smoothr::smooth(gsmooth, method = "densify", max_distance = VertexDistance)

    # Convert the smoothed feature to a point set. This makes querying the step
    # points a little easier.
    #gsmooth_points <- sf::st_cast(gsmooth, "POINT")
    gsmooth_vertices <- sf::st_coordinates(gsmooth)[, 1:2]

    vs <- sf::st_coordinates(g)[, 1:2]
    d <- sapply( seq_len(nrow(vs)-1), function(k) sqrt(sum((vs[k,] - vs[k+1,])^2)) )
    d <- c(0, cumsum(d))

    # Locate the vertex closest to the middle of the feature
    dmax <- max(d)
    imid <- which.min( abs(d - dmax/2) )

    # It's just a jump to the left...
    step_points <- icur <- imid
    repeat {
      dtarget <- d[icur] - line_spacing
      if (dtarget > line_spacing/2 + VertexDistance) {
        icur <- which.min(abs(d - dtarget))
        step_points <- c(step_points, icur)
      } else {
        break
      }
    }

    # ...and then a jump to the right
    icur <- imid
    repeat {
      dtarget <- d[icur] + line_spacing
      if (dtarget < dmax - (line_spacing/2 + VertexDistance)) {
        icur <- which.min(abs(d - dtarget))
        step_points <- c(step_points, icur)
      } else {
        break
      }
    }

    # Get step vertices, ordered by the direction of back-burning line vertices
    step_points <- sort(step_points)
    step_vertices <- vs[step_points, , drop=FALSE]

    # Create the pair of sample line segments at each step vertex, either side
    # of the back-burning line
    sample_lines <- lapply(seq_len(length(step_points)), function(istep) {
      vcur <- step_vertices[istep, ]

      # Use the nearest vertices from the smoothed line feature
      # to set the angle of the normal vector at the current step.
      dsmooth2 <- apply(gsmooth_vertices, 1, function(vxy) sum((vxy - vcur)^2))
      inear <- which.min(dsmooth2)

      ibefore <- max(inear-1, 1)
      p0 <- gsmooth_vertices[ibefore, ]

      iafter <- min(inear+1, nrow(gsmooth_vertices))
      p1 <- gsmooth_vertices[iafter, ]

      dxy <- p1 - p0
      len <- sqrt(sum(dxy^2))
      lenfac <- line_half_length / len

      pnorm1 <- vcur + c(dxy[2], -dxy[1]) * lenfac
      pnorm2 <- vcur + c(-dxy[2], dxy[1]) * lenfac

      left_seg <- sf::st_linestring(rbind(vcur, pnorm1))
      right_seg <- sf::st_linestring(rbind(vcur, pnorm2))

      # Geom list with the two line segments, labelled 'R' (right) and
      # 'L' (left) relative to the order of vertices of the target feature
      segments <- sf::st_sfc(left_seg, right_seg, crs = CRS)
      sf::st_sf(refid = FeatureID, index = istep, segment = c('R', 'L'), geom = segments)
    })

    do.call(rbind, sample_lines)
  })

  # Combine sets of sample lines into a single sf data frame
  res <- do.call(rbind, res)

  # Attributes for client code use
  attr(res, "line_spacing") <- line_spacing
  attr(res, "line_half_length") <- line_half_length

  # Return sample lines
  res
}
