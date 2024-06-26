#' Sample date-time values from fire progression polygons
#'
#' @param prog An \code{sf} spatial data frame of fire progression polygons. The
#'   data must have a coordinate reference system defined and the progression
#'   polygons must not overlap (both of these conditions will be checked).
#'
#' @param pts An \code{sf} spatial data frame of sample points as generated by
#'   function \code{make_sample_points()}. It is assumed that each sample point
#'   has attributes to relate it to a back-burning line.
#'
#' @param prog_datetime The name of the column in the input \code{progressions} data
#'   frame with date-time values.
#'
#' @param formats WRITE ME!!!
#'
#' @param prog_tz WRITE ME!!!
#'
#' @return A copy of the input \code{sample_points} spatial data frame with an
#'   additional column \code{'prog_datetime'} containing the date-time values.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(backburning)
#'
#' # Load a layer of back-burning lines and generate sample points at the
#' # desired spacing
#' dat_bb <- st_read("path/to/bb_lines")
#' dat_pts <- make_sample_points(dat_bb, bb_id = "OID",
#'                               line_spacing = 1000,
#'                               line_half_length = 5000,
#'                               point_spacing = 500)
#'
#' # Load a corresponding layer of fire progression polygons and get
#' # the date-time value for each sample point
#' dat_prog <- st_read("path/to/progressions")
#'
#' dat_pts <- get_progression_times(dat_prog, dat_pts, prog_datetime = "LASTTIME")
#' }
#'
#' @export
#
get_progression_times <- function(prog, pts,
                                  prog_datetime,
                                  formats = c("ymdHMS", "ymdHM", "ymd"),
                                  prog_tz = "Australia/NSW") {

  # Name to use for output date-time values appended to the sample points data frame
  OutTimeName <- "prog_datetime"

  checkmate::assert_class(prog, "sf")

  if (!all(sf::st_geometry_type(prog) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("All progression features must be either POLYGON or MULTIPOLYGON geometries")
  }

  CRSprog <- sf::st_crs(prog)
  if (is.na(CRSprog)) stop("A cooordinate reference system must be set for the progression layer")

  checkmate::assert_string(prog_datetime, min.chars = 1)
  if (!(prog_datetime %in% colnames(prog))) {
    msg <- glue::glue("The value of argument prog_datetime {prog_datetime} is not a column in the progressions data")
    stop(msg)
  }

  # Check the type and format of the data-time values
  dtcl <- class(prog[[prog_datetime]])

  # Note: using `any` here because POSIX date-times will have a length-2 class vector, e.g. c("POSIXct", "POSIXt")
  if (any(dtcl %in% c("POSIXt", "Date"))) {
    FnDateTime <- base::identity

  } else if (dtcl == "character") {
    if (missing(prog_tz) || is.null(prog_tz)) {
      tz <- "Australia/NSW"
      warning("No time zone specified for argument tz. Defaulting to ", tz,
              immediate. = TRUE)
    }

    FnDateTime <- function(x) {
      lubridate::parse_date_time(x,
                                 orders = formats,
                                 quiet = TRUE,
                                 tz = prog_tz)
    }

  } else {
    msg <- glue::glue("{dtcl} is not a supported class for date-time values.
                      Valid options are POSIXt, Date or character")
    stop(msg)
  }

  # Check sample points layer
  checkmate::assert_class(pts, "sf")

  if (!all(sf::st_geometry_type(pts) == "POINT")) {
    stop("All sample points must be geometry type POINT")
  }

  CRSpts <- sf::st_crs(pts)
  if (is.na(CRSpts)) stop("A cooordinate reference system must be set for the sample points layer")

  # Make a temporary copy of the points and align coordinate reference system to
  # that of the progressions if necessary
  pts_query <- sf::st_sf(index__ = 1:nrow(pts), geom = st_geometry(pts))
  if (CRSpts != CRSprog) {
    pts_query <- st_transform(pts_query, CRSprog)
  }

  # Intersect with progression polygons and check that there are no points that
  # intersect more than one polygon
  xii <- sf::st_intersects(pts_query, prog)
  xn <- lengths(xii)
  if (any(xn > 1)) {
    n <- sum(xn > 1)

    if (n > 1) {
      msg <- glue::glue("{n} points intersect more than one progression polygon")
    } else {
      msg <- glue::glue("{n} point intersects more than one progression polygon")
    }
    msg <- glue::glue("{msg}.
                      Only the first intersection will be returned.")

    warning(msg, immediate. = TRUE)
  }

  # Note: can't just use sapply, vapply etc here because of the (known) problem of
  # vector concatenation dropping S3 class attributes which causes POSIX values to
  # be reduced to numeric. Instead we use the hacky work-around of do.call(c, ...)
  # that is described here: https://stackoverflow.com/a/70736094/40246
  #
  times <- do.call(c,
                   lapply(xii, FUN = function(xi) {
                     # Only use first value in xi if more than one
                     FnDateTime( prog[[prog_datetime]][xi[1]] )
                   })
  )


  # Append progression times to the point data in the original CRS
  pts[[OutTimeName]] <- times

  # Put point geometry column last as per usual 'sf' convention
  gname <- attr(pts, which = "sf_column", exact = TRUE)
  if (!is.null(gname)) {
    cn <- c(setdiff(colnames(pts), gname), gname)
    pts <- pts[, cn]
  }

  pts
}


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
#' @param increasing_distance (character) One of \code{('all', 'parent',
#'   'none')}. May be abbreviated. If \code{'all'} (default), each sample point
#'   must lie at a greater distance to all back-burning lines than the preceding
#'   point AND must be closer to the parent back-burning line than any other; if
#'   \code{'parent'}, a similar increasing distance rule is applied but only in
#'   relation to the parent back-burning line that the sample points are
#'   associated with; if \code{'none'}, no distance rule applies.
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
#'   \code{line_index} Integer index (from 1) of the perpendicular sample line
#'   on which the point was positioned.
#'
#'   \code{segment} One of \code{'X'} (point at intersection with back-burning
#'   line); \code{'L'} (point on left of line); \code{'R'} (point on right of
#'   line). Note that left and right are relative to the digitizing direction
#'   (i.e. order of vertices) of the back-burning line.
#'
#'   \code{dist_line} Distance of the point along the sample line segment on
#'   which it was positioned, relative to the parent back-burning line.
#'
#'   \code{dist_parent} Shortest distance from the point to the parent
#'   back-burning line. This can be less than \code{dist_line} if the parent
#'   line is wiggly.
#'
#'   \code{dist_all} Shortest distance from the point to \strong{any}
#'   back-burning line. This column will only be present if the argument
#'   \code{increasing_distance == 'all'}.
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
                               increasing_distance = c("all", "parent", "none"),
                               smoothing_bw = 1000) {

  # Most argument checking will be done by the helper line function.
  # Here we just validate the arguments specific to points.
  #
  checkmate::assert_number(point_spacing, lower = 1, finite = TRUE)

  checkmate::assert_string(increasing_distance)
  increasing_distance = match.arg(increasing_distance)

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

  # Cumulative distance along the line segment for each point.
  dat_points <- dat_points %>%
    dplyr::group_by(across(c(refid, line_index, segment))) %>% # Truly horrible syntax but it seems to work...

    dplyr::mutate(point_index = dplyr::row_number(),
                  dist_line = (point_index-1) * point_spacing) %>%

    dplyr::ungroup() %>%
    dplyr::select(-point_index)

  # Shortest distance of each point to the parent back-burning line.
  # This can be less than 'dist_line' (calculated above) if the parent line
  # is wiggly.
  dat_points$dist_parent <- NA_real_
  for (xid in unique(bb_lines[[bb_id]])) {
    ibb <- bb_lines[[bb_id]] == xid
    ipoints <- dat_points$refid == xid
    dat_points$dist_parent[ipoints] <- as.numeric( sf::st_distance(dat_points[ipoints, ], bb_lines[ibb, ]) )
  }

  dat_points$dist_parent <- zapsmall(dat_points$dist_parent)

  # Shortest distance to *any* back-burning line (i.e. not just the parent)
  # required if filtering points on all distances.
  #
  if (increasing_distance == 'all') {
    dat_points$dist_all <- NA_real_

    for (xid in unique(bb_lines[[bb_id]])) {
      ibb <- bb_lines[[bb_id]] == xid
      ipoints <- dat_points$refid == xid

      # Find any other lines that are close enough to the parent line to
      # be relevant for point distances (trying to reduce run time)
      #
      xbuf <- sf::st_buffer(bb_lines[ibb, ], dist = 2 * line_half_length)
      xi <- which( lengths(sf::st_intersects(bb_lines, xbuf)) > 0)

      d <- sf::st_distance(dat_points[ipoints, ], bb_lines[xi, ])  # matrix where ncol = number of relevant bb lines
      d <- apply(d, MARGIN = 1, min)
      dat_points$dist_all[ipoints] <- as.numeric(d)
    }

    dat_points$dist_all <- zapsmall(dat_points$dist_all)
  }

  # If a point distance filter is requested filter points on each sample line
  # segment to ensure that the nearest distance to either the reference
  # back-burning line (increasing_distance == 'parent') or all back-burning
  # lines (increasing_distance == 'all') for each point is greater than for the
  # preceding point as we proceed along the segment.
  #
  if (increasing_distance != 'none') {
    # Function to check the distance values of points along an individual sample
    # line segment. Returns a logical vector indicating which points to keep.
    # Either considers just the distance to the parent feature (if dall is NULL)
    # or the distance to both parent and other features (if dall is not NULL)
    #
    fn_keep <- function(dparent, dall = NULL) {
      n <- length(dparent)
      delta_parent <- diff(dparent)
      keep <- c(TRUE, delta_parent > 0)  # will work even if length(d) < 2 so no deltas

      if (!is.null(dall)) {  # considering both distance to parent and other features
        delta_all <- diff(dall)
        keep <- keep & c(TRUE, delta_all > 0) & dall >= dparent
      }

      x <- which(!keep)
      if (length(x) > 0) keep[min(x):n] <- FALSE

      keep
    }

    dat_points <- dat_points %>%
      dplyr::group_by(across(c(refid, line_index, segment)))

    if (increasing_distance == 'parent') {
      # Filter on values of distance to parent feature
      dat_points <- dplyr::mutate(dat_points, keep = fn_keep(dist_parent))
    } else {
      # Filter on values of distance to any feature
      dat_points <- dplyr::mutate(dat_points, keep = fn_keep(dist_parent, dist_all))
    }

    dat_points <- dat_points %>%
      dplyr::ungroup() %>%
      dplyr::filter(keep) %>%
      dplyr::select(-keep)
  }

  # Label intersection points (dist_line == 0) and discard duplicate
  # within each line index
  ii <- with(dat_points, which(segment == 'R' & dist_line == 0))
  dat_points$segment[ii] <- 'X'

  ii <- with(dat_points, which(segment == 'L' & dist_line == 0))
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
#'   \code{line_index} Integer index (from 1) of the perpendicular sample line
#'   on which the point was positioned.
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
  if (!(bb_id %in% colnames(bb_lines))) {
    msg <- glue::glue("The value of argument bb_id {bb_id} is not a column in the input line features data frame")
    stop(msg)
  }

  # Line IDs must be unique
  num_ids <- length(unique(bb_lines[[bb_id]]))
  if (num_ids != nrow(bb_lines)) {
    msg <- glue::glue("The values of the specified ID column {bb_id} are not all unique")
    stop(msg)
  }

  checkmate::assert_number(line_spacing, finite = TRUE, lower = 1)
  checkmate::assert_number(line_half_length, finite = TRUE, lower = 1)

  checkmate::assert_number(smoothing_bw, finite = TRUE, lower = 1)

  res <- lapply(seq_len(nrow(bb_lines)), function(line_index) {
    # Get the identifier for this line feature
    FeatureID <- bb_lines[[bb_id]][line_index]

    # Get the line feature and densify its vertices
    g <- sf::st_geometry(bb_lines[line_index, ])

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
      sf::st_sf(refid = FeatureID, line_index = istep, segment = c('R', 'L'), geom = segments)
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
