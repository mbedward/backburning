#' Generate sampling lines along one or more back-burning line features
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
#' @param step_length A single numeric value for the distance in metres
#'   between sampling lines along each back-burning line feature.
#'
#' @param sample_length A single numeric value specifying the maximum distance
#'   in metres that a sample line will extend either side of a back-burning
#'   line.
#'
#' @param stop_if_returning (logical) If \code{TRUE} (default), each sample line
#'   must be moving away from its associated back-burning line along its whole
#'   length. This is intended to avoid cases where a sample line approaches or
#'   crosses another section of a wiggly back-burning line.
#'
#' @param smoothing_bw A single numeric value for the bandwidth (metres) of the
#'   Gaussian kernel filter used to smooth each back-burning feature. The
#'   default value of 1000m seems to give good results.
#'
#' @return An \code{sf} spatial data frame containing the sampling lines.
#'
#' @export
#
make_sampling_lines <- function(bb_lines,
                                bb_id = "OID",
                                step_length = 500,
                                sample_length = 5000,
                                stop_if_returning = TRUE,
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

  checkmate::assert_number(step_length, finite = TRUE, lower = 1)
  checkmate::assert_number(sample_length, finite = TRUE, lower = 1)

  checkmate::assert_flag(stop_if_returning)

  checkmate::assert_number(smoothing_bw, finite = TRUE, lower = 1)

  res <- lapply(seq_len(nrow(bb_lines)), function(index) {
    # Get the line feature and densify its vertices
    g <- sf::st_geometry(bb_lines[index, ])

    # Prepare a smoothed version of the line feature using Gaussian kernel smoothing
    gsmooth <- smoothr::smooth(g, method = "ksmooth", bandwidth = smoothing_bw)

    # Densify the vertices of the original and smoothed features
    g <- smoothr::smooth(g, method = "densify", max_distance = 1)
    gsmooth <- smoothr::smooth(gsmooth, method = "densify", max_distance = 1)

    # Convert the smoothed feature to a point set. This makes querying the step
    # points a little easier.
    gsmooth_points <- sf::st_cast(gsmooth, "POINT")

    vs <- sf::st_coordinates(g)[, 1:2]
    d <- sapply( seq_len(nrow(vs)-1), function(k) sqrt(sum((vs[k,] - vs[k+1,])^2)) )
    d <- c(0, cumsum(d))

    if (max(d) < step_length) {
      # This feature is too short
      return(NULL)
    }

    isteps <- floor(d / step_length)
    step_points <- match(seq_len(max(isteps)), isteps)

    # Guard against hitting the very end of the line (unlikely but possible)
    nsteps <- length(step_points)
    if (step_points[nsteps] == length(d)) {
      step_points <- step_points[-nsteps]
      nsteps <- nsteps - 1

      if (nsteps == 0) {  # line feature was too short for step distance
        return(NULL)
      }
    }

    step_vertices <- vs[step_points, , drop=FALSE]

    norms <- lapply(seq_len(length(step_points)), function(k) {
      vstep <- step_vertices[k, ]

      pstep <- sf::st_point(vstep) |>
        sf::st_sfc(crs = CRS)

      # Use the nearest vertices from the smoothed line feature
      # to set the angle of the normal vector at the current step.
      pnear <- sf::st_distance(pstep, gsmooth_points)
      inear <- which.min(pnear)

      ibefore <- max(inear-1, 1)
      p0 <- sf::st_coordinates(gsmooth_points[ibefore])[, 1:2]

      iafter <- min(inear+1, length(gsmooth_points))
      p1 <- sf::st_coordinates(gsmooth_points[iafter])[, 1:2]

      dxy <- p1 - p0
      len <- sqrt(sum(dxy^2))
      lenfac <- sample_length / len

      pnorm1 <- vstep + c(dxy[2], -dxy[1]) * lenfac
      pnorm2 <- vstep + c(-dxy[2], dxy[1]) * lenfac

      l1 <- sf::st_linestring(rbind(pnorm1, vstep))
      l2 <- sf::st_linestring(rbind(vstep, pnorm2))

      sf::st_sfc(l1, l2, crs = CRS)
    })

    g <- do.call(c, norms)
    sf::st_sf(index = index, geom = g)
  })

  do.call(rbind, res)
}
