# Generate parallel coastal survey lines; see california_transects.md for usage.

require_transect_packages <- function() {
  for (package in c("sf", "geosphere")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install the required package with install.packages('", package, "').",
           call. = FALSE)
    }
  }
}

load_california_land <- function(path = NULL) {
  require_transect_packages()
  if (is.null(path)) {
    path <- file.path(tools::R_user_dir("california_transects", "cache"),
                      "ne_10m_land.geojson")
    if (!file.exists(path)) {
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      temporary <- tempfile(tmpdir = dirname(path), fileext = ".geojson")
      on.exit(unlink(temporary), add = TRUE)
      # Pin the public-domain Natural Earth dataset for reproducible designs.
      url <- paste0(
        "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/",
        "ca96624a56bd078437bca8184e78163e5039ad19/geojson/ne_10m_land.geojson"
      )
      utils::download.file(url, temporary, mode = "wb")
      sf::st_read(temporary, quiet = TRUE)
      if (!file.rename(temporary, path)) {
        stop("Could not save the coastline cache: ", path, call. = FALSE)
      }
    }
  }
  sf::st_read(path, quiet = TRUE)
}

snap_to_ring <- function(point, ring) {
  starts <- ring[-nrow(ring), , drop = FALSE]
  changes <- ring[-1, , drop = FALSE] - starts
  lengths <- sqrt(rowSums(changes^2))
  if (any(lengths == 0)) {
    stop("Mainland boundary contains duplicate consecutive vertices.", call. = FALSE)
  }
  offsets <- matrix(point, nrow(starts), 2, byrow = TRUE) - starts
  fractions <- pmax(0, pmin(1, rowSums(offsets * changes) / lengths^2))
  candidates <- starts + changes * fractions
  distances <- sqrt(rowSums(
    (candidates - matrix(point, nrow(starts), 2, byrow = TRUE))^2
  ))
  index <- which.min(distances)
  list(point = candidates[index, ], distance = distances[index],
       chainage = c(0, cumsum(lengths))[index] + fractions[index] * lengths[index])
}

short_coast_section <- function(ring, north, south) {
  lengths <- sqrt(rowSums((ring[-1, ] - ring[-nrow(ring), ])^2))
  positions <- c(0, cumsum(lengths))
  perimeter <- tail(positions, 1)
  forward <- function(start, end) {
    distance <- (end$chainage - start$chainage) %% perimeter
    offsets <- (head(positions, -1) - start$chainage) %% perimeter
    indices <- which(offsets > 1e-7 & offsets < distance - 1e-7)
    indices <- indices[order(offsets[indices])]
    rbind(start$point, ring[indices, , drop = FALSE], end$point)
  }
  if ((south$chainage - north$chainage) %% perimeter <= perimeter / 2) {
    forward(north, south)
  } else {
    section <- forward(south, north)
    section[nrow(section):1, , drop = FALSE]
  }
}

generate_california_transects <- function(north, south, length_nm, spacing_nm,
                                          land = load_california_land(),
                                          snap_tolerance_m = 5000) {
  require_transect_packages()
  check_positive <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) {
      stop(name, " must be a finite positive number.", call. = FALSE)
    }
  }
  check_point <- function(point, name) {
    if (!is.numeric(point) || length(point) != 2L || any(!is.finite(point)) ||
        point[1] < -125 || point[1] > -116 ||
        point[2] < 32.4 || point[2] > 42.1) {
      stop(name, " must be c(longitude, latitude) on the California mainland ",
           "coast (longitude -125 to -116; latitude 32.4 to 42.1).",
           call. = FALSE)
    }
  }
  check_point(north, "north")
  check_point(south, "south")
  if (north[2] <= south[2]) {
    stop("north must be north of south.", call. = FALSE)
  }
  check_positive(length_nm, "length_nm")
  check_positive(spacing_nm, "spacing_nm")
  check_positive(snap_tolerance_m, "snap_tolerance_m")
  if (!inherits(land, c("sf", "sfc")) || is.na(sf::st_crs(land))) {
    stop("land must be sf/sfc land polygons with a defined CRS.", call. = FALSE)
  }
  if (any(!sf::st_geometry_type(land) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("land must contain only POLYGON or MULTIPOLYGON features.", call. = FALSE)
  }

  # California Albers provides one metric frame across both California UTM zones.
  land <- sf::st_transform(sf::st_geometry(land), 4326)
  region <- sf::st_bbox(c(xmin = -130, ymin = 30, xmax = -113, ymax = 44),
                        crs = sf::st_crs(4326))
  land <- sf::st_crop(land, region)
  land <- sf::st_union(sf::st_make_valid(sf::st_transform(land, 3310)))
  parts <- suppressWarnings(sf::st_cast(land, "POLYGON"))
  parts <- parts[!sf::st_is_empty(parts)]
  if (!length(parts)) {
    stop("The land data do not cover California.", call. = FALSE)
  }
  mainland <- parts[which.max(as.numeric(sf::st_area(parts)))]
  ring <- mainland[[1]][[1]][, 1:2, drop = FALSE]
  ring <- ring[c(TRUE, rowSums((ring[-1, ] - ring[-nrow(ring), ])^2) > 0), ,
               drop = FALSE]
  anchors <- sf::st_sfc(sf::st_point(north), sf::st_point(south), crs = 4326)
  anchors <- sf::st_coordinates(sf::st_transform(anchors, 3310))
  north_snap <- snap_to_ring(anchors[1, ], ring)
  south_snap <- snap_to_ring(anchors[2, ], ring)
  if (max(north_snap$distance, south_snap$distance) > snap_tolerance_m) {
    stop("An endpoint is more than ", snap_tolerance_m,
         " m from the mainland shoreline. Supply coastal coordinates or ",
         "better land data.", call. = FALSE)
  }
  baseline <- south_snap$point - north_snap$point
  extent_m <- sqrt(sum(baseline^2))
  if (extent_m < 1) {
    stop("Snapped endpoints must be at least one metre apart.", call. = FALSE)
  }
  along <- baseline / extent_m
  offshore <- c(along[2], -along[1])
  section <- short_coast_section(ring, north_snap, south_snap)
  coast <- sf::st_sfc(sf::st_linestring(section), crs = 3310)
  count <- floor((extent_m + 1e-7) / (spacing_nm * 1852)) + 1
  if (!is.finite(count) || count > 10000) {
    stop("More than 10,000 lines requested; increase spacing_nm.", call. = FALSE)
  }
  stations <- (seq_len(count) - 1) * spacing_nm * 1852
  stations[abs(stations - extent_m) < 1e-7] <- extent_m
  cross_positions <- as.vector(section %*% offshore)
  length_m <- length_nm * 1852
  to_lonlat <- function(point) {
    as.numeric(sf::st_coordinates(sf::st_transform(
      sf::st_sfc(sf::st_point(as.numeric(point)), crs = 3310), 4326
    ))[1, 1:2])
  }

  rows <- vector("list", count)
  tracks <- vector("list", count)
  for (i in seq_len(count)) {
    origin <- north_snap$point + stations[i] * along
    cross_origin <- sum(origin * offshore)
    probe <- sf::st_sfc(sf::st_linestring(rbind(
      origin + (min(cross_positions) - cross_origin - 1) * offshore,
      origin + (max(cross_positions) - cross_origin + 1) * offshore
    )), crs = 3310)
    hits <- sf::st_intersection(coast, probe)
    hits <- hits[!sf::st_is_empty(hits)]
    if (any(!sf::st_geometry_type(hits) %in% c("POINT", "MULTIPOINT"))) {
      stop("Coastline overlaps line ", i, "; adjust study endpoints.",
           call. = FALSE)
    }
    candidates <- if (length(hits)) {
      sf::st_coordinates(hits)[, 1:2, drop = FALSE]
    } else {
      matrix(numeric(), ncol = 2)
    }
    # GEOS can miss an endpoint touch after floating-point projection arithmetic.
    if (i == 1) {
      candidates <- rbind(candidates, north_snap$point)
    }
    if (abs(stations[i] - extent_m) < 1e-7) {
      candidates <- rbind(candidates, south_snap$point)
    }
    if (!nrow(candidates)) {
      stop("No shoreline intersection for line ", i, ".", call. = FALSE)
    }
    # Use the ocean-facing intersection, not an inland bay crossing.
    shore <- candidates[which.max(as.vector(candidates %*% offshore)), ]
    shore_lonlat <- to_lonlat(shore)
    distance_error <- function(distance) {
      geosphere::distGeo(shore_lonlat, to_lonlat(shore + distance * offshore)) -
        length_m
    }
    upper <- length_m * 2
    if (!is.finite(distance_error(upper)) || distance_error(upper) <= 0) {
      stop("Requested length exceeds the supported local projection for line ",
           i, ".", call. = FALSE)
    }
    distance <- stats::uniroot(distance_error, c(0, upper), tol = 1e-6)$root
    end <- shore + distance * offshore
    end_lonlat <- to_lonlat(end)
    tracks[[i]] <- sf::st_linestring(rbind(shore, end))
    rows[[i]] <- data.frame(
      ID = i, shore_lat = shore_lonlat[2], shore_lon = shore_lonlat[1],
      offshore_lat = end_lonlat[2], offshore_lon = end_lonlat[1],
      length_nm = geosphere::distGeo(shore_lonlat, end_lonlat) / 1852,
      baseline_distance_nm = stations[i] / 1852
    )
  }
  lines <- sf::st_sfc(tracks, crs = 3310)
  land_crossing <- vapply(seq_along(lines), function(i) {
    intersection <- sf::st_intersection(lines[i], land)
    sum(as.numeric(sf::st_length(intersection))) > 1
  }, logical(1))
  endpoints <- do.call(rbind, rows)
  endpoints$crosses_land <- land_crossing
  if (any(land_crossing)) {
    warning("Lines crossing land (not shortened): ",
            paste(which(land_crossing), collapse = ", "),
            ". Review the design before use.", call. = FALSE)
  }
  list(endpoints = endpoints, lines = sf::st_sf(endpoints, geometry = lines),
       coast = coast, land = land,
       baseline = sf::st_sfc(sf::st_linestring(rbind(north_snap$point,
                                                   south_snap$point)), crs = 3310))
}

transects_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste(
    "Usage: Rscript california_transects.R",
    "--north-lon=-118.10 --north-lat=33.73",
    "--south-lon=-117.24 --south-lat=32.67",
    "--length-nm=20 --spacing-nm=2",
    "[--output=tracklines.csv] [--land-file=land.gpkg]"
  )
  if (identical(args, "--help")) {
    cat(usage, "\n")
    return(invisible(NULL))
  }
  if (!length(args) || any(!grepl("^--[^=]+=.+$", args))) {
    stop(usage, call. = FALSE)
  }
  keys <- sub("=.*$", "", sub("^--", "", args))
  values <- sub("^--[^=]+=", "", args)
  names(values) <- keys
  required <- c("north-lon", "north-lat", "south-lon", "south-lat",
                "length-nm", "spacing-nm")
  allowed <- c(required, "output", "land-file")
  if (anyDuplicated(keys) || any(!keys %in% allowed) ||
      any(!required %in% keys)) {
    stop("Unknown, duplicated, or missing arguments.\n", usage, call. = FALSE)
  }
  numeric_values <- suppressWarnings(as.numeric(values[required]))
  if (any(!is.finite(numeric_values))) {
    stop("All coordinate, length, and spacing arguments must be numeric.",
         call. = FALSE)
  }
  land <- load_california_land(if ("land-file" %in% keys) values[["land-file"]])
  result <- generate_california_transects(
    north = numeric_values[1:2], south = numeric_values[3:4],
    length_nm = numeric_values[5], spacing_nm = numeric_values[6], land = land
  )
  output <- if ("output" %in% keys) values[["output"]] else "tracklines.csv"
  utils::write.csv(result$endpoints, output, row.names = FALSE)
  message("Wrote ", nrow(result$endpoints), " lines to ", output)
  invisible(result)
}

if (sys.nframe() == 0L) {
  transects_main()
}
