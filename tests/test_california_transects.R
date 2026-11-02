source("california_transects.R")
require_transect_packages()

expect_error <- function(expression, pattern) {
  error <- tryCatch({
    force(expression)
    NULL
  }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

lonlat <- function(point) {
  as.numeric(sf::st_coordinates(sf::st_transform(
    sf::st_sfc(sf::st_point(point), crs = 3310), 4326
  ))[1, 1:2])
}

north_xy <- c(0, 100000)
south_xy <- c(20000, 0)
north <- lonlat(north_xy)
south <- lonlat(south_xy)
mainland <- sf::st_polygon(list(rbind(
  c(-10000, 150000), c(120000, 150000), c(120000, -50000),
  c(30000, -50000), c(-10000, 150000)
)))
land <- sf::st_sfc(mainland, crs = 3310)
design <- generate_california_transects(north, south, 20, 2, land)
xy <- sf::st_coordinates(design$lines)[, 1:2]
starts <- xy[seq(1, nrow(xy), 2), , drop = FALSE]
ends <- xy[seq(2, nrow(xy), 2), , drop = FALSE]
along <- (south_xy - north_xy) / sqrt(sum((south_xy - north_xy)^2))
offshore <- c(along[2], -along[1])
vectors <- ends - starts
stopifnot(
  nrow(design$endpoints) == floor(sqrt(sum((south_xy - north_xy)^2)) / 3704) + 1,
  max(abs(as.vector(vectors %*% along))) < 1e-5,
  all(as.vector(vectors %*% offshore) > 0),
  max(abs(diff(as.vector(starts %*% along)) - 3704)) < 1e-5,
  max(abs(design$endpoints$length_nm * 1852 - 37040)) < 1e-4,
  max(abs(as.numeric(sf::st_distance(
    sf::st_sfc(lapply(seq_len(nrow(starts)), function(i) sf::st_point(starts[i, ])),
               crs = 3310), design$coast
  )))) < 1e-5,
  !any(design$endpoints$crosses_land),
  max(design$endpoints$baseline_distance_nm) * 1852 <=
    sqrt(sum((south_xy - north_xy)^2)),
  identical(names(design$endpoints), c(
    "ID", "shore_lat", "shore_lon", "offshore_lat", "offshore_lon",
    "length_nm", "baseline_distance_nm", "crosses_land"
  ))
)

# An island must be reported without changing the requested line length.
centre <- north_xy + along * 37040 + offshore * 10000
island_xy <- rbind(
  centre - along * 1000 - offshore * 1000,
  centre + along * 1000 - offshore * 1000,
  centre + along * 1000 + offshore * 1000,
  centre - along * 1000 + offshore * 1000,
  centre - along * 1000 - offshore * 1000
)
with_island <- sf::st_sfc(mainland, sf::st_polygon(list(island_xy)), crs = 3310)
warnings <- character()
island_design <- withCallingHandlers(
  generate_california_transects(north, south, 20, 2, with_island),
  warning = function(warning) {
    warnings <<- c(warnings, conditionMessage(warning))
    invokeRestart("muffleWarning")
  }
)
stopifnot(
  island_design$endpoints$crosses_land[11],
  any(grepl("Lines crossing land", warnings)),
  max(abs(island_design$endpoints$length_nm - 20)) < 1e-7
)

single <- generate_california_transects(north, south, 1, 100, land)
stopifnot(nrow(single$endpoints) == 1, single$endpoints$ID == 1)
exact_spacing <- sqrt(sum((south_xy - north_xy)^2)) / 1852 / 2
exact <- generate_california_transects(north, south, 1, exact_spacing, land)
stopifnot(nrow(exact$endpoints) == 3,
          max(abs(as.numeric(sf::st_coordinates(exact$lines)[5, 1:2]) -
                    south_xy)) < 1e-5)

# A north-opening bay has three intersections; use the outer coast, not the bay.
bay_ring <- rbind(
  c(0, 150000), c(120000, 150000), c(120000, -50000), c(0, -50000),
  c(0, 50000), c(20000, 50000), c(20000, 90000), c(-10000, 90000),
  c(-10000, 30000), c(-20000, 30000), c(-20000, 120000),
  c(0, 120000), c(0, 150000)
)
bay_land <- sf::st_sfc(sf::st_polygon(list(bay_ring)), crs = 3310)
bay <- generate_california_transects(
  lonlat(c(0, 140000)), lonlat(c(0, 0)), 1, 70000 / 1852, bay_land
)
stopifnot(nrow(bay$endpoints) == 3,
          max(abs(as.numeric(sf::st_coordinates(bay$lines)[3, 1:2]) -
                    c(-20000, 70000))) < 1e-5,
          !any(bay$endpoints$crosses_land))

expect_error(generate_california_transects(north, south, 0, 2, land), "positive")
expect_error(generate_california_transects(north, south, 20, NA_real_, land),
             "positive")
expect_error(generate_california_transects(south, north, 20, 2, land),
             "north must be north")
expect_error(generate_california_transects(c(NA_real_, 35), south, 20, 2, land),
             "longitude, latitude")
expect_error(generate_california_transects(north, south, 20, 0.00001, land),
             "10,000")
expect_error(generate_california_transects(
  lonlat(north_xy + c(30000, 0)), south, 20, 2, land
), "more than 5000")
expect_error(generate_california_transects(
  north, south, 20, 2, sf::st_set_crs(land, NA)
), "defined CRS")
expect_error(generate_california_transects(
  north, south, 20, 2,
  sf::st_sfc(sf::st_linestring(rbind(north_xy, south_xy)), crs = 3310)
), "POLYGON")
expect_error(transects_main(c("--length-nm=20")), "missing arguments")
expect_error(transects_main(c(
  "--north-lon=wrong", "--north-lat=33.73",
  "--south-lon=-117.24", "--south-lat=32.67",
  "--length-nm=20", "--spacing-nm=2"
)), "must be numeric")

# Exercise the command interface, local land loading, and the actual CSV schema.
land_file <- tempfile(fileext = ".gpkg")
csv_file <- tempfile(fileext = ".csv")
sf::st_write(sf::st_sf(geometry = land), land_file, quiet = TRUE)
transects_main(c(
  paste0("--north-lon=", format(north[1], digits = 16)),
  paste0("--north-lat=", format(north[2], digits = 16)),
  paste0("--south-lon=", format(south[1], digits = 16)),
  paste0("--south-lat=", format(south[2], digits = 16)),
  "--length-nm=20", "--spacing-nm=2",
  paste0("--land-file=", land_file), paste0("--output=", csv_file)
))
csv <- read.csv(csv_file)
stopifnot(identical(names(csv), names(design$endpoints)),
          nrow(csv) == nrow(design$endpoints),
          max(abs(csv$length_nm - 20)) < 1e-7,
          max(abs(csv$shore_lon - design$endpoints$shore_lon)) < 1e-10)
unlink(c(land_file, csv_file))
cat("California transect regression checks passed.\n")
