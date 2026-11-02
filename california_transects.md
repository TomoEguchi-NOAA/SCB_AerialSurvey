# California coastal survey track lines

`california_transects.R` is a standalone alternative to
`draw_transect_lines.Rmd`; the original notebook is unchanged. It works along
the California mainland coast, including study areas spanning UTM zones.

## Install and run

Install R (4.0 or later) and the two packages once:

```r
install.packages(c("sf", "geosphere"))
```

From the repository directory:

```sh
Rscript california_transects.R \
  --north-lon=-118.10 --north-lat=33.73 \
  --south-lon=-117.24 --south-lat=32.67 \
  --length-nm=20 --spacing-nm=2 \
  --output=tracklines.csv
```

Use decimal degrees (WGS84), negative longitude for west, and coastal coordinates
for both study endpoints. `--length-nm` is **X**, the distance from the mainland
shoreline to each offshore endpoint, not an offset from shore. `--spacing-nm` is
**Y**, the separation between parallel lines. One nautical mile is 1,852 metres.
Run with `--help` to display the command syntax.

The first run downloads public-domain [Natural Earth 1:10m land polygons](https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-land/)
from a pinned revision of its GitHub mirror and caches them in
`tools::R_user_dir("california_transects", "cache")`. Subsequent runs work offline.
For higher-resolution work or offline first use, pass `--land-file=/path/to/land.gpkg`
(GeoJSON and shapefiles also work). Supply **land polygons**, including the
California mainland and relevant islands, with a defined coordinate reference
system. A coastline-only LINESTRING file is not supported.

## Geometry and interpretation

- Both study endpoints are snapped to the mainland shoreline (maximum 5 km).
  The largest connected land polygon in the California region is treated as
  mainland. Custom data must include that mainland, not just islands.
- The straight line between the snapped endpoints defines the **average
  coastline direction**, not a regression through local coastline vertices.
  The shorter mainland boundary path between them defines the study coastline.
- All track lines are parallel and perpendicular to that average direction in
  California Albers (EPSG:3310). Their spacing is exactly Y nautical miles
  **in that projected frame**, not distance along the winding shoreline.
  Ground spacing and bearings vary slightly because of projection distortion.
- Stations begin at the northern boundary and advance southward by Y. A final
  partial interval is omitted, so a line at the southern boundary occurs only
  if the baseline length is a multiple of Y. No lines are generated beyond it.
  For safety, requests exceeding 10,000 lines are rejected.
- Where a station intersects the study coastline more than once, the
  ocean-facing intersection is used, avoiding inland bay crossings. The
  offshore direction is the west-facing perpendicular to the north-to-south
  baseline; this is a mainland California design, not an island survey tool.
- Every line starts on the supplied shoreline and extends X nautical miles
  offshore. Length is solved using WGS84 ellipsoidal endpoint distance while
  retaining parallel lines in the projected frame. Projection-defined straight
  tracks are not constant-bearing rhumb lines or geodesic paths.
- Lines crossing more than one metre of supplied land, including islands or
  another part of the mainland, are **flagged and warned about, not shortened**.
  Review these lines on a map and revise the study extent or length as needed.
  Endpoints alone do not describe an obstacle-avoiding flight route.

The default coastline is generalized cartographic data, not a tidal shoreline
or navigation chart. A shoreline start has no operational clearance from land.
Always review a design against current, appropriately detailed local data before
survey planning or flight use. Small-scale coastal features absent from the
input data cannot be detected.

## CSV output

One row per line, ordered north to south:

| Column | Meaning |
|---|---|
| `ID` | Sequential line identifier |
| `shore_lat`, `shore_lon` | Mainland endpoint, WGS84 decimal degrees |
| `offshore_lat`, `offshore_lon` | Offshore endpoint, WGS84 decimal degrees |
| `length_nm` | Ellipsoidal distance between endpoints, nautical miles |
| `baseline_distance_nm` | Station distance along the average coastline baseline |
| `crosses_land` | Whether the line intersects supplied land beyond a shoreline touch |

Coordinate values are not rounded to display precision.

## Use from R and inspect a map

```r
source("california_transects.R")
design <- generate_california_transects(
  north = c(-118.10, 33.73), # longitude, latitude
  south = c(-117.24, 32.67),
  length_nm = 20,
  spacing_nm = 2
)
write.csv(design$endpoints, "tracklines.csv", row.names = FALSE)

plot(design$lines$geometry, col = "blue", axes = TRUE)
plot(design$land, col = "grey80", border = "grey50", add = TRUE)
plot(design$lines$geometry, col = "blue", add = TRUE)
plot(design$baseline, col = "red", lty = 2, add = TRUE)
```

All returned map geometries use EPSG:3310; the CSV uses WGS84. To use your own
land polygons from R, supply `land = sf::st_read("land.gpkg")`. The function
argument `snap_tolerance_m` changes the default 5,000 m snapping limit.

Run the standalone regression checks with:

```sh
Rscript tests/test_california_transects.R
```
