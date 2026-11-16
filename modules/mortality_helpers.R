# =============================================================================
# mortality_helpers.R
# Functions for SCF mortality parameter calibration.
#
# Covers:
#   Site mortality  (Eq 7): dNBR ~ Clay + ET + EWS + CWD + FineFuels + LadderFuels
#                           Gamma GLM with inverse link
#   Cohort mortality (Eq 8/9): P(mort) ~ BarkThickness + dNBR
#                              Binomial GLM with logit link
#
# Note: compute_effective_wind() is defined in spread_helpers.R and must be
# sourced before this file.  Do NOT redefine it here.
# =============================================================================

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, ggplot2, tidyr, purrr, stringr,
               lubridate, glue, broom, scales, httr, ncdf4, jsonlite)


# =============================================================================
# Clay raster helpers
# =============================================================================

# -----------------------------------------------------------------------------
# process_clay_raster
#
# Purpose: Build a clay fraction raster (0-1) over cal_vect by combining
#   SOLUS100 (primary) with SoilGrids (gap-fill).
#
# Args:
#   solus_clay_dir  path to SOLUS100 clay tiles; expects claytotal_{depth}_p.tif
#   sg_clay_dir     path to SoilGrids clay tiles; expects clay_{depth}.tif
#   cal_vect        SpatVector of calibration boundary, projected to working_crs
#   working_crs     CRS string for output
#   progress_fn     optional function(msg) for status messages
#
# Returns: SpatRaster of clay fraction (0-1) over cal_vect extent.
# -----------------------------------------------------------------------------

process_clay_raster <- function(solus_clay_dir, sg_clay_dir,
                                 cal_vect, working_crs,
                                 progress_fn = NULL) {

  depths <- c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm")

  # ---- Helper: load a set of depth files and average -------------------------
  load_depth_stack <- function(dir_path, pattern_fn, unit_divisor, label) {
    files <- vapply(depths, function(d) {
      file.path(dir_path, pattern_fn(d))
    }, character(1))

    present <- file.exists(files)
    if (!any(present))
      stop(sprintf(
        "%s directory '%s' has no matching files. Expected names like: %s",
        label, dir_path, basename(files[1])
      ))
    if (!all(present))
      warning(sprintf(
        "%s: %d of %d depth files found (missing: %s). Mean computed from present files.",
        label, sum(present), length(depths),
        paste(depths[!present], collapse = ", ")
      ))

    stack <- rast(files[present])
    mean_r <- app(stack, mean, na.rm = TRUE)
    mean_r / unit_divisor
  }

  # ---- SOLUS100 --------------------------------------------------------------
  if (!is.null(progress_fn)) progress_fn("Loading SOLUS100 clay layers...")
  message("process_clay_raster: loading SOLUS100 from ", solus_clay_dir)

  solus_r <- load_depth_stack(
    solus_clay_dir,
    pattern_fn    = function(d) sprintf("claytotal_%s_p.tif", d),
    unit_divisor  = 100,   # g/100g -> fraction 0-1
    label         = "SOLUS100"
  )

  # ---- Project and crop to cal_vect ------------------------------------------
  if (!is.null(progress_fn)) progress_fn("Projecting SOLUS100 to working CRS and cropping...")
  solus_proj <- project(solus_r, working_crs, method = "bilinear")
  solus_cal  <- crop(solus_proj, cal_vect, mask = TRUE)
  message(sprintf("SOLUS100 clay: extent [%s], range [%.3f, %.3f]",
                  paste(as.vector(ext(solus_cal)), collapse = " "),
                  min(values(solus_cal, mat = FALSE), na.rm = TRUE),
                  max(values(solus_cal, mat = FALSE), na.rm = TRUE)))

  # ---- SoilGrids -------------------------------------------------------------
  if (!is.null(progress_fn)) progress_fn("Loading SoilGrids clay layers for gap-fill...")
  message("process_clay_raster: loading SoilGrids from ", sg_clay_dir)

  sg_r <- load_depth_stack(
    sg_clay_dir,
    pattern_fn    = function(d) sprintf("clay_%s.tif", d),
    unit_divisor  = 1000,  # g/kg -> fraction 0-1
    label         = "SoilGrids"
  )

  sg_proj <- project(sg_r, working_crs, method = "bilinear")
  # Resample SoilGrids to SOLUS grid before gap-fill
  sg_cal  <- resample(sg_proj, solus_cal, method = "bilinear")
  sg_cal  <- crop(sg_cal, cal_vect, mask = TRUE)

  # ---- Gap-fill: replace NA SOLUS pixels with SoilGrids ---------------------
  if (!is.null(progress_fn)) progress_fn("Gap-filling SOLUS NA pixels with SoilGrids...")
  na_before <- sum(is.na(values(solus_cal, mat = FALSE)))
  solus_cal[is.na(solus_cal)] <- sg_cal
  na_after  <- sum(is.na(values(solus_cal, mat = FALSE)))
  message(sprintf("Clay gap-fill: %d NA pixels filled (%.1f%%); %d remain NA",
                  na_before - na_after,
                  100 * (na_before - na_after) / max(ncell(solus_cal), 1),
                  na_after))

  # ---- Clamp to valid range --------------------------------------------------
  solus_cal[solus_cal < 0] <- 0
  solus_cal[solus_cal > 1] <- 1

  names(solus_cal) <- "clay_fraction"
  solus_cal
}


# =============================================================================
# MTBS dNBR pixel loading
# =============================================================================

# -----------------------------------------------------------------------------
# load_mtbs_dnbr_pixels
#
# Purpose: Scan an MTBS dNBR directory, load per-fire GeoTIFFs, sample pixels,
#   and return a tidy tibble for mortality model fitting.
#
# Args:
#   mtbs_dnbr_dir  directory containing *_dnbr.tif files (not *_dnbr6.tif)
#   cal_vect       SpatVector of calibration boundary (projected to working_crs)
#   working_crs    CRS string
#   sample_frac    fraction of pixels to keep per fire (default 0.05)
#   min_dnbr       minimum dNBR to retain (filters non-fire/pre-fire; default 100)
#   max_dnbr       maximum dNBR (clips sensor saturation; default 2000)
#   progress_fn    optional function(msg)
#
# Returns: tibble(fire_id, fire_year, dNBR, x, y)
# -----------------------------------------------------------------------------

load_mtbs_dnbr_pixels <- function(mtbs_dnbr_dir, cal_vect, working_crs,
                                   sample_frac = 0.05,
                                   min_dnbr = 1, max_dnbr = 2000,
                                   dnbr_offset = 0,
                                   progress_fn = NULL) {

  # min_dnbr is a LEFT TRUNCATION of the Gamma response.  It must stay low.
  # The old default of 100 discarded the entire low-severity tail, which inflates
  # the fitted intercept (and therefore the baseline dNBR SCF predicts) and biases
  # every slope.  Gamma only needs strictly positive support, so 1 is enough.
  if (min_dnbr > 50)
    warning(sprintf(paste0(
      "min_dnbr = %g truncates the low-severity tail of the response. ",
      "Gamma only requires dNBR > 0; values above ~50 bias SiteMortalityB0 ",
      "upward and flatten the fitted slopes."), min_dnbr))

  # ---- Discover files --------------------------------------------------------
  all_files <- list.files(mtbs_dnbr_dir, pattern = "_dnbr\\.tif$",
                           full.names = TRUE, recursive = FALSE)
  # Exclude 6-class severity products
  all_files <- all_files[!grepl("_dnbr6\\.tif$", all_files)]

  if (length(all_files) == 0)
    stop(sprintf(paste0(
      "No *_dnbr.tif files found in '%s'. Check directory path and that files ",
      "are not named *_dnbr6.tif."
    ), mtbs_dnbr_dir))

  cat(sprintf("load_mtbs_dnbr_pixels: found %d dNBR files in %s\n",
              length(all_files), mtbs_dnbr_dir))

  # ---- WGS84 version of cal_vect for overlap check ---------------------------
  cal_wgs84 <- project(cal_vect, "EPSG:4326")
  cat(sprintf("  cal_vect WGS84 extent: xmin=%.3f xmax=%.3f ymin=%.3f ymax=%.3f\n",
              xmin(ext(cal_wgs84)), xmax(ext(cal_wgs84)),
              ymin(ext(cal_wgs84)), ymax(ext(cal_wgs84))))

  # ---- Loop over files -------------------------------------------------------
  pixel_list <- vector("list", length(all_files))
  n_loaded    <- 0L
  n_thematic  <- 0L

  for (i in seq_along(all_files)) {
    fpath <- all_files[i]
    fname <- basename(fpath)

    if (!is.null(progress_fn))
      progress_fn(sprintf("Processing file %d of %d ...", i, length(all_files)))

    cat(sprintf("  [%d/%d] %s\n", i, length(all_files), fname))

    # -- Extract fire_id and fire_year from filename ---------------------------
    fire_id     <- sub("_dnbr\\.tif$", "", fname)
    date_blocks <- regmatches(fname, gregexpr("[0-9]{8}", fname))[[1]]
    if (length(date_blocks) == 0) {
      cat(sprintf("    SKIP: cannot parse date from filename\n"))
      next
    }
    fire_year <- as.integer(substr(date_blocks[1], 1, 4))

    # First 8-digit block is the ignition date for per-fire MTBS bundles, but is
    # only the mosaic window start for the CONUS annual mosaics.  Keep it either
    # way; extract_site_predictors() decides whether it is usable as a burn date.
    fire_date <- suppressWarnings(as.Date(date_blocks[1], format = "%Y%m%d"))
    cat(sprintf("    fire_year=%d  fire_date=%s\n", fire_year,
                if (is.na(fire_date)) "NA" else format(fire_date)))

    # -- Load raster -----------------------------------------------------------
    r <- tryCatch(rast(fpath), error = function(e) {
      cat(sprintf("    SKIP: cannot read raster — %s\n", conditionMessage(e)))
      NULL
    })
    if (is.null(r)) next
    cat(sprintf("    raster CRS: %s  extent: %.3f %.3f %.3f %.3f\n",
                crs(r, describe = TRUE)$code,
                xmin(r), xmax(r), ymin(r), ymax(r)))

    # -- Overlap check using WGS84 extents ------------------------------------
    r_wgs84 <- tryCatch(
      project(r, "EPSG:4326", method = "near"),
      error = function(e) {
        cat(sprintf("    SKIP: project to WGS84 failed — %s\n", conditionMessage(e)))
        NULL
      }
    )
    if (is.null(r_wgs84)) next

    # Direct bounding-box overlap test — avoids terra::intersect() quirks
    re <- ext(r_wgs84)
    ce <- ext(cal_wgs84)
    overlap <- (xmax(re) > xmin(ce)) && (xmin(re) < xmax(ce)) &&
               (ymax(re) > ymin(ce)) && (ymin(re) < ymax(ce))

    if (!overlap) {
      cat(sprintf("    SKIP: no overlap with cal_vect (raster WGS84: %.3f %.3f %.3f %.3f)\n",
                  xmin(re), xmax(re), ymin(re), ymax(re)))
      next
    }

    # -- Diagnostics: raw values before any crop ------------------------------
    raw_vals  <- values(r, mat = FALSE)
    raw_nonNA <- sum(!is.na(raw_vals))
    raw_unique <- sort(unique(raw_vals[!is.na(raw_vals)]))
    raw_unique_show <- if (length(raw_unique) > 10)
      c(head(raw_unique, 5), "...", tail(raw_unique, 5)) else raw_unique
    message(sprintf("  %s: CRS=%s, raw non-NA=%d, unique values: %s",
                    fname, crs(r, describe = TRUE)$code,
                    raw_nonNA, paste(raw_unique_show, collapse = " ")))

    # -- Auto-reclassify raw MTBS severity classes (0-6) to dNBR midpoints ---
    # Files downloaded before reclassification was applied (or cached copies)
    # contain raw severity classes 1-6 rather than dNBR midpoints 185/355/600.
    # Detect by checking whether ALL non-NA values are <= 6.
    # Use subst() for exact integer matching — classify() interval approach
    # silently produces empty intervals when from == to (e.g. (2,2] = {}).
    is_thematic <- raw_nonNA > 0 && max(raw_unique, na.rm = TRUE) <= 6
    if (is_thematic) {
      cat(sprintf("    THEMATIC product detected - reclassifying classes to dNBR midpoints\n"))
      # Class 1 ("Unburned to Low") is now mapped to ~60 rather than dropped.
      # SCF assigns a dNBR to EVERY cell inside the perimeter, including barely
      # burned ones, so deleting class 1 made the calibration sample
      # systematically more severe than what SCF has to reproduce.
      # 5 = increased greenness and 6 = non-mapping remain NA.
      r <- subst(r,
                 from = c(0,  1,   2,   3,   4,   5,  6),
                 to   = c(NA, 60, 185, 355, 600,  NA, NA))
      n_thematic <- n_thematic + 1L
    }

    # -- Project, crop, mask to cal_vect --------------------------------------
    # Use nearest-neighbor: values are discrete thematic (185/355/600/NA);
    # bilinear would interpolate fire pixels into NA background and lose them.
    r_proj <- tryCatch({
      rp <- project(r, working_crs, method = "near")
      crop(rp, cal_vect, mask = TRUE)
    }, error = function(e) {
      message(sprintf("  Crop/mask failed for %s: %s", fname, conditionMessage(e)))
      NULL
    })
    if (is.null(r_proj)) next

    # -- Extract pixel coordinates and values ---------------------------------
    vals <- values(r_proj, mat = FALSE)
    xy   <- as.data.frame(xyFromCell(r_proj, seq_len(ncell(r_proj))))

    n_total    <- sum(!is.na(vals))
    n_in_range <- sum(!is.na(vals) & vals >= min_dnbr & vals <= max_dnbr)
    message(sprintf("  %s: after crop — %d non-NA pixels, %d in dNBR range [%d,%d]",
                    fname, n_total, n_in_range, min_dnbr, max_dnbr))

    df <- tibble(
      fire_id   = fire_id,
      fire_year = fire_year,
      fire_date = fire_date,
      # MTBS continuous products carry a per-fire dNBR offset in their metadata
      # that should be removed before use.  Zero by default; set dnbr_offset when
      # working from per-fire bundles whose offset you have read.
      dNBR      = vals - dnbr_offset,
      x         = xy[[1]],
      y         = xy[[2]]
    ) %>%
      filter(is.finite(dNBR), dNBR >= min_dnbr, dNBR <= max_dnbr)

    if (nrow(df) == 0) next

    # -- Sample -----------------------------------------------------------------
    n_keep <- max(1L, as.integer(ceiling(nrow(df) * sample_frac)))
    df <- df[sample.int(nrow(df), min(n_keep, nrow(df))), ]

    pixel_list[[i]] <- df
    n_loaded <- n_loaded + 1L
  }

  result <- bind_rows(pixel_list)
  if (nrow(result) == 0)
    stop(sprintf(paste0(
      "No MTBS dNBR pixels survived filtering (min_dnbr=%d, max_dnbr=%d, ",
      "sample_frac=%.2f). Check that cal_vect overlaps the MTBS fire footprints."
    ), min_dnbr, max_dnbr, sample_frac))

  message(sprintf("load_mtbs_dnbr_pixels: %d fires loaded, %s pixels retained",
                  n_loaded, scales::comma(nrow(result))))

  n_lvl <- length(unique(result$dNBR))
  attr(result, "thematic_files") <- n_thematic
  attr(result, "n_dnbr_levels")  <- n_lvl

  if (n_thematic > 0) {
    warning(sprintf(paste0(
      "%d of %d MTBS files were THEMATIC (6-class) products, so the response ",
      "takes only %d distinct values. A Gamma GLM fitted to a handful of class ",
      "midpoints has a meaningless variance structure, compressed predicted ",
      "range, and coefficients driven by the class breaks rather than by real ",
      "severity gradients. Use continuous per-fire MTBS rasters (RdNBR preferred, ",
      "as in Robbins' reference calibration) for anything beyond a smoke test."),
      n_thematic, n_loaded, n_lvl))
  } else if (n_lvl < 25) {
    warning(sprintf(paste0(
      "Only %d distinct dNBR values across %s pixels - the response is close to ",
      "categorical and the Gamma fit will be poorly determined."),
      n_lvl, scales::comma(nrow(result))))
  }

  result
}


# =============================================================================
# MTBS burn severity download via IIPP ImageServer
# =============================================================================

# -----------------------------------------------------------------------------
# fetch_mtbs_from_imageserver
#
# Purpose: Download annual MTBS thematic burn severity rasters for the
#   calibration area from the interagency IIPP ImageServer, reclassify the
#   6-class severity product to approximate dNBR midpoints, and write one
#   GeoTIFF per year to out_dir in a format compatible with
#   load_mtbs_dnbr_pixels().
#
# Data source:
#   IIPP ArcGIS ImageServer (CONUS annual mosaics, 30m, 1984-present):
#   https://imagery.geoplatform.gov/iipp/rest/services/Fire_Aviation/
#     USFS_EDW_MTBS_CONUS/ImageServer
#
# ImageServer structure:
#   One raster per year (e.g. "mtbs_CONUS_2000") covering all CONUS fires in
#   that year. The mosaic stores thematic burn severity classes 0-6:
#     0 = background / no fire mapped
#     1 = unburned to low
#     2 = low        → approximate dNBR midpoint: 185
#     3 = moderate   → approximate dNBR midpoint: 355
#     4 = high       → approximate dNBR midpoint: 600
#     5 = increased greenness (excluded from mortality model)
#     6 = non-mapping mask (excluded)
#   Classes 0, 1, 5, 6 are set to NA; classes 2-4 are mapped to dNBR
#   midpoints so that load_mtbs_dnbr_pixels() min_dnbr filter still works.
#
# Output file naming convention compatible with load_mtbs_dnbr_pixels():
#   mtbs_CONUS_{YYYYMMDD}_{YYYYMMDD}_dnbr.tif
#   The first 8-digit block is Jan 1 of the year; load_mtbs_dnbr_pixels()
#   parses fire_year from it.
#
# Args:
#   cal_vect    SpatVector of calibration boundary (any CRS)
#   working_crs CRS string (e.g. "EPSG:32610")
#   year_start  integer start year
#   year_end    integer end year
#   out_dir     output directory for annual dNBR GeoTIFFs (created if needed)
#   max_px      maximum pixels per side of exported image (default 4096);
#               lower = faster download, lower spatial resolution
#   progress_fn optional function(msg)
#
# Returns: list(summary = tibble(year, objectid, status, out_path),
#               n_ok, n_failed)
# -----------------------------------------------------------------------------

fetch_mtbs_from_imageserver <- function(cal_vect, working_crs,
                                         year_start, year_end, out_dir,
                                         max_px      = 4096L,
                                         progress_fn = NULL) {

  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("jsonlite package required: install.packages('jsonlite')")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  base_url <- paste0(
    "https://imagery.geoplatform.gov/iipp/rest/services/",
    "Fire_Aviation/USFS_EDW_MTBS_CONUS/ImageServer"
  )

  # ---- Calibration boundary in WGS84 ----------------------------------------
  cal_wgs84 <- project(cal_vect, "EPSG:4326")
  bb <- as.vector(ext(cal_wgs84))   # xmin, xmax, ymin, ymax

  # ---- Step 1: Query catalog for objectids in year range ---------------------
  if (!is.null(progress_fn))
    progress_fn("Querying MTBS ImageServer catalog...")

  where_clause <- sprintf("Year >= %d AND Year <= %d", year_start, year_end)
  geom_json    <- sprintf(
    '{"xmin":%f,"ymin":%f,"xmax":%f,"ymax":%f,"spatialReference":{"wkid":4326}}',
    bb[1], bb[3], bb[2], bb[4]
  )

  query_resp <- tryCatch(
    httr::GET(
      paste0(base_url, "/query"),
      query = list(
        where        = where_clause,
        geometry     = geom_json,
        geometryType = "esriGeometryEnvelope",
        inSR         = "4326",
        spatialRel   = "esriSpatialRelIntersects",
        outFields    = "objectid,name,Year,start_time,end_time",
        returnGeometry = "false",
        f            = "json"
      ),
      httr::timeout(60)
    ),
    error = function(e) stop("ImageServer catalog query failed: ", conditionMessage(e))
  )

  if (httr::http_error(query_resp))
    stop("ImageServer query returned HTTP ", httr::status_code(query_resp))

  catalog  <- jsonlite::fromJSON(httr::content(query_resp, "text", encoding = "UTF-8"),
                                  simplifyDataFrame = TRUE)

  if (is.null(catalog$features) || nrow(catalog$features) == 0)
    stop(sprintf(
      "No MTBS annual mosaics found for years %d\u2013%d in the ImageServer.",
      year_start, year_end
    ))

  attrs <- catalog$features$attributes
  message(sprintf("fetch_mtbs_from_imageserver: found %d annual mosaics", nrow(attrs)))

  # ---- Step 2: Compute export size from cal_vect extent at 30m ---------------
  # Project bbox to metres for resolution calculation
  cal_meter <- project(cal_wgs84, "EPSG:3857")
  bb_m      <- as.vector(ext(cal_meter))
  width_m   <- bb_m[2] - bb_m[1]
  height_m  <- bb_m[4] - bb_m[3]
  px_w      <- min(max_px, as.integer(ceiling(width_m  / 30)))
  px_h      <- min(max_px, as.integer(ceiling(height_m / 30)))
  px_w      <- max(px_w, 128L)
  px_h      <- max(px_h, 128L)

  actual_res_m <- max(width_m / px_w, height_m / px_h)
  message(sprintf("  Export size: %d x %d px (~%.0f m/pixel)",
                  px_w, px_h, actual_res_m))

  # ---- Severity class -> approximate dNBR mapping ----------------------------
  # Classes 0,1,5,6 -> NA (no fire or non-burn)
  # Class 2 (low)      -> 185 dNBR
  # Class 3 (moderate) -> 355 dNBR
  # Class 4 (high)     -> 600 dNBR
  # subst() used (not classify()) — classify() interval matching silently
  # produces empty sets when from == to (e.g. the interval (2,2] = {}).
  reclass_from <- c(0,  1,   2,   3,   4,   5,  6)
  reclass_to   <- c(NA, NA, 185, 355, 600,  NA, NA)

  # ---- Step 3: Export + reclassify one raster per year ----------------------
  n_years <- nrow(attrs)
  summary_rows <- vector("list", n_years)

  for (i in seq_len(n_years)) {
    row       <- attrs[i, ]
    oid       <- as.integer(row$objectid)
    yr        <- as.integer(row$Year)
    start_str <- trimws(as.character(row$start_time))   # e.g. "20000101"
    end_str   <- trimws(as.character(row$end_time))

    # Build output filename parseable by load_mtbs_dnbr_pixels()
    if (nchar(start_str) < 8) start_str <- sprintf("%d0101", yr)
    if (nchar(end_str)   < 8) end_str   <- sprintf("%d1231", yr)
    out_path <- file.path(out_dir,
                           sprintf("mtbs_CONUS_%s_%s_dnbr.tif", start_str, end_str))

    status_i <- "pending"

    if (!is.null(progress_fn))
      progress_fn(sprintf("MTBS ImageServer: exporting year %d (%d/%d)...",
                          yr, i, n_years))
    message(sprintf("fetch_mtbs_from_imageserver: year %d (OID %d)", yr, oid))

    # Skip if already present
    if (file.exists(out_path) && file.info(out_path)$size > 1e4) {
      message(sprintf("  year %d: cached at %s", yr, basename(out_path)))
      status_i <- "cached"
      summary_rows[[i]] <- tibble(year=yr, objectid=oid, status=status_i,
                                   out_path=out_path)
      next
    }

    # Mosaic rule: lock to this year's raster
    mosaic_rule    <- sprintf(
      '{"mosaicMethod":"esriMosaicLockRaster","lockRasterIds":[%d]}', oid
    )
    rendering_rule <- '{"rasterFunction":"None"}'  # raw class values, no symbology

    # exportImage -> binary GeoTIFF
    tmp_tif <- tempfile(fileext = ".tif")
    export_resp <- tryCatch(
      httr::GET(
        paste0(base_url, "/exportImage"),
        query = list(
          bbox          = sprintf("%f,%f,%f,%f", bb[1], bb[3], bb[2], bb[4]),
          bboxSR        = "4326",
          size          = sprintf("%d,%d", px_w, px_h),
          imageSR       = "4326",
          format        = "tiff",
          pixelType     = "U8",
          noData        = "0",
          mosaicRule    = mosaic_rule,
          renderingRule = rendering_rule,
          f             = "image"
        ),
        httr::write_disk(tmp_tif, overwrite = TRUE),
        httr::timeout(300)
      ),
      error = function(e) {
        message(sprintf("  year %d: exportImage error — %s", yr, conditionMessage(e)))
        NULL
      }
    )

    if (is.null(export_resp) || httr::http_error(export_resp)) {
      message(sprintf("  year %d: HTTP %s",
                      yr, if (!is.null(export_resp)) httr::status_code(export_resp) else "error"))
      if (file.exists(tmp_tif)) file.remove(tmp_tif)
      status_i <- "failed"
      summary_rows[[i]] <- tibble(year=yr, objectid=oid, status=status_i,
                                   out_path=NA_character_)
      next
    }

    # Read, reclassify, project, crop, write
    # IMPORTANT: terra creates a lazy file reference; do NOT delete tmp_tif
    # until all terra operations are complete (round, classify, project, crop).
    r_raw <- tryCatch(rast(tmp_tif), error = function(e) {
      message(sprintf("  year %d: cannot read exported TIF — %s", yr, conditionMessage(e)))
      NULL
    })

    if (is.null(r_raw)) {
      if (file.exists(tmp_tif)) file.remove(tmp_tif)
      status_i <- "failed"
      summary_rows[[i]] <- tibble(year=yr, objectid=oid, status=status_i,
                                   out_path=NA_character_)
      next
    }

    # Classify: any value outside 0-6 -> NA first, then remap.
    # tmp_tif is still on disk here — terra needs it for lazy reads.
    r_int  <- round(r_raw)
    r_int[r_int < 0 | r_int > 6] <- NA
    r_dnbr <- subst(r_int, from = reclass_from, to = reclass_to)

    # Skip entirely-NA years (no fires in cal_vect)
    if (all(is.na(values(r_dnbr, mat = FALSE)))) {
      message(sprintf("  year %d: no burned pixels in cal_vect — skipping", yr))
      if (file.exists(tmp_tif)) file.remove(tmp_tif)
      status_i <- "no_fires"
      summary_rows[[i]] <- tibble(year=yr, objectid=oid, status=status_i,
                                   out_path=NA_character_)
      next
    }

    # Project to working_crs and crop to cal_vect
    r_proj   <- project(r_dnbr, working_crs, method = "near")
    cal_proj <- project(cal_wgs84, working_crs)
    r_crop   <- crop(r_proj, cal_proj, mask = TRUE)

    writeRaster(r_crop, out_path, overwrite = TRUE,
                datatype = "INT2S",    # signed 16-bit: covers 0-600 without loss
                gdal = c("COMPRESS=LZW", "PREDICTOR=2"))

    # Safe to remove temp file now that writeRaster has finished
    if (file.exists(tmp_tif)) file.remove(tmp_tif)

    message(sprintf("  year %d: written -> %s (%.1f MB)",
                    yr, basename(out_path), file.info(out_path)$size / 1e6))
    status_i <- "downloaded"
    summary_rows[[i]] <- tibble(year=yr, objectid=oid, status=status_i,
                                 out_path=out_path)
  }

  summary_tbl <- bind_rows(summary_rows)
  n_ok     <- sum(summary_tbl$status %in% c("downloaded", "cached"))
  n_failed <- sum(summary_tbl$status == "failed")

  message(sprintf(
    "fetch_mtbs_from_imageserver: %d years OK, %d failed, %d no fires in region",
    n_ok, n_failed,
    sum(summary_tbl$status == "no_fires")
  ))

  list(summary = summary_tbl, n_ok = n_ok, n_failed = n_failed)
}


# =============================================================================
# TerraClimate ET and CWD
# =============================================================================

# -----------------------------------------------------------------------------
# fetch_terraclimate_et_cwd
#
# Purpose: Fetch TerraClimate potential evapotranspiration (pet) and climatic
#   water deficit (def) by downloading per-year NetCDF files directly from the
#   TerraClimate THREDDS HTTP fileServer using httr.
#
# NOTE: SCF FireEvent.cs uses SiteVars.PotentialEvapotranspiration (PET) for
#   SiteMortalityB2, NOT actual ET (AET). PET is downloaded here to match.
#
# Why not OPeNDAP / vsicurl:
#   - The TerraClimate OPeNDAP server (dodsC) returns a malformed DAP2 DDS
#     response ("syntax error, unexpected $end") — ncdf4 and climateR both fail.
#   - terra's /vsicurl/ for NetCDF requires Linux userfaultfd and does not work
#     on macOS.
#   Direct HTTP download with httr is the only reliable cross-platform path.
#
# File URL pattern:
#   http://thredds.northwestknowledge.net:8080/thredds/fileServer/
#     TERRACLIMATE_ALL/data/TerraClimate_{var}_{year}.nc
#   ~50-100 MB per year per variable; cached locally to avoid re-downloading.
#
# Args:
#   cal_vect    SpatVector of calibration boundary (any projected CRS)
#   working_crs CRS string (e.g. "EPSG:32610")
#   years       integer vector — only these years are downloaded/processed
#   cache_dir   directory for cached NetCDF files; set a persistent path so
#               files survive across sessions (default: tempdir subdir)
#   progress_fn optional function(msg) for progress reporting
#
# Returns: list(et_r = SpatRaster, cwd_r = SpatRaster) with layers named
#   by year.  Returns NULL on complete failure.
# -----------------------------------------------------------------------------

fetch_terraclimate_et_cwd <- function(cal_vect, working_crs, years,
                                       cache_dir   = file.path(tempdir(),
                                                               "terraclimate_cache"),
                                       progress_fn = NULL) {

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  tc_base <- paste0(
    "http://thredds.northwestknowledge.net:8080/thredds/fileServer/",
    "TERRACLIMATE_ALL/data"
  )

  # ---- Helper: ensure local file exists, download if needed -----------------
  get_nc_path <- function(var_name, yr) {
    fname    <- sprintf("TerraClimate_%s_%d.nc", var_name, yr)
    local_nc <- file.path(cache_dir, fname)

    # Use cached file if it exists and is plausibly complete (> 1 MB)
    if (file.exists(local_nc) && file.info(local_nc)$size > 1e6) {
      message(sprintf("  %s: using cached file", fname))
      return(local_nc)
    }

    url <- paste0(tc_base, "/", fname)
    message(sprintf("  %s: downloading from THREDDS...", fname))

    # Write to a temp file first; rename on success to avoid corrupt partials
    tmp <- paste0(local_nc, ".part")
    resp <- tryCatch(
      httr::GET(
        url,
        httr::write_disk(tmp, overwrite = TRUE),
        httr::timeout(600)        # 10-min timeout; files are 50-100 MB
      ),
      error = function(e) {
        if (file.exists(tmp)) file.remove(tmp)
        warning(sprintf("Download error for %s: %s", fname, conditionMessage(e)))
        NULL
      }
    )

    if (is.null(resp)) return(NULL)

    if (httr::http_error(resp)) {
      if (file.exists(tmp)) file.remove(tmp)
      warning(sprintf("HTTP %d for %s", httr::status_code(resp), fname))
      return(NULL)
    }

    file.rename(tmp, local_nc)
    message(sprintf("  %s: download complete (%.1f MB)",
                    fname, file.info(local_nc)$size / 1e6))
    local_nc
  }

  # ---- Helper: select the target variable layers from a SpatRaster ----------
  # TerraClimate files are single-variable; terra may name layers
  # "pet_1".."pet_12" (or "aet_1".."aet_12") or just use numeric indices.
  select_var_layers <- function(rr, var_name) {
    nm  <- names(rr)
    idx <- grep(paste0("^", var_name, "(_[0-9]+)?$"), nm, ignore.case = TRUE)
    if (length(idx) >= 12) return(rr[[idx[seq_len(12)]]])
    if (nlyr(rr) >= 12)   return(rr[[seq_len(12)]])
    rr   # partial year at record edges
  }

  # ---- Process one variable across all requested years ----------------------
  process_var <- function(var_name, label) {
    annual_list <- vector("list", length(years))

    for (yi in seq_along(years)) {
      yr <- years[yi]
      if (!is.null(progress_fn))
        progress_fn(sprintf(
          "TerraClimate %s: year %d (%d of %d) — downloading if not cached...",
          label, yr, yi, length(years)
        ))
      message(sprintf("fetch_terraclimate_et_cwd: %s %d", var_name, yr))

      nc_path <- get_nc_path(var_name, yr)
      if (is.null(nc_path)) {
        warning(sprintf("TerraClimate %s %d: skipping (download failed).", var_name, yr))
        next
      }

      r12 <- tryCatch({
        rr <- rast(nc_path)
        select_var_layers(rr, var_name)
      }, error = function(e) {
        warning(sprintf("Cannot read %s %d: %s", var_name, yr, conditionMessage(e)))
        NULL
      })

      if (!is.null(r12))
        annual_list[[yi]] <- app(r12, sum, na.rm = TRUE)
    }

    valid <- !sapply(annual_list, is.null)
    if (!any(valid))
      stop(sprintf("No annual layers built for TerraClimate %s. Check internet connection and cache_dir.", var_name))

    s        <- rast(annual_list[valid])
    names(s) <- as.character(years[valid])
    s
  }

  # ---- Fetch PET (pet) and CWD (def) ----------------------------------------
  # SCF uses PotentialEvapotranspiration (PET) for SiteMortalityB2, not AET.
  et_raw  <- process_var("pet", "potential ET")
  cwd_raw <- process_var("def", "climatic water deficit")

  # ---- Reproject and crop to cal_vect ----------------------------------------
  if (!is.null(progress_fn)) progress_fn("Reprojecting TerraClimate to working CRS...")

  et_proj  <- project(et_raw,  working_crs, method = "bilinear")
  cwd_proj <- project(cwd_raw, working_crs, method = "bilinear")

  et_cal  <- crop(et_proj,  cal_vect, mask = TRUE)
  cwd_cal <- crop(cwd_proj, cal_vect, mask = TRUE)

  message(sprintf("fetch_terraclimate_et_cwd: ET layers=%d, CWD layers=%d",
                  nlyr(et_cal), nlyr(cwd_cal)))
  list(et_r = et_cal, cwd_r = cwd_cal)
}


# =============================================================================
# Annual fire-season wind
# =============================================================================

# -----------------------------------------------------------------------------
# compute_annual_fire_season_wind
#
# Purpose: Summarise rv$wind_daily to annual fire-season mean wind speed and
#   mean direction (circular mean).
#
# Args:
#   wind_daily    tibble(date, WindSpeed_kmh, WindDir_deg) from rv$wind_daily
#   fire_months   integer vector of months to include (default June–September)
#
# Returns: tibble(year, mean_ws_kmh, mean_wd_deg)
# -----------------------------------------------------------------------------

compute_annual_fire_season_wind <- function(wind_daily, fire_months = 6:9) {

  wind_daily %>%
    mutate(
      year  = lubridate::year(date),
      month = lubridate::month(date)
    ) %>%
    filter(month %in% fire_months, is.finite(WindSpeed_kmh)) %>%
    group_by(year) %>%
    summarise(
      mean_ws_kmh = mean(WindSpeed_kmh, na.rm = TRUE),
      # Circular mean: atan2(mean sin, mean cos) -> degrees
      mean_wd_deg = {
        wd_rad <- WindDir_deg * pi / 180
        atan2(mean(sin(wd_rad), na.rm = TRUE),
              mean(cos(wd_rad), na.rm = TRUE)) * 180 / pi
      },
      .groups = "drop"
    ) %>%
    mutate(mean_wd_deg = (mean_wd_deg + 360) %% 360)
}


# =============================================================================
# Ladder fuels (LANDFIRE CBH)
# =============================================================================

# -----------------------------------------------------------------------------
# fetch_landfire_ladder_fuels
#
# Purpose: Build a LADDER FUEL BIOMASS surface in g m-2, on the scale SCF
#   actually uses, from LANDFIRE canopy products.
#
# WHY BIOMASS AND NOT AN INDEX (FireEvent.cs:319-325):
#   SCF computes
#       ladderFuelBiomass = SUM over cohorts of cohort.Data.Biomass
#                           for species in LadderFuelSpeciesList
#                           with Age <= LadderFuelMaxAge
#   That is RAW summed cohort biomass in g m-2, typically 0 - 3000, and it enters
#   Equation 7 unbounded and unscaled.  A 0-1 index fitted here would produce a
#   B6 roughly 1000x too large: at 800 g m-2 the B6 term alone would drive the
#   linear predictor negative, hit the max(eta, 0.0005) clamp in FireEvent.cs:371
#   and pin dNBR to its 2000 ceiling across the landscape.  So the calibration
#   predictor must be biomass, in g m-2.
#
# CONSTRUCTION:
#   The ladder layer is the part of the canopy hanging below `ladder_ceiling_m`,
#   i.e. the sub-canopy stratum that carries surface fire into the crown - the
#   same stratum SCF represents with sub-LadderFuelMaxAge cohorts.
#
#       depth_m   = max(0, min(CH_m, ladder_ceiling_m) - CBH_m)
#       ladder_g  = CBD_kg_m3 * depth_m * 1000
#
#   Where CBH >= ladder_ceiling_m the canopy starts above the ladder zone and
#   ladder fuel is 0.  Where CBH is low and the canopy fills the zone, the full
#   load is counted.  Units: kg m-3 * m = kg m-2, * 1000 -> g m-2.
#
# LANDFIRE native units (LF2024):
#   CBH  metres * 10       CH  metres * 10       CBD  kg m-3 * 100
#   -9999 / negative = nodata; 0 = non-forested (no canopy, no ladder fuel)
#
# Args:
#   cal_vect          SpatVector of calibration boundary (any CRS)
#   working_crs       CRS string
#   spread_template   SpatRaster defining the output grid
#   ladder_ceiling_m  top of the ladder-fuel stratum in metres (default 8).
#                     Should correspond to the height of a cohort at
#                     LadderFuelMaxAge for your ladder species.
#   biomass_basis     "biomass" (default) returns g dry biomass m-2, matching
#                     NECN cohort biomass.  "carbon" halves it for landscapes
#                     whose succession extension reports g C m-2.
#   progress_fn       optional function(msg)
#
# Returns: SpatRaster on spread_template grid, values in g m-2 (>= 0).
# -----------------------------------------------------------------------------

# Internal: pull one LANDFIRE LF2024 CONUS product as a SpatRaster over a bbox.
.fetch_landfire_layer <- function(product, bb, progress_fn = NULL) {

  deg_per_30m <- 30 / 111320
  nx <- min(as.integer(ceiling((bb["xmax"] - bb["xmin"]) / deg_per_30m)), 4096L)
  ny <- min(as.integer(ceiling((bb["ymax"] - bb["ymin"]) / deg_per_30m)), 4096L)

  if (!is.null(progress_fn))
    progress_fn(sprintf("Downloading LANDFIRE %s...", product))
  message(sprintf(".fetch_landfire_layer: %s at %d x %d px", product, nx, ny))

  base_url <- sprintf(paste0("https://lfps.usgs.gov/arcgis/rest/services/",
                             "Landfire_LF2024/LF2024_%s_CONUS/ImageServer/exportImage"),
                      product)

  tmp_tif <- tempfile(fileext = ".tif")
  resp <- tryCatch(
    httr::GET(
      base_url,
      query = list(
        bbox          = sprintf("%.6f,%.6f,%.6f,%.6f",
                                bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"]),
        bboxSR        = "4326",
        imageSR       = "4326",
        size          = sprintf("%d,%d", nx, ny),
        format        = "tiff",
        pixelType     = "S16",
        renderingRule = '{"rasterFunction":"None"}',
        f             = "image"
      ),
      httr::write_disk(tmp_tif, overwrite = TRUE),
      httr::timeout(300)
    ),
    error = function(e)
      stop(sprintf("LANDFIRE %s HTTP request failed: %s", product, conditionMessage(e)))
  )

  if (httr::http_error(resp))
    stop(sprintf("LANDFIRE %s REST API returned HTTP %d.", product,
                 httr::status_code(resp)))

  r <- tryCatch(rast(tmp_tif), error = function(e)
    stop(sprintf("Could not read LANDFIRE %s response as raster: %s",
                 product, conditionMessage(e))))

  rv <- values(r, mat = FALSE)
  message(sprintf("  %s raw range: [%.0f, %.0f]",
                  product, min(rv, na.rm = TRUE), max(rv, na.rm = TRUE)))
  r
}


fetch_landfire_ladder_fuels <- function(cal_vect, working_crs, spread_template,
                                        ladder_ceiling_m = 8,
                                        biomass_basis    = c("biomass", "carbon"),
                                        progress_fn      = NULL) {

  biomass_basis <- match.arg(biomass_basis)

  if (!is.finite(ladder_ceiling_m) || ladder_ceiling_m <= 0)
    stop("ladder_ceiling_m must be a positive height in metres.")

  cal_wgs84 <- project(cal_vect, "EPSG:4326")
  bb        <- as.vector(ext(cal_wgs84))
  names(bb) <- c("xmin", "xmax", "ymin", "ymax")

  # ---- Pull the three canopy products ----------------------------------------
  cbh_r <- .fetch_landfire_layer("CBH", bb, progress_fn)
  cbd_r <- .fetch_landfire_layer("CBD", bb, progress_fn)
  ch_r  <- .fetch_landfire_layer("CH",  bb, progress_fn)

  # Align in case the service returns marginally different grids
  cbd_r <- resample(cbd_r, cbh_r, method = "bilinear")
  ch_r  <- resample(ch_r,  cbh_r, method = "bilinear")

  if (!is.null(progress_fn)) progress_fn("Computing ladder fuel biomass...")

  # ---- Native units -> SI ----------------------------------------------------
  cbh_m    <- cbh_r / 10.0      # decimetres -> m
  ch_m     <- ch_r  / 10.0      # decimetres -> m
  cbd_kgm3 <- cbd_r / 100.0     # kg m-3 * 100 -> kg m-3

  # Nodata (-9999 and friends). 0 is legitimate (non-forested) and is kept.
  cbh_m[cbh_m    < 0] <- NA
  ch_m[ch_m      < 0] <- NA
  cbd_kgm3[cbd_kgm3 < 0] <- NA

  # ---- Ladder-zone depth and load -------------------------------------------
  top      <- min(ch_m, ladder_ceiling_m)          # terra: cellwise min vs scalar
  depth_m  <- top - cbh_m
  depth_m[depth_m < 0] <- 0                        # CBH above the ladder ceiling

  ladder_g <- cbd_kgm3 * depth_m * 1000            # kg m-2 -> g m-2
  if (biomass_basis == "carbon") ladder_g <- ladder_g * 0.5
  ladder_g[ladder_g < 0] <- 0

  lv <- values(ladder_g, mat = FALSE)
  lv <- lv[is.finite(lv)]
  if (length(lv) == 0)
    stop("Ladder fuel biomass calculation produced no valid values. ",
         "Check that the calibration boundary overlaps CONUS LANDFIRE coverage.")

  message(sprintf(paste0("fetch_landfire_ladder_fuels: ladder fuel biomass ",
                         "range [%.0f, %.0f] g m-2, mean %.0f, %.0f%% non-zero ",
                         "(ceiling %.1f m, basis '%s')"),
                  min(lv), max(lv), mean(lv), 100 * mean(lv > 0),
                  ladder_ceiling_m, biomass_basis))

  if (max(lv) < 1)
    warning("Ladder fuel biomass is essentially zero everywhere. ",
            "Either the landscape is non-forested or ladder_ceiling_m (",
            ladder_ceiling_m, " m) is below the local canopy base height. ",
            "SiteMortalityB6 will be unidentifiable — consider setting it to 0.")

  # ---- Reproject and resample to the working grid ----------------------------
  if (!is.null(progress_fn))
    progress_fn("Reprojecting ladder fuels to working CRS...")
  lf_proj <- project(ladder_g, working_crs, method = "bilinear")
  out     <- resample(lf_proj, spread_template, method = "bilinear")
  names(out) <- "ladder_fuel_g_m2"
  out
}



# =============================================================================
# Site predictor extraction
# =============================================================================

# -----------------------------------------------------------------------------
# extract_site_predictors
#
# Purpose: Join spatial predictor rasters to the MTBS pixel tibble by
#   extracting raster values at pixel coordinates.
#
# Args:
#   pixel_df        tibble(fire_id, fire_year, dNBR, x, y) from load_mtbs_dnbr_pixels
#   clay_r          SpatRaster of clay fraction (0-1) over cal_vect
#   et_r            multi-layer SpatRaster of annual ET (mm/yr); layer names = year
#   cwd_r           multi-layer SpatRaster of annual CWD (mm/yr); layer names = year
#   wind_annual_tbl tibble(year, mean_ws_kmh, mean_wd_deg) or NULL
#   slope_r         SpatRaster of slope (degrees) from rv$slope_r (may be small extent)
#   aspect_r        SpatRaster of aspect (degrees) from rv$aspect_r (may be small extent)
#   fine_fuels_r    SpatRaster of fine fuel index (0-1) or NULL
#   ladder_fuels_r  SpatRaster of ladder fuel BIOMASS in g m-2 (NOT a 0-1
#                   index) from fetch_landfire_ladder_fuels(), or NULL
#   working_crs     CRS string (pixel_df x,y are in this CRS)
#
# Returns: pixel_df with additional columns Clay, ET_mm, CWD_mm, EWS_kmh,
#   FineFuels, LadderFuels; rows with any NA predictor dropped.
#   Predictors that are entirely unavailable are left NA rather than filled with
#   a constant, so fit_site_mortality() drops them and writes 0.0 into the slot.
# -----------------------------------------------------------------------------

extract_site_predictors <- function(pixel_df, clay_r, et_r, cwd_r,
                                     wind_annual_tbl, slope_r, aspect_r,
                                     fine_fuels_r, ladder_fuels_r,
                                     working_crs = NULL,
                                     wind_daily  = NULL) {

  # ---- Build SpatVector of pixel locations -----------------------------------
  pv <- vect(pixel_df[, c("x", "y")], geom = c("x", "y"),
             crs = if (!is.null(working_crs)) working_crs else crs(clay_r))

  # ---- Clay ------------------------------------------------------------------
  message("extract_site_predictors: extracting Clay...")
  clay_vals <- terra::extract(clay_r, pv)[, 2]
  pixel_df$Clay <- clay_vals

  # ---- ET (previous water year) — vectorised extract -------------------------
  # Strategy: extract all year-layers at once, then row-index into the result
  # matrix using fire_year - 1.  This is O(1) terra::extract calls vs O(n_pixels).
  extract_annual_by_year <- function(r_annual, prev_years, label) {
    year_names <- names(r_annual)
    # Extract all layers for all pixels in one call -> matrix [n_pixels x n_years]
    ex_mat <- terra::extract(r_annual, pv)    # first col = ID, rest = year layers
    ex_mat <- ex_mat[, -1, drop = FALSE]      # drop ID column
    colnames(ex_mat) <- year_names

    vals <- vapply(seq_len(nrow(pixel_df)), function(i) {
      yr  <- as.character(prev_years[i])
      col <- which(colnames(ex_mat) == yr)
      if (length(col) == 0) NA_real_ else ex_mat[i, col]
    }, numeric(1))
    message(sprintf("extract_site_predictors: %s extracted, %.1f%% non-NA",
                    label, 100 * mean(is.finite(vals))))
    vals
  }

  if (!is.null(et_r)) {
    pixel_df$ET_mm <- extract_annual_by_year(et_r,  pixel_df$fire_year - 1, "ET")
  } else {
    message("extract_site_predictors: ET raster not available, using NA")
    pixel_df$ET_mm <- NA_real_
  }

  # ---- CWD (previous water year) — vectorised extract ----------------------
  if (!is.null(cwd_r)) {
    pixel_df$CWD_mm <- extract_annual_by_year(cwd_r, pixel_df$fire_year - 1, "CWD")
  } else {
    message("extract_site_predictors: CWD raster not available, using NA")
    pixel_df$CWD_mm <- NA_real_
  }

  # ---- Effective Wind Speed --------------------------------------------------
  # SCF evaluates EWS with the wind of the DAY THE CELL BURNS
  # (FireEvent.cs:604-606 reads DailyWindSpeed[day]).  Calibrating against an
  # annual fire-season mean shrinks the predictor variance by roughly an order of
  # magnitude, so B3 gets estimated over a narrow x-range and then extrapolated
  # across the full daily distribution at runtime.  Because site severity also
  # feeds back into combustion buoyancy (FireEvent.cs:606-610), an over-steep B3
  # is exactly how the landscape ends up pinned at the dNBR ceiling.
  #
  # So: use per-fire ignition-date wind when pixel_df carries a usable fire_date,
  # and fall back to the annual fire-season mean only when it does not.
  message("extract_site_predictors: computing EWS...")

  have_dates <- "fire_date" %in% names(pixel_df) &&
                sum(!is.na(pixel_df$fire_date)) > 0
  # Annual CONUS mosaics all report 1 January - that is a mosaic window bound,
  # not a burn date, so treat it as unusable.
  if (have_dates) {
    d <- pixel_df$fire_date[!is.na(pixel_df$fire_date)]
    if (all(format(d, "%m-%d") == "01-01")) {
      message("  fire_date is 1 Jan for every file (annual mosaic bound, not a ",
              "burn date) - falling back to annual fire-season mean wind")
      have_dates <- FALSE
    }
  }

  use_daily <- have_dates && !is.null(wind_daily) && nrow(wind_daily) > 0

  if (use_daily) {
    # NB: temporary column names must NOT start with "." - dplyr 1.1.x fails to
    # resolve bare column references inside a mutate that creates a dot-prefixed
    # name, and as.Date(date) then silently picks up base::date and errors.
    date_col <- intersect(c("date", "Date", "DATE"), names(wind_daily))[1]
    if (is.na(date_col))
      stop("wind_daily has no date column (expected one of: date, Date, DATE). ",
           "Columns present: ", paste(names(wind_daily), collapse = ", "))

    wd_tbl <- wind_daily
    wd_tbl$burn_day <- as.Date(wd_tbl[[date_col]])
    wd_tbl <- wd_tbl %>%
      group_by(burn_day) %>%
      summarise(mean_ws_kmh = mean(WindSpeed_kmh, na.rm = TRUE),
                mean_wd_deg = mean(WindDir_deg,   na.rm = TRUE),
                .groups = "drop")

    pixel_df$burn_day <- as.Date(pixel_df$fire_date)
    pixel_df <- pixel_df %>% left_join(wd_tbl, by = "burn_day")

    n_hit <- sum(is.finite(pixel_df$mean_ws_kmh))
    message(sprintf("  day-of-burn wind matched for %.1f%% of pixels (%d fires)",
                    100 * mean(is.finite(pixel_df$mean_ws_kmh)),
                    length(unique(pixel_df$fire_id))))

    # Backfill any unmatched burn dates from the annual means
    if (n_hit < nrow(pixel_df) && !is.null(wind_annual_tbl)) {
      pixel_df <- pixel_df %>%
        left_join(wind_annual_tbl %>%
                    select(year, ws_yr = mean_ws_kmh, wd_yr = mean_wd_deg),
                  by = c("fire_year" = "year")) %>%
        mutate(mean_ws_kmh = coalesce(mean_ws_kmh, ws_yr),
               mean_wd_deg = coalesce(mean_wd_deg, wd_yr)) %>%
        select(-ws_yr, -wd_yr)
    }
    pixel_df <- pixel_df %>% select(-burn_day)
  }

  if (!use_daily && !is.null(wind_annual_tbl)) {
    warning("Site mortality EWS is being built from ANNUAL fire-season mean wind ",
            "because per-fire burn dates are unavailable. SCF uses day-of-burn ",
            "wind, so SiteMortalityB3 will be fitted over a much narrower range ",
            "than it is applied to. Load continuous per-fire MTBS rasters to get ",
            "ignition dates.")
    pixel_df <- pixel_df %>%
      left_join(wind_annual_tbl %>% select(year, mean_ws_kmh, mean_wd_deg),
                by = c("fire_year" = "year"))
  }

  if (use_daily || !is.null(wind_annual_tbl)) {

    # Try slope/aspect extraction; fall back to flat terrain if small extent
    slope_vals  <- rep(0, nrow(pixel_df))
    aspect_vals <- rep(0, nrow(pixel_df))
    terrain_ok  <- FALSE

    if (!is.null(slope_r) && !is.null(aspect_r)) {
      slope_ext  <- ext(slope_r)
      pv_ext     <- ext(pv)
      # Check if any pixel falls within the slope raster extent
      if (xmin(pv_ext) <= xmax(slope_ext) && xmax(pv_ext) >= xmin(slope_ext) &&
          ymin(pv_ext) <= ymax(slope_ext) && ymax(pv_ext) >= ymin(slope_ext)) {
        sv <- terra::extract(slope_r,  pv)[, 2]
        av <- terra::extract(aspect_r, pv)[, 2]
        in_extent <- !is.na(sv)
        if (any(in_extent)) {
          slope_vals[in_extent]  <- sv[in_extent]
          aspect_vals[in_extent] <- av[in_extent]
          terrain_ok <- TRUE
          pct <- 100 * mean(in_extent)
          message(sprintf("  Terrain coverage: %.1f%% of pixels within slope_r extent", pct))
          if (pct < 50)
            warning(
              "slope_r extent covers <50% of MTBS pixels. ",
              "Remaining pixels use flat terrain (slope=0) for EWS. ",
              "Consider providing a slope raster over the full calibration region."
            )
        }
      }
    }

    if (!terrain_ok)
      message("  No terrain raster or no overlap — using raw wind speed as EWS (slope=0)")

    # compute_effective_wind is defined in spread_helpers.R
    pixel_df$EWS_kmh <- compute_effective_wind(
      wind_speed          = pixel_df$mean_ws_kmh,
      wind_dir            = pixel_df$mean_wd_deg,
      slope_deg           = slope_vals,
      aspect_deg          = aspect_vals,
      combustion_buoyancy = 10
    )

    pixel_df <- pixel_df %>% select(-mean_ws_kmh, -mean_wd_deg)
    attr(pixel_df, "ews_wind_basis") <- if (use_daily) "day-of-burn" else "annual mean"

  } else {
    message("extract_site_predictors: no wind data, EWS set to NA")
    pixel_df$EWS_kmh <- NA_real_
    attr(pixel_df, "ews_wind_basis") <- "none"
  }

  # ---- Fine fuels ------------------------------------------------------------
  # NOTE: constant placeholders are NOT used any more.  A predictor with zero
  # variance cannot be estimated, yet glm() will still return a coefficient for
  # it that silently absorbs part of the intercept - and the intercept is what
  # sets SCF's baseline severity.  Absent predictors are left NA so that
  # fit_site_mortality() drops them and writes an honest 0.0 into the B slot.
  if (!is.null(fine_fuels_r)) {
    message("extract_site_predictors: extracting fine fuels...")
    pixel_df$FineFuels <- terra::extract(fine_fuels_r, pv)[, 2]
  } else {
    message("extract_site_predictors: fine fuels not available - SiteMortalityB5 will be 0")
    pixel_df$FineFuels <- NA_real_
  }

  # ---- Ladder fuels (g m-2 biomass, matching SCF) ----------------------------
  if (!is.null(ladder_fuels_r)) {
    message("extract_site_predictors: extracting ladder fuel biomass...")
    pixel_df$LadderFuels <- terra::extract(ladder_fuels_r, pv)[, 2]
    lv <- pixel_df$LadderFuels[is.finite(pixel_df$LadderFuels)]
    if (length(lv) > 0) {
      message(sprintf("  ladder fuel biomass at pixels: range [%.0f, %.0f] g m-2, mean %.0f",
                      min(lv), max(lv), mean(lv)))
      if (max(lv) <= 1.5)
        warning("Ladder fuel values look like a 0-1 INDEX, not biomass in g m-2. ",
                "SCF Eq. 7 multiplies B6 by raw summed cohort biomass, so an index ",
                "produces a B6 about three orders of magnitude too large. Rebuild ",
                "with fetch_landfire_ladder_fuels(), or set the source to 'none'.")
    }
  } else {
    message("extract_site_predictors: ladder fuels not available - SiteMortalityB6 will be 0")
    pixel_df$LadderFuels <- NA_real_
  }

  # ---- Drop rows with any NA predictor ---------------------------------------
  # Any predictor may be entirely absent; require only those that actually carry
  # data, so a missing fuels layer no longer wipes out the whole sample.
  candidates <- c("Clay", "ET_mm", "CWD_mm", "EWS_kmh", "FineFuels", "LadderFuels")
  required_preds <- candidates[vapply(candidates, function(v)
    v %in% names(pixel_df) && any(is.finite(pixel_df[[v]])), logical(1))]
  n_before <- nrow(pixel_df)
  wind_basis <- attr(pixel_df, "ews_wind_basis")
  pixel_df <- pixel_df %>%
    filter(if_all(all_of(required_preds), is.finite))
  message(sprintf("extract_site_predictors: %d of %d rows retained after NA drop",
                  nrow(pixel_df), n_before))

  # dplyr verbs drop attributes, so re-attach after the last one.
  attr(pixel_df, "ews_wind_basis") <- wind_basis
  # Calibration EWS uses Cb = 10 (SCF's low-intensity default). SCF escalates Cb
  # to 25 above Intensity 3 and 50 above Intensity 6 (FireEvent.cs:606-610), so
  # B3 is applied to systematically larger EWS in severe fire than it was fit on.
  attr(pixel_df, "combustion_buoyancy") <- 10

  pixel_df
}


# =============================================================================
# Site mortality model fitting
# =============================================================================

# -----------------------------------------------------------------------------
# fit_site_mortality
#
# Purpose: Fit a Gamma GLM (inverse link) for dNBR ~ site predictors.
#   Predictors are scaled for convergence; coefficients are unscaled back to
#   original units before returning.
#
# Args:
#   pixel_df  tibble with columns dNBR, Clay, ET_mm, EWS_kmh, CWD_mm,
#             FineFuels, LadderFuels (ET_mm/CWD_mm may be NA)
#
# Returns: list(model, coef_df, data, null_deviance, residual_deviance, n)
# -----------------------------------------------------------------------------

fit_site_mortality <- function(pixel_df) {

  # ---- Determine which predictors are available (all-NA = skip) --------------
  candidate_preds <- c("Clay", "ET_mm", "EWS_kmh", "CWD_mm",
                       "FineFuels", "LadderFuels")

  # A predictor needs data AND variance.  A constant column is not estimable:
  # glm() would hand back a coefficient that quietly absorbs part of the
  # intercept, and the intercept is what sets SCF's baseline severity.
  available_preds <- candidate_preds[vapply(candidate_preds, function(v) {
    if (!v %in% names(pixel_df)) return(FALSE)
    x <- pixel_df[[v]][is.finite(pixel_df[[v]])]
    if (length(x) <= 10) return(FALSE)
    if (sd(x) <= .Machine$double.eps * max(1, abs(mean(x)))) {
      message(sprintf("fit_site_mortality: dropping '%s' - constant (%.4g)", v, x[1]))
      return(FALSE)
    }
    TRUE
  }, logical(1))]

  dropped <- setdiff(candidate_preds, available_preds)
  if (length(dropped) > 0)
    message("fit_site_mortality: not fitted (B slot written as 0.0): ",
            paste(dropped, collapse = ", "))

  if (length(available_preds) == 0)
    stop("No predictor columns with sufficient non-NA data found in pixel_df.")

  # ---- Build analysis dataset ------------------------------------------------
  keep_cols <- c("dNBR", available_preds)
  dat <- pixel_df %>%
    select(all_of(keep_cols)) %>%
    filter(dNBR > 0, if_all(everything(), is.finite))

  if (nrow(dat) < 20)
    stop(sprintf(
      "Only %d rows available for site mortality GLM after filtering. ",
      "Minimum 20 required. Check MTBS pixel loading and predictor coverage.",
      nrow(dat)
    ))

  message(sprintf("fit_site_mortality: fitting Gamma GLM on %d pixels, predictors: %s",
                  nrow(dat), paste(available_preds, collapse = ", ")))

  # ---- Scale predictors for convergence (record params for unscaling) --------
  scale_params <- lapply(available_preds, function(v) {
    list(mean = mean(dat[[v]], na.rm = TRUE), sd = sd(dat[[v]], na.rm = TRUE))
  })
  names(scale_params) <- available_preds

  dat_scaled <- dat
  for (v in available_preds) {
    s <- scale_params[[v]]
    dat_scaled[[v]] <- if (s$sd > 0) (dat[[v]] - s$mean) / s$sd else dat[[v]] - s$mean
  }

  # ---- Build formula ---------------------------------------------------------
  fmla <- as.formula(paste("dNBR ~", paste(available_preds, collapse = " + ")))

  # ---- Fit Gamma(inverse link) -----------------------------------------------
  fit <- tryCatch({
    glm(fmla, data = dat_scaled,
        family  = Gamma(link = "inverse"),
        mustart = rep(mean(dat$dNBR), nrow(dat_scaled)))
  }, warning = function(w) {
    message("Gamma(inverse) warning: ", conditionMessage(w))
    suppressWarnings(
      glm(fmla, data = dat_scaled,
          family  = Gamma(link = "inverse"),
          mustart = rep(mean(dat$dNBR), nrow(dat_scaled)))
    )
  }, error = function(e) {
    message("Gamma(inverse) failed (", conditionMessage(e),
            "). Falling back to Gamma(log).")
    warning("Site mortality: Gamma inverse-link did not converge — using log link. ",
            "Consider verifying predictors and dNBR range.")
    glm(fmla, data = dat_scaled, family = Gamma(link = "log"))
  })

  # ---- Unscale coefficients --------------------------------------------------
  # For linear predictor eta = b0 + b1*x_scaled where x_scaled = (x-mu)/sd:
  # eta = b0 + b1*(x-mu)/sd = (b0 - b1*mu/sd) + (b1/sd)*x
  # -> intercept_orig = b0 - sum(b_j * mu_j / sd_j)
  # -> slope_orig_j   = b_j / sd_j
  coef_scaled <- coef(fit)
  coef_orig   <- coef_scaled  # copy to preserve names

  for (v in available_preds) {
    s         <- scale_params[[v]]
    b_scaled  <- coef_scaled[v]
    if (is.na(b_scaled)) next
    b_orig    <- if (s$sd > 0) b_scaled / s$sd else b_scaled
    coef_orig[v]              <- b_orig
    coef_orig["(Intercept)"] <- coef_orig["(Intercept)"] - b_scaled * s$mean /
                                  max(s$sd, .Machine$double.eps)
  }

  # ---- Standard errors (delta method approximation: SE_orig = SE_scaled/sd) --
  se_scaled <- sqrt(diag(vcov(fit)))
  se_orig   <- se_scaled
  for (v in available_preds) {
    s <- scale_params[[v]]
    if (!is.na(se_orig[v]) && s$sd > 0)
      se_orig[v] <- se_scaled[v] / s$sd
  }

  # ---- p-values from z-scores ------------------------------------------------
  z_vals  <- coef_orig / se_orig
  p_vals  <- 2 * pnorm(abs(z_vals), lower.tail = FALSE)

  # ---- SCF parameter names ---------------------------------------------------
  # SCF uses FIXED Bx indices regardless of which predictors are available.
  # FireEvent.cs reads SiteMortalityB1 as Clay, B2 as ET, B3 as EWS, etc.
  # Missing predictors MUST appear as 0.00 in the parameter file — they cannot
  # be renumbered.  Build a complete 7-row coef_df (B0-B6) with fitted values
  # for available predictors and estimate=0, se=NA, p=NA for absent ones.
  fixed_scf_map <- c(
    "(Intercept)" = "SiteMortalityB0",
    "Clay"        = "SiteMortalityB1",
    "ET_mm"       = "SiteMortalityB2",
    "EWS_kmh"     = "SiteMortalityB3",
    "CWD_mm"      = "SiteMortalityB4",
    "FineFuels"   = "SiteMortalityB5",
    "LadderFuels" = "SiteMortalityB6"
  )

  all_terms <- names(fixed_scf_map)  # canonical order: Intercept, then B1-B6

  coef_df <- tibble(
    term          = all_terms,
    estimate      = vapply(all_terms, function(t)
                      if (t %in% names(coef_orig)) coef_orig[[t]] else 0.0, numeric(1)),
    std_error     = vapply(all_terms, function(t)
                      if (t %in% names(se_orig))   se_orig[[t]]  else NA_real_, numeric(1)),
    p_value       = vapply(all_terms, function(t)
                      if (t %in% names(p_vals))    p_vals[[t]]   else NA_real_, numeric(1)),
    scf_parameter = unname(fixed_scf_map[all_terms]),
    fitted        = all_terms %in% names(coef_orig)   # flag fitted vs zero-padded
  )

  # ---------------------------------------------------------------------------
  # Diagnostic A — linear-predictor positivity / clamp exposure
  #
  # SCF computes siteMortality = max(eta, 0.0005)^-1 (FireEvent.cs:371).  The
  # floor is a cliff, not a soft bound: any cell whose linear predictor reaches
  # zero is assigned dNBR 2000, the maximum.  If that is firing anywhere in the
  # calibration sample it will fire far more often at runtime, where predictors
  # range wider than the fitted data.
  # ---------------------------------------------------------------------------
  eta_hat  <- as.numeric(predict(fit, type = "link"))
  n_clamp  <- sum(eta_hat <= 0.0005, na.rm = TRUE)
  pct_clamp <- 100 * n_clamp / length(eta_hat)

  if (pct_clamp > 0)
    warning(sprintf(paste0(
      "%.2f%% of calibration pixels have a linear predictor at or below the SCF ",
      "floor of 0.0005, so SCF would assign them the maximum dNBR of 2000. This ",
      "almost always means a predictor is on the wrong scale - check that ladder ",
      "fuels are in g m-2 and that fine fuels match MaximumFineFuels."), pct_clamp))

  # ---------------------------------------------------------------------------
  # Diagnostic B — predicted vs observed dNBR RANGE
  #
  # A GLM's fitted values compress variance.  If the site model spans only a
  # narrow band while observed severity spans a wide one, then at runtime the
  # cohort model sees a near-constant SiteMortality and CohortMortalityB2 stops
  # doing anything - mortality becomes species-driven only.
  # ---------------------------------------------------------------------------
  pred_dnbr <- as.numeric(predict(fit, type = "response"))
  obs_q  <- quantile(dat$dNBR,  c(0.05, 0.5, 0.95), na.rm = TRUE)
  pred_q <- quantile(pred_dnbr, c(0.05, 0.5, 0.95), na.rm = TRUE)
  span_ratio <- as.numeric(diff(pred_q[c(1, 3)]) / max(diff(obs_q[c(1, 3)]), 1e-9))

  message(sprintf(
    "fit_site_mortality: observed dNBR 5/50/95 = %.0f/%.0f/%.0f; predicted = %.0f/%.0f/%.0f (span ratio %.2f)",
    obs_q[1], obs_q[2], obs_q[3], pred_q[1], pred_q[2], pred_q[3], span_ratio))

  if (span_ratio < 0.25)
    warning(sprintf(paste0(
      "The site model reproduces only %.0f%% of the observed dNBR spread ",
      "(predicted 5-95%%: %.0f-%.0f vs observed %.0f-%.0f). SCF will feed a ",
      "nearly constant SiteMortality into cohort mortality, so CohortMortalityB2 ",
      "will have little effect and severity will be driven almost entirely by ",
      "species bark thickness."),
      100 * span_ratio, pred_q[1], pred_q[3], obs_q[1], obs_q[3]))

  # ---------------------------------------------------------------------------
  # Diagnostic D — sign guard
  #
  # Inverse link: eta = 1/mu, so anything that INCREASES severity must carry a
  # NEGATIVE coefficient, and B0 must be positive (B0 = 3e-3 -> baseline ~333).
  # ---------------------------------------------------------------------------
  b0 <- coef_df$estimate[coef_df$term == "(Intercept)"]
  sign_warnings <- character(0)
  if (length(b0) == 1 && b0 <= 0)
    sign_warnings <- c(sign_warnings, sprintf(
      "SiteMortalityB0 = %.3e is not positive; SCF needs a positive intercept (1/B0 is baseline dNBR).", b0))

  expect_negative <- c(EWS_kmh = "SiteMortalityB3", FineFuels = "SiteMortalityB5",
                       LadderFuels = "SiteMortalityB6")
  for (v in names(expect_negative)) {
    if (!v %in% available_preds) next
    est <- coef_df$estimate[coef_df$term == v]
    if (length(est) == 1 && est > 0)
      sign_warnings <- c(sign_warnings, sprintf(
        "%s (%s) = %+.3e is positive, implying MORE %s gives LESS severity. Under the inverse link this term should normally be negative.",
        expect_negative[[v]], v, est, v))
  }
  if (length(sign_warnings) > 0)
    warning(paste(c("Site mortality sign check:", sign_warnings), collapse = "\n  - "))

  list(
    model             = fit,
    coef_df           = coef_df,
    data              = dat,
    null_deviance     = fit$null.deviance,
    residual_deviance = fit$deviance,
    n                 = nrow(dat),
    link              = fit$family$link,
    available_preds   = available_preds,
    dropped_preds     = dropped,
    pct_clamped       = pct_clamp,
    obs_dnbr_q        = obs_q,
    pred_dnbr_q       = pred_q,
    span_ratio        = span_ratio,
    sign_warnings     = sign_warnings,
    ews_wind_basis    = if (is.null(attr(pixel_df, "ews_wind_basis")))
                          NA_character_ else attr(pixel_df, "ews_wind_basis")
  )
}


# =============================================================================
# Diagnostic plots
# =============================================================================

# -----------------------------------------------------------------------------
# plot_site_mortality_obs_pred
#
# Purpose: Scatter plot of observed vs predicted dNBR with severity classes.
# Args:  site_fit — return value of fit_site_mortality()
# Returns: ggplot object
# -----------------------------------------------------------------------------

plot_site_mortality_obs_pred <- function(site_fit) {

  dat <- site_fit$data
  dat$predicted <- predict(site_fit$model,
                            newdata = dat,
                            type    = "response")

  # Severity classification (MTBS standard)
  sev_breaks <- c(-Inf, 100, 270, 440, 660, Inf)
  sev_labels <- c("Unburned/Low", "Low", "Moderate", "High", "Very High")
  dat$severity <- cut(dat$dNBR, breaks = sev_breaks, labels = sev_labels)

  # Metrics
  r2   <- cor(dat$dNBR, dat$predicted, use = "complete.obs")^2
  rmse <- sqrt(mean((dat$dNBR - dat$predicted)^2, na.rm = TRUE))

  ggplot(dat, aes(x = dNBR, y = predicted, colour = severity)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
    scale_colour_manual(
      values = c("Unburned/Low" = "grey70", "Low" = "#fee08b",
                 "Moderate"     = "#fc8d59", "High" = "#d73027",
                 "Very High"    = "#762a83"),
      name = "Severity Class"
    ) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = sprintf("R\u00B2 = %.3f\nRMSE = %.1f", r2, rmse),
             size = 3.5) +
    labs(
      title = "Site Mortality: Observed vs Predicted dNBR",
      x     = "Observed dNBR",
      y     = "Predicted dNBR (Gamma GLM)"
    ) +
    theme_bw(base_size = 11)
}


# -----------------------------------------------------------------------------
# plot_site_mortality_coef
#
# Purpose: Coefficient plot (point + CI) for site mortality predictors.
# Args:  site_fit — return value of fit_site_mortality()
# Returns: ggplot object
# -----------------------------------------------------------------------------

plot_site_mortality_coef <- function(site_fit) {

  coef_df <- site_fit$coef_df %>%
    filter(term != "(Intercept)", isTRUE(fitted)) %>%
    mutate(
      ci_lo = estimate - 1.96 * std_error,
      ci_hi = estimate + 1.96 * std_error,
      sign  = ifelse(estimate >= 0, "Increases severity", "Decreases severity")
    )

  ggplot(coef_df, aes(x = estimate, y = reorder(scf_parameter, estimate),
                      colour = sign)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25) +
    geom_point(size = 3) +
    scale_colour_manual(
      values = c("Increases severity" = "#d73027", "Decreases severity" = "#4575b4"),
      name   = NULL
    ) +
    labs(
      title = "Site Mortality Coefficients (95% CI)",
      subtitle = sprintf("Gamma GLM (%s link) — unscaled, original predictor units",
                         site_fit$link),
      x = "Coefficient estimate",
      y = "SCF parameter"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
}


# =============================================================================
# SCF parameter snippet
# =============================================================================

# -----------------------------------------------------------------------------
# write_mortality_snippet
#
# Purpose: Format mortality coefficients for direct copy into an SCF parameter
#   file. Mirrors the style of write_scf_snippet() in spread_helpers.R.
#
# Args:
#   site_fit    return value of fit_site_mortality()
#   cohort_fit  return value of fit_cohort_mortality() or NULL
#
# Returns: character string (multi-line)
# -----------------------------------------------------------------------------

write_mortality_snippet <- function(site_fit, cohort_fit = NULL,
                                    et_r = NULL, cwd_r = NULL,
                                    max_fine_fuels   = NULL,
                                    ladder_max_age   = NULL,
                                    ladder_species   = NULL,
                                    ladder_ceiling_m = NULL) {

  now <- format(Sys.time(), "%Y-%m-%d %H:%M")

  fmt_par <- function(name, value, note) {
    sprintf("%-30s  %+.8f  %s", paste0(name, " ="), value, note)
  }

  lines <- c(
    "# =============================================================================",
    paste0("# SCF Mortality Parameters - generated ", now),
    "# =============================================================================",
    "#",
    "# SITE MORTALITY (Equation 7, FireEvent.cs:363-379)",
    "#   eta  = B0 + B1*Clay + B2*PET + B3*EWS + B4*CWD + B5*FineFuels + B6*LadderFuels",
    "#   dNBR = max(eta, 0.0005) ^ -1        <- Gamma inverse link",
    "#",
    "#   The max(eta, 0.0005) floor is a CLIFF, not a soft bound: any cell whose",
    "#   linear predictor reaches zero is assigned the maximum dNBR of 2000.",
    "#   Because of the inverse link, SIGNS ARE INVERTED - a variable that",
    "#   increases severity carries a NEGATIVE coefficient, and B0 must be",
    "#   positive (1/B0 is the baseline dNBR).",
    "#",
    "# COHORT MORTALITY (Equations 8-9, FireEvent.cs:448-450)",
    "#   P = exp(B0 + B1*BarkThickness + B2*SiteMortality);  p = P/(1+P)",
    "#   One Bernoulli draw per cohort; a kill removes the ENTIRE cohort biomass.",
    "#   BarkThickness is SINGLE-SIDED cm (Cansler/FOFEM BT_coef x DBH).",
    "# =============================================================================",
    ""
  )

  # ---- Site mortality --------------------------------------------------------
  lines <- c(lines,
    "# ---------------------------------------------------------------------------",
    "# Site Mortality",
    paste0("#   Model: Gamma GLM, link = ", site_fit$link),
    paste0("#   N pixels: ", scales::comma(site_fit$n)),
    paste0("#   Predictors fitted: ", paste(site_fit$available_preds, collapse = ", ")),
    sprintf("#   Residual deviance: %.4f  (null: %.4f)",
            site_fit$residual_deviance, site_fit$null_deviance),
    "# ---------------------------------------------------------------------------",
    ""
  )

  if (!is.null(site_fit$dropped_preds) && length(site_fit$dropped_preds) > 0)
    lines <- c(lines,
      paste0("# NOT FITTED (no data or no variance), written as 0.0: ",
             paste(site_fit$dropped_preds, collapse = ", ")),
      "")

  for (i in seq_len(nrow(site_fit$coef_df))) {
    row <- site_fit$coef_df[i, ]
    unit_note <- switch(row$term,
      "(Intercept)" = sprintf("# intercept -> baseline dNBR = %.0f",
                              if (row$estimate > 0) 1 / row$estimate else NA_real_),
      "Clay"        = "# clay FRACTION 0-1 (SCF User Guide 2.18: ClayMap values 0-1)",
      "ET_mm"       = "# previous-year POTENTIAL ET, mm/yr (SCF Succession.PET)",
      "EWS_kmh"     = "# effective wind speed km/h (Nelson 2002, Cb=10)",
      "CWD_mm"      = "# previous-year climatic water deficit, mm/yr (SCF Succession.CWD)",
      "FineFuels"   = "# fine fuel fraction 0-1 -- SEE MaximumFineFuels NOTE BELOW",
      "LadderFuels" = "# ladder fuel BIOMASS, g m-2 (SCF sums raw cohort biomass)",
      ""
    )
    if (!isTRUE(row$fitted) && row$term != "(Intercept)")
      unit_note <- paste(unit_note, "  [not fitted - set to 0]")
    lines <- c(lines, fmt_par(row$scf_parameter, row$estimate, unit_note))
  }

  # ---- Companion parameters that change what the B's mean --------------------
  lines <- c(lines, "",
    "# ---------------------------------------------------------------------------",
    "# COMPANION PARAMETERS - these change the meaning of the coefficients above",
    "# ---------------------------------------------------------------------------")

  lines <- c(lines, "",
    "# MaximumFineFuels rescales the fine fuel term in BOTH Equation 5 (spread)",
    "# and Equation 7 (severity) -- User Guide 2.28. SCF computes",
    "#     fineFuelPercent = min(SiteVars.FineFuels / MaximumFineFuels, 1.0)",
    "# Calibration normalised FBFM40 by 9.0 t/ac (GR9 max), so MaximumFineFuels",
    "# must be set to a comparable load or SiteMortalityB5 (and SpreadProbabilityB2)",
    "# are on the wrong scale.")
  lines <- c(lines,
    if (!is.null(max_fine_fuels) && is.finite(max_fine_fuels))
      sprintf("%-30s  %.2f", "MaximumFineFuels =", max_fine_fuels)
    else
      "# MaximumFineFuels           <SET THIS>   # ~900-1000 g m-2 to match the 9 t/ac normalisation")

  # TimeZero fallbacks
  lines <- c(lines, "",
    "# SCF falls back to these at year 0, before NECN has produced PET/CWD",
    "# (FireEvent.cs:344-360). Leaving them at 0.0 drops two terms from Eq. 7 in",
    "# the first simulated fire year, which -- with negative B2/B4 -- lowers eta",
    "# and therefore RAISES first-year severity.")

  mean_of <- function(r) {
    if (is.null(r)) return(NA_real_)
    v <- try(terra::values(r, mat = FALSE), silent = TRUE)
    if (inherits(v, "try-error")) return(NA_real_)
    v <- v[is.finite(v)]
    if (length(v) == 0) NA_real_ else mean(v)
  }
  tz_pet <- mean_of(et_r)
  tz_cwd <- mean_of(cwd_r)
  if (!is.finite(tz_pet) && "ET_mm"  %in% names(site_fit$data)) tz_pet <- mean(site_fit$data$ET_mm,  na.rm = TRUE)
  if (!is.finite(tz_cwd) && "CWD_mm" %in% names(site_fit$data)) tz_cwd <- mean(site_fit$data$CWD_mm, na.rm = TRUE)

  lines <- c(lines,
    if (is.finite(tz_pet)) sprintf("%-30s  %.2f   # landscape mean PET, mm/yr", "TimeZeroPET =", tz_pet)
    else "# TimeZeroPET                <SET THIS>   # landscape mean PET, mm/yr",
    if (is.finite(tz_cwd)) sprintf("%-30s  %.2f   # landscape mean CWD, mm/yr", "TimeZeroCWD =", tz_cwd)
    else "# TimeZeroCWD                <SET THIS>   # landscape mean CWD, mm/yr")

  # Ladder fuel definition
  b6 <- site_fit$coef_df$estimate[site_fit$coef_df$scf_parameter == "SiteMortalityB6"]
  lines <- c(lines, "",
    "# LadderFuelMaxAge / LadderFuelSpeciesList define the biomass SCF sums into",
    "# the B6 term. They MUST match the stratum B6 was calibrated on.")
  if (!is.null(ladder_ceiling_m) && is.finite(ladder_ceiling_m))
    lines <- c(lines, sprintf(
      "# Calibrated against LANDFIRE canopy biomass below %.1f m, in g m-2.",
      ladder_ceiling_m))
  if (length(b6) == 1 && b6 == 0) {
    lines <- c(lines,
      "# SiteMortalityB6 is 0.0, so the ladder fuel term is inactive. The two",
      "# parameters below are still required by the parser.")
  }
  lines <- c(lines,
    if (!is.null(ladder_max_age) && is.finite(ladder_max_age))
      sprintf("%-30s  %d", "LadderFuelMaxAge =", as.integer(ladder_max_age))
    else "# LadderFuelMaxAge           <SET THIS>   # age below which a cohort is ladder fuel",
    if (!is.null(ladder_species) && length(ladder_species) > 0)
      paste0("LadderFuelSpeciesList  ", paste(ladder_species, collapse = "  "))
    else "# LadderFuelSpeciesList      <SET THIS>   # species codes acting as ladder fuel")

  # ---- Diagnostics -----------------------------------------------------------
  lines <- c(lines, "",
    "# ---------------------------------------------------------------------------",
    "# Diagnostics",
    "# ---------------------------------------------------------------------------")

  if (!is.null(site_fit$pct_clamped))
    lines <- c(lines, sprintf(
      "#   Clamp exposure: %.2f%% of calibration pixels have eta <= 0.0005 (SCF would give them dNBR 2000)",
      site_fit$pct_clamped))
  if (!is.null(site_fit$obs_dnbr_q) && !is.null(site_fit$pred_dnbr_q))
    lines <- c(lines,
      sprintf("#   Observed dNBR  5/50/95: %.0f / %.0f / %.0f",
              site_fit$obs_dnbr_q[1], site_fit$obs_dnbr_q[2], site_fit$obs_dnbr_q[3]),
      sprintf("#   Predicted dNBR 5/50/95: %.0f / %.0f / %.0f  (span ratio %.2f)",
              site_fit$pred_dnbr_q[1], site_fit$pred_dnbr_q[2], site_fit$pred_dnbr_q[3],
              site_fit$span_ratio),
      "#   A low span ratio means SCF feeds a near-constant SiteMortality into",
      "#   cohort mortality, leaving CohortMortalityB2 with little to act on.")
  if (!is.null(site_fit$ews_wind_basis) && !is.na(site_fit$ews_wind_basis))
    lines <- c(lines, sprintf("#   EWS wind basis: %s (SCF uses day-of-burn wind)",
                              site_fit$ews_wind_basis))
  if (!is.null(site_fit$sign_warnings) && length(site_fit$sign_warnings) > 0)
    lines <- c(lines, "#   SIGN WARNINGS:",
               paste0("#     - ", site_fit$sign_warnings))

  # ---- Cohort mortality ------------------------------------------------------
  if (!is.null(cohort_fit)) {
    lines <- c(lines, "",
      "# ---------------------------------------------------------------------------",
      "# Cohort Mortality (Equations 8-9)",
      "#   P(mortality) = P/(1+P),  P = exp(B0 + B1*BarkThickness + B2*SiteMortality)",
      "#   Model: Binomial GLM, logit link",
      paste0("#   N trees: ", scales::comma(cohort_fit$n),
             if (!is.null(cohort_fit$n_fires))
               paste0(" across ", cohort_fit$n_fires, " fires") else ""),
      paste0("#   SEs: ", if (isTRUE(cohort_fit$cluster_robust))
                            "cluster-robust by fire" else
                            "NAIVE - p-values are anticonservative"),
      "# ---------------------------------------------------------------------------",
      ""
    )
    for (i in seq_len(nrow(cohort_fit$coef_df))) {
      row <- cohort_fit$coef_df[i, ]
      unit_note <- switch(row$term,
        "(Intercept)"   = "# intercept",
        "BarkThickness" = "# SINGLE-SIDED bark thickness, cm -- must match Species_CSV_File MaximumBarkThickness",
        "dNBR"          = "# site-level dNBR 0-2000 (SCF passes its own predicted dNBR here)",
        ""
      )
      lines <- c(lines, fmt_par(row$scf_parameter, row$estimate, unit_note))
    }
    if (!is.null(cohort_fit$n_dnbr_levels))
      lines <- c(lines, "",
        sprintf("#   Distinct dNBR values available to B2: %d", cohort_fit$n_dnbr_levels))
    if (!is.null(cohort_fit$sign_warnings) && length(cohort_fit$sign_warnings) > 0)
      lines <- c(lines, "#   SIGN WARNINGS:",
                 paste0("#     - ", cohort_fit$sign_warnings))
  }

  paste(lines, collapse = "\n")
}


# =============================================================================
# FTM cohort mortality data loading
# =============================================================================

# -----------------------------------------------------------------------------
# load_ftm_cohort_data
#
# Purpose: Load the Cansler et al. (2020) Fire Tree Mortality database and
#   return a tree-level tibble ready for fit_cohort_mortality().
#
# Reads three CSVs from ftm_dir:
#   FTM_fires.csv            -- fire locations (Latitude, Longitude, yr_fire)
#   FTM_trees.csv            -- tree records (Genus_species, DBH_cm, yr1status)
#   Species_BarkThickness.csv -- BT_coef per species for BT = BT_coef * DBH_cm
#
# Mortality status coding (Cansler et al. 2020 FTM metadata, RDS-2020-0001-2):
#   yr1status 0 = alive at first post-fire survey -> dead = 0
#   yr1status 1 = dead at first post-fire survey  -> dead = 1
#   NA           = not assessed in that timestep (2nd edition: no longer filled-in)
#   Trees with NA yr1status are dropped.
#
# Return value: list(tree_df, yr1status_counts)
#   tree_df          — standard tibble suitable for fit_cohort_mortality()
#   yr1status_counts — tibble(yr1status, n) of raw counts before any status
#                      filtering, used by diagnostic plot functions
#
# dNBR extraction: for each FTM fire, load the matching annual MTBS GeoTIFF
#   from mtbs_dnbr_dir (matched by yr_fire), buffer the fire location by
#   buffer_m metres, and take the median dNBR of pixels in that buffer.
#   Fires with no matching MTBS file or no pixels in buffer are dropped.
#
# Args:
#   ftm_dir       directory containing FTM_fires.csv, FTM_trees.csv,
#                 Species_BarkThickness.csv
#   cal_vect      SpatVector of calibration boundary (any CRS)
#   working_crs   CRS string (e.g. "EPSG:32610")
#   mtbs_dnbr_dir directory of annual *_dnbr.tif files from
#                 fetch_mtbs_from_imageserver() or load_mtbs_dnbr_pixels()
#   min_dbh_cm    minimum DBH to include (default 2 cm)
#   buffer_m      radius in metres around fire location for dNBR extraction
#                 (default 1000 m)
#   progress_fn   optional function(msg)
#
# Returns: tibble(YrFireName, yr_fire, dead, BarkThickness, dNBR,
#                 Genus_species, DBH_cm, BT_coef, Latitude, Longitude)
# -----------------------------------------------------------------------------

load_ftm_cohort_data <- function(ftm_dir, cal_vect, working_crs,
                                  mtbs_dnbr_dir,
                                  min_dbh_cm  = 2.0,
                                  buffer_m    = 1000,
                                  progress_fn = NULL) {

  # ---- Verify expected files exist -------------------------------------------
  fires_path <- file.path(ftm_dir, "FTM_fires.csv")
  trees_path <- file.path(ftm_dir, "FTM_trees.csv")
  bark_path  <- file.path(ftm_dir, "Species_BarkThickness.csv")

  missing <- c(fires_path, trees_path, bark_path)[
    !file.exists(c(fires_path, trees_path, bark_path))]
  if (length(missing) > 0)
    stop(sprintf(
      "Missing FTM files in '%s': %s",
      ftm_dir, paste(basename(missing), collapse = ", ")
    ))

  if (!is.null(progress_fn)) progress_fn("Reading FTM CSV files...")
  latin1 <- readr::locale(encoding = "latin1")
  fires <- readr::read_csv(fires_path, show_col_types = FALSE, locale = latin1)
  trees <- readr::read_csv(trees_path, show_col_types = FALSE, locale = latin1)
  bark  <- readr::read_csv(bark_path,  show_col_types = FALSE, locale = latin1)

  # Sanitise any remaining non-UTF-8 bytes in character columns
  # (some species names contain Latin-1 diacritics not caught by locale above)
  clean_utf8 <- function(x) iconv(x, from = "latin1", to = "UTF-8", sub = "")

  message(sprintf("load_ftm_cohort_data: %d fires, %d trees read",
                  nrow(fires), nrow(trees)))

  # ---- Filter fires to cal_vect ----------------------------------------------
  if (!is.null(progress_fn)) progress_fn("Filtering FTM fires to calibration area...")

  fires_valid <- fires %>%
    filter(!is.na(Latitude), !is.na(Longitude))

  fires_sv   <- vect(fires_valid, geom = c("Longitude", "Latitude"),
                     crs = "EPSG:4326")
  cal_wgs84  <- project(cal_vect, "EPSG:4326")

  # Use a generous 50 km buffer around cal_vect to capture nearby fires
  cal_buf    <- buffer(cal_wgs84, 50000)
  in_area    <- as.logical(relate(fires_sv, cal_buf, relation = "intersects"))
  fires_cal  <- fires_valid[in_area, ]

  if (nrow(fires_cal) == 0)
    stop(paste0(
      "No FTM fires found within 50 km of cal_vect. ",
      "Check that cal_vect covers the study area and that Latitude/Longitude ",
      "columns are present in FTM_fires.csv."
    ))

  message(sprintf("  %d of %d fires fall within/near cal_vect",
                  nrow(fires_cal), nrow(fires_valid)))

  # ---- Join trees to filtered fires ------------------------------------------
  if (!is.null(progress_fn)) progress_fn("Joining trees to fires...")

  # Deduplicate fires_cal: same fire can appear in multiple datasets
  fires_cal <- fires_cal %>%
    distinct(YrFireName, .keep_all = TRUE)

  trees_cal <- inner_join(trees, fires_cal, by = "YrFireName",
                           relationship = "many-to-many")

  if (nrow(trees_cal) == 0)
    stop("No trees matched the in-area fires. Check YrFireName join key.")

  # ---- Determine mortality status from first post-fire year ------------------
  # FTM metadata (RDS-2020-0001-2) definition for all yrXstatus columns:
  #   0 = alive,  1 = dead,  NA = not assessed in that timestep
  if (!"yr1status" %in% names(trees_cal))
    stop("yr1status column not found in FTM_trees.csv.")

  # Save raw distribution BEFORE filtering for diagnostic plots
  yr1status_counts <- trees_cal %>%
    count(yr1status, name = "n") %>%     # NA is its own bucket
    arrange(yr1status)

  message("  yr1status raw distribution (before filtering):")
  for (k in seq_len(nrow(yr1status_counts)))
    message(sprintf("    yr1status=%s  n=%d",
                    ifelse(is.na(yr1status_counts$yr1status[k]), "NA",
                           as.character(yr1status_counts$yr1status[k])),
                    yr1status_counts$n[k]))

  # Keep trees with a measured status (0=alive or 1=dead); drop NA (not assessed)
  trees_cal <- trees_cal %>%
    filter(!is.na(yr1status), yr1status %in% c(0L, 1L)) %>%
    mutate(dead = as.integer(yr1status == 1L))

  # Warn if no dead trees found (unusual but possible for some regions/fires)
  n_dead_initial <- sum(trees_cal$dead)
  if (n_dead_initial == 0) {
    warning(paste0(
      "No yr1status=1 (dead) trees found within/near cal_vect. ",
      "The cohort mortality GLM requires both alive AND dead trees. ",
      "Cohort mortality fitting will fail. ",
      "This may indicate the fires in your calibration region had very low mortality, ",
      "or that dead trees were not re-sampled within the FTM dataset."
    ))
  }

  # ---- Filter by DBH ---------------------------------------------------------
  trees_cal <- trees_cal %>%
    filter(!is.na(DBH_cm), DBH_cm >= min_dbh_cm)

  message(sprintf("  %d trees with yr1status 0 or 1 (dead=%d, alive=%d), DBH >= %.1f cm",
                  nrow(trees_cal), sum(trees_cal$dead), sum(trees_cal$dead == 0),
                  min_dbh_cm))

  # ---- Join bark thickness coefficient ---------------------------------------
  if (!is.null(progress_fn)) progress_fn("Joining bark thickness coefficients...")

  # Normalise whitespace and case for join
  # Normalise genus-species strings for joining:
  #   FTM_trees uses underscores:        "Abies_amabilis", "Abies_NA"
  #   Species_BarkThickness uses spaces: "Abies amabilis", "Abies " (genus-only)
  # Steps: clean UTF-8 → replace underscores with spaces → strip trailing " NA"
  #        (genus-only IDs) → trimws
  norm_gs <- function(x) {
    x <- trimws(clean_utf8(as.character(x)))
    x <- gsub("_", " ", x)            # underscore → space
    x <- gsub("\\s+NA$", "", x)       # "Abies NA" → "Abies"
    trimws(x)
  }

  bark_clean <- bark %>%
    rename(Genus_species = Genus_Species) %>%
    mutate(Genus_species = norm_gs(Genus_species)) %>%
    select(Genus_species, BT_coef)

  trees_cal <- trees_cal %>%
    mutate(Genus_species = norm_gs(Genus_species)) %>%
    left_join(bark_clean, by = "Genus_species") %>%
    mutate(BarkThickness = BT_coef * DBH_cm) %>%
    filter(!is.na(BarkThickness), BarkThickness > 0)

  n_no_bark <- sum(is.na(trees_cal$BT_coef))
  if (n_no_bark > 0) {
    unmatched <- trees_cal %>%
      filter(is.na(BT_coef)) %>%
      count(Genus_species, sort = TRUE) %>%
      head(10)
    message(sprintf("  %d trees have NO bark thickness match. Top unmatched species:",
                    n_no_bark))
    for (k in seq_len(nrow(unmatched)))
      message(sprintf("    '%s'  (n=%d)", unmatched$Genus_species[k], unmatched$n[k]))
  }
  message(sprintf("  %d trees have bark thickness coefficients; %d dropped (no match)",
                  sum(!is.na(trees_cal$BT_coef)), n_no_bark))

  if (nrow(trees_cal) == 0)
    stop(paste0(
      "No trees remain after bark thickness join. ",
      "Check that Genus_species values in FTM_trees.csv match ",
      "Genus_Species in Species_BarkThickness.csv."
    ))

  # ---- Extract dNBR from MTBS at each fire location --------------------------
  if (!is.null(progress_fn)) progress_fn("Extracting dNBR at fire locations from MTBS...")

  # Find annual MTBS TIF files
  all_tifs <- list.files(mtbs_dnbr_dir, pattern = "_dnbr\\.tif$",
                          full.names = TRUE)
  all_tifs <- all_tifs[!grepl("_dnbr6\\.tif$", all_tifs)]

  if (length(all_tifs) == 0)
    stop(sprintf(
      "No *_dnbr.tif files found in '%s'. Load MTBS data first.", mtbs_dnbr_dir
    ))

  # Parse year from filename (first 8-digit block -> YYYY)
  tif_years <- vapply(basename(all_tifs), function(f) {
    b <- regmatches(f, gregexpr("[0-9]{8}", f))[[1]]
    if (length(b) == 0) NA_integer_ else as.integer(substr(b[1], 1, 4))
  }, integer(1))

  # Unique fires for extraction
  fires_uniq <- fires_cal %>%
    select(YrFireName, yr_fire, Latitude, Longitude) %>%
    distinct()

  fire_dnbr_list <- vector("list", nrow(fires_uniq))

  for (i in seq_len(nrow(fires_uniq))) {
    frow <- fires_uniq[i, ]
    yr   <- as.integer(frow$yr_fire)

    # Match MTBS TIF for this year
    idx <- which(tif_years == yr)
    if (length(idx) == 0) {
      message(sprintf("  %s (yr %d): no MTBS file found", frow$YrFireName, yr))
      fire_dnbr_list[[i]] <- tibble(YrFireName = frow$YrFireName,
                                     dNBR = NA_real_)
      next
    }

    r <- tryCatch(rast(all_tifs[idx[1]]), error = function(e) {
      message(sprintf("  %s: cannot read MTBS TIF — %s",
                      frow$YrFireName, conditionMessage(e)))
      NULL
    })
    if (is.null(r)) {
      fire_dnbr_list[[i]] <- tibble(YrFireName = frow$YrFireName,
                                     dNBR = NA_real_)
      next
    }

    # Reclassify raw severity classes if needed (same logic as load_mtbs_dnbr_pixels)
    rv  <- values(r, mat = FALSE)
    nna <- sum(!is.na(rv))
    if (nna > 0 && max(rv, na.rm = TRUE) <= 6) {
      r <- subst(r, from = c(0L, 1L, 2L, 3L, 4L, 5L, 6L),
                    to   = c(NA,  NA, 185, 355, 600,  NA, NA))
    }

    # Buffer around fire location and extract median dNBR
    pt     <- vect(data.frame(x = frow$Longitude, y = frow$Latitude),
                   geom = c("x", "y"), crs = "EPSG:4326")
    pt_r   <- project(pt, crs(r))
    pt_buf <- buffer(pt_r, buffer_m)

    pix <- tryCatch(
      terra::extract(r, pt_buf, fun = NULL, ID = FALSE),
      error = function(e) NULL
    )

    if (is.null(pix) || nrow(pix) == 0 || ncol(pix) == 0) {
      message(sprintf("  %s: no pixels in %d m buffer", frow$YrFireName, buffer_m))
      fire_dnbr_list[[i]] <- tibble(YrFireName = frow$YrFireName,
                                     dNBR = NA_real_)
      next
    }

    vals_buf <- pix[[1]]
    med_dnbr <- median(vals_buf[is.finite(vals_buf)], na.rm = TRUE)

    message(sprintf("  %s (yr %d): %d pixels in buffer, median dNBR = %.0f",
                    frow$YrFireName, yr,
                    sum(is.finite(vals_buf)), coalesce(med_dnbr, NA_real_)))

    fire_dnbr_list[[i]] <- tibble(YrFireName = frow$YrFireName,
                                   dNBR = med_dnbr)
  }

  fire_dnbr <- bind_rows(fire_dnbr_list)

  # ---- Join dNBR back to tree data and finalise ------------------------------
  result <- trees_cal %>%
    left_join(fire_dnbr, by = "YrFireName") %>%
    filter(is.finite(dNBR)) %>%
    select(YrFireName, yr_fire, dead, BarkThickness, dNBR,
           Genus_species, DBH_cm, BT_coef, Latitude, Longitude)

  n_fires_matched <- length(unique(result$YrFireName))
  message(sprintf(
    "load_ftm_cohort_data: %d trees from %d fires with dNBR (dead=%d, alive=%d)",
    nrow(result), n_fires_matched, sum(result$dead), sum(result$dead == 0)
  ))

  if (nrow(result) == 0)
    stop(paste0(
      "No trees remain after dNBR extraction. Possible causes: ",
      "(1) MTBS files don't cover the FTM fire years, ",
      "(2) fire locations fall outside the MTBS raster extent, or ",
      "(3) all MTBS pixels within ", buffer_m, " m of fire locations are NA."
    ))

  list(
    tree_df          = result,
    yr1status_counts = yr1status_counts
  )
}


# =============================================================================
# Cohort mortality diagnostic plots
# =============================================================================

# -----------------------------------------------------------------------------
# plot_yr1status_distribution
#
# Purpose: Bar chart of raw yr1status values from FTM trees in the calibration
#   area. This is the most diagnostic figure: FTM codes yr1status as
#   0 = alive, 1 = dead, NA = not assessed (RDS-2020-0001-2). Both 0 and 1 must
#   be present or the binomial GLM cannot be fitted.
#
# Args:
#   yr1status_counts — tibble(yr1status, n) from load_ftm_cohort_data()$yr1status_counts
#
# Returns: ggplot object
# -----------------------------------------------------------------------------

plot_yr1status_distribution <- function(yr1status_counts) {

  # FTM coding: 0=alive, 1=dead, NA=not assessed
  status_colors <- c("0" = "#2c7bb6", "1" = "#d7191c", "NA" = "#aaaaaa")

  df <- yr1status_counts %>%
    mutate(
      key   = ifelse(is.na(yr1status), "NA", as.character(yr1status)),
      label = dplyr::case_when(
        key == "0"  ~ "0: Alive",
        key == "1"  ~ "1: Dead",
        key == "NA" ~ "NA: Not assessed\n(not measured that year)",
        TRUE        ~ paste0("Other: ", yr1status)
      )
    )

  ggplot(df, aes(x = label, y = n, fill = key)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = scales::comma(n)), vjust = -0.4, size = 3.5) +
    scale_fill_manual(
      values = status_colors,
      guide  = "none"
    ) +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.12))) +
    labs(
      title    = "FTM yr1status Distribution (Calibration Region)",
      subtitle = paste0("FTM coding: 0=alive, 1=dead, NA=not assessed.\n",
                        "Both alive (0) and dead (1) trees must be present to fit cohort mortality."),
      x        = NULL,
      y        = "Number of trees"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(size = 9))
}


# -----------------------------------------------------------------------------
# plot_bt_by_status
#
# Purpose: Violin + box plot of BarkThickness by dead/alive status.
#   Ecologically: dead trees should have THINNER bark (less fire protection).
#   If dead trees show THICKER bark the dead/alive coding is inverted.
#
# Args:
#   tree_df — tibble from load_ftm_cohort_data()$tree_df (requires dead, BarkThickness)
#
# Returns: ggplot object
# -----------------------------------------------------------------------------

plot_bt_by_status <- function(tree_df) {

  n_alive <- sum(tree_df$dead == 0, na.rm = TRUE)
  n_dead  <- sum(tree_df$dead == 1, na.rm = TRUE)

  # FTM coding: dead=0 means yr1status=0 (alive), dead=1 means yr1status=1 (dead)
  df <- tree_df %>%
    filter(is.finite(dead), is.finite(BarkThickness)) %>%
    mutate(status = factor(dead, levels = c(0L, 1L),
                           labels = c(
                             sprintf("Alive (yr1=0)\nn=%s", scales::comma(n_alive)),
                             sprintf("Dead (yr1=1)\nn=%s",  scales::comma(n_dead))
                           )))

  fill_vals <- setNames(
    c("#2c7bb6", "#d7191c"),
    c(sprintf("Alive (yr1=0)\nn=%s", scales::comma(n_alive)),
      sprintf("Dead (yr1=1)\nn=%s",  scales::comma(n_dead)))
  )

  p <- ggplot(df, aes(x = status, y = BarkThickness, fill = status))

  # Only add violin when both groups have >=3 observations (ggplot requirement)
  if (n_alive >= 3 && n_dead >= 3) {
    p <- p + geom_violin(alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75))
  }

  p +
    geom_boxplot(width = 0.15, fill = "white",
                 outlier.size = 0.4, outlier.alpha = 0.2) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    labs(
      title    = "Bark Thickness by Survival Status",
      subtitle = if (n_dead == 0)
        "No dead trees (yr1status=1) found in calibration region \u2014 check yr1status distribution"
      else
        "Dead trees should have THINNER bark than alive trees",
      x        = NULL,
      y        = "Bark Thickness (cm)"
    ) +
    theme_bw(base_size = 11)
}


# -----------------------------------------------------------------------------
# plot_mortality_by_bt_quintile
#
# Purpose: Observed mortality rate by BarkThickness quintile.
#   Should show DECREASING mortality with increasing bark thickness (Q1 → Q5).
#
# Args:
#   tree_df — tibble from load_ftm_cohort_data()$tree_df
#
# Returns: ggplot object
# -----------------------------------------------------------------------------

plot_mortality_by_bt_quintile <- function(tree_df) {

  df <- tree_df %>%
    filter(is.finite(dead), is.finite(BarkThickness)) %>%
    mutate(bt_q = ntile(BarkThickness, 5)) %>%
    group_by(bt_q) %>%
    summarise(
      n         = n(),
      n_dead    = sum(dead),
      mort_rate = mean(dead),
      bt_lo     = min(BarkThickness),
      bt_hi     = max(BarkThickness),
      .groups   = "drop"
    ) %>%
    mutate(
      q_label = sprintf("Q%d\n(%.1f\u2013%.1f cm)\nn=%s", bt_q, bt_lo, bt_hi,
                        scales::comma(n))
    )

  ggplot(df, aes(x = bt_q, y = mort_rate)) +
    geom_col(fill = "#fc8d59", width = 0.65, colour = "white") +
    geom_text(aes(label = sprintf("%.0f%%", 100 * mort_rate)),
              vjust = -0.5, size = 3.5) +
    scale_x_continuous(breaks = 1:5, labels = df$q_label) +
    scale_y_continuous(labels = scales::percent,
                       limits = c(0, max(df$mort_rate) * 1.15 + 0.05),
                       expand = expansion(mult = c(0, 0))) +
    labs(
      title    = "Observed Mortality Rate by Bark Thickness Quintile",
      subtitle = "Should DECREASE from Q1 (thinnest) to Q5 (thickest bark = more protection)",
      x        = "Bark Thickness Quintile",
      y        = "Mortality Rate"
    ) +
    theme_bw(base_size = 11)
}


# -----------------------------------------------------------------------------
# plot_mortality_by_severity
#
# Purpose: Observed mortality rate by MTBS fire severity class (dNBR).
#   Should show INCREASING mortality with higher severity.
#
# Args:
#   tree_df — tibble from load_ftm_cohort_data()$tree_df
#
# Returns: ggplot object
# -----------------------------------------------------------------------------

plot_mortality_by_severity <- function(tree_df) {

  sev_breaks <- c(-Inf, 270, 440, Inf)
  sev_labels <- c("Low (<270)", "Moderate (270\u2013440)", "High (>440)")
  sev_colors <- c("Low (<270)" = "#fee08b",
                  "Moderate (270\u2013440)" = "#fc8d59",
                  "High (>440)" = "#d73027")

  df <- tree_df %>%
    filter(is.finite(dead), is.finite(dNBR)) %>%
    mutate(severity = cut(dNBR, breaks = sev_breaks, labels = sev_labels)) %>%
    filter(!is.na(severity)) %>%
    group_by(severity) %>%
    summarise(
      n         = n(),
      mort_rate = mean(dead),
      .groups   = "drop"
    ) %>%
    mutate(
      sev_label = sprintf("%s\nn=%s", severity, scales::comma(n))
    )

  ggplot(df, aes(x = severity, y = mort_rate, fill = severity)) +
    geom_col(width = 0.65, colour = "white") +
    geom_text(aes(label = sprintf("%.0f%%", 100 * mort_rate)),
              vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = sev_colors, guide = "none") +
    scale_x_discrete(labels = df$sev_label) +
    scale_y_continuous(labels = scales::percent,
                       limits = c(0, max(df$mort_rate) * 1.15 + 0.05),
                       expand = expansion(mult = c(0, 0))) +
    labs(
      title    = "Observed Mortality Rate by Fire Severity Class",
      subtitle = "Should INCREASE from Low to High severity (higher dNBR \u2192 more mortality)",
      x        = "MTBS Severity Class (dNBR)",
      y        = "Mortality Rate"
    ) +
    theme_bw(base_size = 11)
}


# =============================================================================
# Cohort mortality GLM fitting
# =============================================================================

# -----------------------------------------------------------------------------
# fit_cohort_mortality
#
# Purpose: Fit a binomial GLM to predict individual tree mortality from
#   bark thickness and fire severity (dNBR).
#
# SCF Eq. 8/9:
#   P(mortality) = logistic(B0 + B1*BarkThickness + B2*dNBR)
#   B0 -> CohortMortalityB0
#   B1 -> CohortMortalityB1   (expect NEGATIVE: thicker bark protects against mortality)
#   B2 -> CohortMortalityB2   (expect POSITIVE: higher dNBR = more severe fire = more mortality)
#
# FTM dead column: 0 = alive (yr1status=0), 1 = dead (yr1status=1)
#
# Args:
#   tree_df  tibble from load_ftm_cohort_data() with columns:
#            dead (0/1), BarkThickness (cm), dNBR
#
# Returns: list(model, coef_df, data, n, n_dead, n_alive,
#               null_deviance, residual_deviance)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# cluster_robust_vcov
#
# Cluster-robust (CR0 / Liang-Zeger) covariance matrix for a fitted glm, with the
# usual finite-sample correction.  Implemented directly so the app does not need
# the 'sandwich' package; it reproduces sandwich::vcovCL(type = "HC0",
# cadjust = TRUE) to numerical precision.
#
#   V = B * ( G/(G-1) * (N-1)/(N-K) * SUM_g s_g s_g' ) * B
#
# where B = vcov(fit) is the model "bread", and s_g is the sum of the score
# (estimating-function) contributions within cluster g.  For a GLM the score
# contribution of observation i is x_i * (y_i - mu_i) * w_i / dispersion.
#
# Args:
#   fit      a glm object
#   cluster  vector of cluster ids, length = nobs(fit)
# Returns: covariance matrix with the same dim/names as vcov(fit)
# -----------------------------------------------------------------------------
cluster_robust_vcov <- function(fit, cluster) {

  X <- stats::model.matrix(fit)
  if (length(cluster) != nrow(X))
    stop(sprintf("cluster length (%d) does not match model rows (%d). ",
                 length(cluster), nrow(X)),
         "This usually means rows were dropped by na.action after the cluster ",
         "vector was built.")

  # Score / estimating-function contributions, one row per observation
  r  <- stats::residuals(fit, type = "response")
  w  <- stats::weights(fit, type = "prior")
  if (is.null(w)) w <- rep(1, nrow(X))
  disp <- summary(fit)$dispersion
  if (!is.finite(disp) || disp <= 0) disp <- 1
  ef <- X * as.numeric(r * w) / disp

  # Sum contributions within each cluster
  cl  <- factor(cluster)
  s_g <- rowsum(ef, group = cl, reorder = FALSE)

  G <- nlevels(cl)
  N <- nrow(X)
  K <- ncol(X)
  adj <- (G / (G - 1)) * ((N - 1) / (N - K))

  bread <- stats::vcov(fit)
  meat  <- crossprod(s_g) * adj
  V     <- bread %*% meat %*% bread

  dimnames(V) <- dimnames(bread)
  V
}


fit_cohort_mortality <- function(tree_df) {

  dat <- tree_df %>%
    filter(is.finite(dead), is.finite(BarkThickness), is.finite(dNBR))

  if (nrow(dat) < 20)
    stop(sprintf(
      "Too few trees for cohort mortality GLM (n = %d; need >= 20). ",
      nrow(dat)
    ))

  if (length(unique(dat$dead)) < 2)
    stop(sprintf(paste0(
      "Only one outcome class present (all dead=%d). ",
      "Both alive and dead trees are required to fit the binomial GLM."
    ), unique(dat$dead)))

  fit <- tryCatch(
    glm(dead ~ BarkThickness + dNBR, data = dat,
        family = binomial(link = "logit")),
    error   = function(e) stop("Cohort mortality GLM failed: ", conditionMessage(e)),
    warning = function(w) {
      message("fit_cohort_mortality warning: ", conditionMessage(w))
      glm(dead ~ BarkThickness + dNBR, data = dat,
          family = binomial(link = "logit"))
    }
  )

  # ---------------------------------------------------------------------------
  # CLUSTER-ROBUST STANDARD ERRORS
  #
  # FTM publishes fire-level coordinates only, so every tree in a fire carries
  # the SAME dNBR (see load_ftm_cohort_data - one median per fire).  Thousands of
  # rows therefore share one value of the key predictor.  Treating them as
  # independent makes the naive SE on CohortMortalityB2 far too small and its
  # p-value look decisive no matter how weak the evidence.  Cluster by fire.
  #
  # This affects INFERENCE only - the point estimates SCF consumes are unchanged.
  # ---------------------------------------------------------------------------
  n_clusters <- length(unique(dat$YrFireName))
  cluster_ok <- FALSE

  tidy_fit <- broom::tidy(fit)

  if (n_clusters > 1) {
    vc <- tryCatch(
      cluster_robust_vcov(fit, cluster = dat$YrFireName),
      error = function(e) { message("cluster_robust_vcov failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(vc)) {
      se_cl <- sqrt(diag(vc))
      est   <- coef(fit)
      keep  <- intersect(names(est), names(se_cl))
      z     <- est[keep] / se_cl[keep]
      tidy_fit <- tibble(
        term      = keep,
        estimate  = as.numeric(est[keep]),
        std.error = as.numeric(se_cl[keep]),
        p.value   = 2 * pnorm(abs(as.numeric(z)), lower.tail = FALSE)
      )
      cluster_ok <- TRUE
      message(sprintf(
        "fit_cohort_mortality: cluster-robust SEs over %d fires (naive SEs would be ~%.1fx smaller)",
        n_clusters,
        mean(sqrt(diag(vcov(fit)))[keep] / se_cl[keep], na.rm = TRUE)))
    }
  }

  if (!cluster_ok) {
    warning(sprintf(paste0(
      "Cohort mortality SEs are NAIVE (independent-observation). Every tree in a ",
      "fire shares one dNBR value, so across %d fires the reported p-values are ",
      "anticonservative - do not read them as evidence for B2."), n_clusters))
  }

  coef_df <- tidy_fit %>%
    rename(std_error = std.error, p_value = p.value) %>%
    select(term, estimate, std_error, p_value) %>%
    mutate(scf_parameter = dplyr::case_when(
      term == "(Intercept)"   ~ "CohortMortalityB0",
      term == "BarkThickness" ~ "CohortMortalityB1",
      term == "dNBR"          ~ "CohortMortalityB2",
      TRUE                    ~ paste0("CohortMortalityB", dplyr::row_number())
    ))

  message(sprintf(
    "fit_cohort_mortality: n=%d (dead=%d, alive=%d)  B0=%.3f  B1=%.4f  B2=%.4f",
    nrow(dat), sum(dat$dead), sum(dat$dead == 0),
    coef(fit)["(Intercept)"],
    coef(fit)["BarkThickness"],
    coef(fit)["dNBR"]
  ))

  # ---- Sign guard ------------------------------------------------------------
  # Logit link: thicker bark must REDUCE mortality (B1 < 0) and higher site
  # severity must RAISE it (B2 > 0).
  b1 <- coef(fit)[["BarkThickness"]]
  b2 <- coef(fit)[["dNBR"]]
  sign_warnings <- character(0)
  if (is.finite(b1) && b1 >= 0)
    sign_warnings <- c(sign_warnings, sprintf(
      "CohortMortalityB1 = %+.4f is not negative - implies thicker bark increases mortality.", b1))
  if (is.finite(b2) && b2 <= 0)
    sign_warnings <- c(sign_warnings, sprintf(
      "CohortMortalityB2 = %+.5f is not positive - implies higher site severity reduces mortality.", b2))
  if (length(sign_warnings) > 0)
    warning(paste(c("Cohort mortality sign check:", sign_warnings), collapse = "\n  - "))

  # How much distinct severity information actually exists
  n_dnbr_levels <- length(unique(dat$dNBR))
  if (n_dnbr_levels < 5)
    warning(sprintf(paste0(
      "Only %d distinct dNBR values across %d fires. CohortMortalityB2 is being ",
      "estimated from almost no severity contrast and should be treated as ",
      "indicative only."), n_dnbr_levels, n_clusters))

  list(
    model             = fit,
    coef_df           = coef_df,
    data              = dat,
    n                 = nrow(dat),
    n_dead            = sum(dat$dead),
    n_alive           = sum(dat$dead == 0),
    n_fires           = n_clusters,
    n_dnbr_levels     = n_dnbr_levels,
    cluster_robust    = cluster_ok,
    sign_warnings     = sign_warnings,
    null_deviance     = fit$null.deviance,
    residual_deviance = fit$deviance
  )
}
