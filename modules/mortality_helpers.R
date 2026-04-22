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
                                   min_dnbr = 100, max_dnbr = 2000,
                                   progress_fn = NULL) {

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
  n_loaded <- 0L

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
    cat(sprintf("    fire_year=%d\n", fire_year))

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
    if (raw_nonNA > 0 && max(raw_unique, na.rm = TRUE) <= 6) {
      cat(sprintf("    reclassifying severity classes to dNBR midpoints\n"))
      r <- subst(r,
                 from = c(0,  1,   2,   3,   4,   5,  6),
                 to   = c(NA, NA, 185, 355, 600,  NA, NA))
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
      dNBR      = vals,
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
# fetch_landfire_cbh
#
# Purpose: Download LANDFIRE LF2024 Canopy Base Height (CBH) via REST API,
#   convert from decimetres to metres, compute a ladder fuel index (0-1), and
#   resample to spread_template.
#
# Args:
#   cal_vect        SpatVector of calibration boundary (any CRS)
#   working_crs     CRS string
#   spread_template SpatRaster defining the output grid
#   progress_fn     optional function(msg)
#
# Returns: SpatRaster on spread_template grid, values 0-1 (high = more ladder fuels)
# -----------------------------------------------------------------------------

fetch_landfire_cbh <- function(cal_vect, working_crs, spread_template,
                                progress_fn = NULL) {

  # ---- Bounding box in WGS84 ------------------------------------------------
  cal_wgs84 <- project(cal_vect, "EPSG:4326")
  bb        <- as.vector(ext(cal_wgs84))   # xmin, xmax, ymin, ymax

  deg_per_30m <- 30 / 111320
  nx <- min(as.integer(ceiling((bb["xmax"] - bb["xmin"]) / deg_per_30m)), 4096L)
  ny <- min(as.integer(ceiling((bb["ymax"] - bb["ymin"]) / deg_per_30m)), 4096L)
  message(sprintf(
    "fetch_landfire_cbh: requesting CBH %d x %d pixels over bbox [%.4f %.4f %.4f %.4f]",
    nx, ny, bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"]
  ))

  base_url <- paste0(
    "https://lfps.usgs.gov/arcgis/rest/services/",
    "Landfire_LF2024/LF2024_CBH_CONUS/ImageServer/exportImage"
  )

  if (!is.null(progress_fn)) progress_fn("Downloading LANDFIRE CBH from USGS server...")

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
      stop("LANDFIRE CBH HTTP request failed: ", conditionMessage(e))
  )

  if (httr::http_error(resp))
    stop("LANDFIRE CBH REST API returned HTTP ", httr::status_code(resp),
         ". Check internet connection or try again later.")

  # ---- Read, validate, convert units -----------------------------------------
  cbh_r <- tryCatch(
    rast(tmp_tif),
    error = function(e)
      stop("Could not read LANDFIRE CBH response as raster: ", conditionMessage(e))
  )

  raw_vals <- values(cbh_r, mat = FALSE)
  message(sprintf("LANDFIRE CBH raw pixel range: [%.0f, %.0f] (decimetres)",
                  min(raw_vals, na.rm = TRUE), max(raw_vals, na.rm = TRUE)))

  # CBH in LANDFIRE = tenths of a metre (decimetres); 0 = non-forested
  cbh_m <- cbh_r / 10.0
  cbh_m[cbh_m < 0] <- NA   # -9999 nodata

  if (!is.null(progress_fn)) progress_fn("Computing ladder fuel index from CBH...")

  # Low CBH -> high ladder fuel risk
  # LadderFuels = 1 / (CBH_m + 0.5) then normalize 0-1 by max
  lf_r     <- 1 / (cbh_m + 0.5)
  lf_vals  <- values(lf_r, mat = FALSE)
  lf_max   <- max(lf_vals, na.rm = TRUE)
  if (!is.finite(lf_max) || lf_max <= 0)
    stop("Ladder fuel index calculation produced no valid values.")
  lf_r <- lf_r / lf_max

  # ---- Reproject and resample to spread grid ----------------------------------
  if (!is.null(progress_fn)) progress_fn("Reprojecting and resampling CBH to working CRS...")
  lf_proj <- project(lf_r, working_crs, method = "bilinear")
  resample(lf_proj, spread_template, method = "bilinear")
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
#   ladder_fuels_r  SpatRaster of ladder fuel index (0-1) or NULL
#   working_crs     CRS string (pixel_df x,y are in this CRS)
#
# Returns: pixel_df with additional columns Clay, ET_mm, CWD_mm, EWS_kmh,
#   FineFuels, LadderFuels; rows with any NA predictor dropped.
# -----------------------------------------------------------------------------

extract_site_predictors <- function(pixel_df, clay_r, et_r, cwd_r,
                                     wind_annual_tbl, slope_r, aspect_r,
                                     fine_fuels_r, ladder_fuels_r,
                                     working_crs = NULL) {

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
  message("extract_site_predictors: computing EWS...")
  if (!is.null(wind_annual_tbl)) {
    pixel_df <- pixel_df %>%
      left_join(wind_annual_tbl %>% select(year, mean_ws_kmh, mean_wd_deg),
                by = c("fire_year" = "year"))

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

  } else {
    message("extract_site_predictors: no wind data, EWS set to NA")
    pixel_df$EWS_kmh <- NA_real_
  }

  # ---- Fine fuels ------------------------------------------------------------
  if (!is.null(fine_fuels_r)) {
    message("extract_site_predictors: extracting fine fuels...")
    pixel_df$FineFuels <- terra::extract(fine_fuels_r, pv)[, 2]
  } else {
    message("extract_site_predictors: fine fuels not available, using placeholder 1.0")
    pixel_df$FineFuels <- 1.0
  }

  # ---- Ladder fuels ----------------------------------------------------------
  if (!is.null(ladder_fuels_r)) {
    message("extract_site_predictors: extracting ladder fuels...")
    pixel_df$LadderFuels <- terra::extract(ladder_fuels_r, pv)[, 2]
  } else {
    message("extract_site_predictors: ladder fuels not available, using placeholder 0.5")
    pixel_df$LadderFuels <- 0.5
  }

  # ---- Drop rows with any NA predictor ---------------------------------------
  pred_cols <- c("Clay", "ET_mm", "CWD_mm", "EWS_kmh", "FineFuels", "LadderFuels")
  # ET_mm and CWD_mm may be all-NA if not fetched; keep them but don't require them
  required_preds <- intersect(c("Clay", "EWS_kmh", "FineFuels", "LadderFuels"),
                               names(pixel_df))
  n_before <- nrow(pixel_df)
  pixel_df <- pixel_df %>%
    filter(if_all(all_of(required_preds), is.finite))
  message(sprintf("extract_site_predictors: %d of %d rows retained after NA drop",
                  nrow(pixel_df), n_before))

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
  available_preds <- candidate_preds[vapply(candidate_preds, function(v) {
    v %in% names(pixel_df) && sum(is.finite(pixel_df[[v]])) > 10
  }, logical(1))]

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

  list(
    model             = fit,
    coef_df           = coef_df,
    data              = dat,
    null_deviance     = fit$null.deviance,
    residual_deviance = fit$deviance,
    n                 = nrow(dat),
    link              = fit$family$link,
    available_preds   = available_preds
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

write_mortality_snippet <- function(site_fit, cohort_fit = NULL) {

  now   <- format(Sys.time(), "%Y-%m-%d %H:%M")
  lines <- c(
    "# =============================================================================",
    paste0("# SCF Mortality Parameters — generated ", now),
    "# =============================================================================",
    "",
    "# ---------------------------------------------------------------------------",
    "# Site Mortality (Eq. 7)",
    "#   dNBR = 1 / (B0 + B1*Clay + B2*ET + B3*EWS + B4*CWD + B5*FineFuels + B6*LadderFuels)",
    "#   Model: Gamma GLM with inverse link",
    paste0("#   Link used: ", site_fit$link),
    paste0("#   N pixels: ", scales::comma(site_fit$n)),
    paste0("#   Predictors fitted: ", paste(site_fit$available_preds, collapse = ", ")),
    sprintf("#   Residual deviance: %.4f  (null: %.4f)",
            site_fit$residual_deviance, site_fit$null_deviance),
    "# ---------------------------------------------------------------------------",
    ""
  )

  for (i in seq_len(nrow(site_fit$coef_df))) {
    row <- site_fit$coef_df[i, ]
    unit_note <- switch(row$term,
      "(Intercept)" = "# intercept",
      "Clay"        = "# clay fraction 0-1",
      "ET_mm"       = "# potential ET mm/yr (previous water year; SCF PotentialEvapotranspiration)",
      "EWS_kmh"     = "# effective wind speed km/h (Nelson 2002; SCF User Guide §2.51: climate library uses km/h)",
      "CWD_mm"      = "# climatic water deficit mm/yr (previous water year)",
      "FineFuels"   = "# fine fuel index 0-1 (LANDFIRE FBFM40)",
      "LadderFuels" = "# ladder fuel index 0-1 (1 / (CBH_m + 0.5), normalized)",
      ""
    )
    # Mark zero-padded (not fitted) parameters clearly
    if (!isTRUE(row$fitted) && row$term != "(Intercept)") {
      unit_note <- paste(unit_note, "  [not fitted — set to 0]")
    }
    lines <- c(lines,
               sprintf("%-30s  %+.8f  %s",
                       paste0(row$scf_parameter, " ="),
                       row$estimate,
                       unit_note))
  }

  if (!is.null(cohort_fit)) {
    lines <- c(
      lines, "",
      "# ---------------------------------------------------------------------------",
      "# Cohort Mortality (Eq. 8/9)",
      "#   P(mortality) = exp(P)/(1+P)  where P = B0 + B1*BarkThickness + B2*dNBR",
      "#   Model: Binomial GLM with logit link",
      paste0("#   N trees: ", scales::comma(cohort_fit$n)),
      "# ---------------------------------------------------------------------------",
      ""
    )
    for (i in seq_len(nrow(cohort_fit$coef_df))) {
      row <- cohort_fit$coef_df[i, ]
      unit_note <- switch(row$term,
        "(Intercept)"  = "# intercept",
        "BarkThickness" = "# bark thickness cm (from species table)",
        "dNBR"         = "# site dNBR (from MTBS)",
        ""
      )
      lines <- c(lines,
                 sprintf("%-30s  %+.8f  %s",
                         paste0(row$scf_parameter, " ="),
                         row$estimate,
                         unit_note))
    }
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
# Mortality status coding (Cansler et al. 2020 FTM standard):
#   yr1status 0 = not returned to (unobserved — NOT confirmed dead; dropped)
#   yr1status 1 = alive at first post-fire survey  -> dead = 0
#   yr1status 2 = dead at first post-fire survey   -> dead = 1
#   Trees with NA, 0, or other yr1status values are dropped.
#   NOTE: Some calibration regions contain only yr1status=0 and 1 (no confirmed
#   dead trees). In that case a warning is issued and fit_cohort_mortality() will
#   fail with a meaningful message.
#
# Return value: list(tree_df, yr1status_counts)
#   tree_df          — standard tibble suitable for fit_cohort_mortality()
#   yr1status_counts — tibble(yr1status, n) of raw counts before filtering,
#                      used by diagnostic plot functions
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
  # Cansler et al. (2020) FTM standard coding:
  #   yr1status 0 = not returned to (unobserved — NOT confirmed dead)
  #   yr1status 1 = alive   -> dead = 0
  #   yr1status 2 = dead    -> dead = 1
  if (!"yr1status" %in% names(trees_cal))
    stop("yr1status column not found in FTM_trees.csv.")

  # Save raw distribution BEFORE filtering for diagnostic plots
  yr1status_counts <- trees_cal %>%
    filter(!is.na(yr1status)) %>%
    count(yr1status, name = "n") %>%
    arrange(yr1status)

  message("  yr1status raw distribution (before filtering):")
  for (k in seq_len(nrow(yr1status_counts)))
    message(sprintf("    yr1status=%d  n=%d", yr1status_counts$yr1status[k], yr1status_counts$n[k]))

  # Filter to confirmed alive (1) or dead (2) only; drop yr1status=0 "not returned to"
  trees_cal <- trees_cal %>%
    filter(!is.na(yr1status), yr1status %in% c(1L, 2L)) %>%
    mutate(dead = as.integer(yr1status == 2L))

  # Warn clearly if no confirmed dead trees are present in calibration region
  n_dead_initial <- sum(trees_cal$dead)
  if (n_dead_initial == 0) {
    warning(paste0(
      "No yr1status=2 (confirmed dead) trees found within/near cal_vect. ",
      "The cohort mortality GLM requires both alive AND dead trees. ",
      "Cohort mortality fitting will fail. ",
      "Possible cause: this calibration region's FTM subset contains only ",
      "yr1status=0 ('not returned to') and yr1status=1 ('alive') trees — ",
      "confirmed dead trees were not re-sampled in this region."
    ))
  }

  # ---- Filter by DBH ---------------------------------------------------------
  trees_cal <- trees_cal %>%
    filter(!is.na(DBH_cm), DBH_cm >= min_dbh_cm)

  message(sprintf("  %d trees with yr1status 1 or 2 (dead=%d, alive=%d), DBH >= %.1f cm",
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
#   area. This is the most diagnostic figure: if yr1status=2 (confirmed dead)
#   is absent, cohort mortality cannot be estimated.
#
# Args:
#   yr1status_counts — tibble(yr1status, n) from load_ftm_cohort_data()$yr1status_counts
#
# Returns: ggplot object
# -----------------------------------------------------------------------------

plot_yr1status_distribution <- function(yr1status_counts) {

  status_colors <- c("0" = "#aaaaaa", "1" = "#2c7bb6", "2" = "#d7191c")

  df <- yr1status_counts %>%
    mutate(
      key   = as.character(yr1status),
      label = dplyr::case_when(
        key == "0" ~ "0: Not returned to\n(unobserved \u2014 NOT dead)",
        key == "1" ~ "1: Alive",
        key == "2" ~ "2: Dead (confirmed)",
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
      subtitle = paste0("yr1status=2 (confirmed dead) MUST be present for cohort mortality GLM.\n",
                        "yr1status=0 ('not returned to') is NOT dead — these are excluded."),
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

  df <- tree_df %>%
    filter(is.finite(dead), is.finite(BarkThickness)) %>%
    mutate(status = factor(dead, levels = c(0L, 1L),
                           labels = c(
                             sprintf("Alive (yr1=1)\nn=%s", scales::comma(n_alive)),
                             sprintf("Dead (yr1=2)\nn=%s",  scales::comma(n_dead))
                           )))

  fill_vals <- setNames(
    c("#2c7bb6", "#d7191c"),
    c(sprintf("Alive (yr1=1)\nn=%s", scales::comma(n_alive)),
      sprintf("Dead (yr1=2)\nn=%s",  scales::comma(n_dead)))
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
        "No confirmed dead trees (yr1status=2) in calibration region \u2014 check yr1status distribution"
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
#   B1 -> CohortMortalityB1   (expect negative: thicker bark protects)
#   B2 -> CohortMortalityB2   (expect positive: higher dNBR kills)
#
# Args:
#   tree_df  tibble from load_ftm_cohort_data() with columns:
#            dead (0/1), BarkThickness (cm), dNBR
#
# Returns: list(model, coef_df, data, n, n_dead, n_alive,
#               null_deviance, residual_deviance)
# -----------------------------------------------------------------------------

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

  tidy_fit <- broom::tidy(fit)

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

  list(
    model             = fit,
    coef_df           = coef_df,
    data              = dat,
    n                 = nrow(dat),
    n_dead            = sum(dat$dead),
    n_alive           = sum(dat$dead == 0),
    null_deviance     = fit$null.deviance,
    residual_deviance = fit$deviance
  )
}
