# =============================================================================
# species_helpers.R
# FIA-based bark parameter estimation for SCF Fire_Spp_Table.csv
#
# SCF Eq. 10:  BarkThickness = (MaxBarkThickness_s * Age) / (Age + AgeDBH_s)
#   MaxBarkThickness_s : asymptotic bark thickness (cm) as DBH → max
#   AgeDBH_s           : age at which bark thickness = MaxBarkThickness/2
#
# FIA TREE table fields used:
#   SPCD       : species code (integer)
#   DIA        : DBH outside bark (inches)
#   BARK_THICK : bark thickness, tenths of an inch, one-side measurement
#                (e.g. value 6 = 0.6 inches one side)
#   TOTAGE     : total tree age at breast height (present for increment-cored trees)
#   STATUSCD   : 1 = live, 2 = dead  →  keep only 1
# =============================================================================

pacman::p_load(dplyr, tibble, readr, stringr, purrr)

# rFIA is available on CRAN; install once if missing
if (!requireNamespace("rFIA", quietly = TRUE)) {
  message("Installing rFIA from CRAN...")
  install.packages("rFIA")
}


# =============================================================================
# load_species_lookup
# Read user-supplied LANDIS → FIA SPCD lookup CSV.
# Expected columns (case-insensitive): SpeciesName, FIA_SPCD, SCIENTIFIC_NAME
# Returns a tibble with standardised column names.
# =============================================================================
load_species_lookup <- function(csv_path) {
  if (!file.exists(csv_path)) stop("Lookup CSV not found: ", csv_path)

  df <- read_csv(csv_path, show_col_types = FALSE, name_repair = "minimal")

  # Normalise column names to lower-case for robust matching
  names(df) <- tolower(names(df))

  required <- c("speciesname", "fia_spcd")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Lookup CSV is missing required columns: ",
         paste(missing, collapse = ", "),
         "\nExpected: SpeciesName, FIA_SPCD, SCIENTIFIC_NAME")
  }

  df <- df %>%
    rename(
      SpeciesName     = speciesname,
      FIA_SPCD        = fia_spcd
    ) %>%
    mutate(
      FIA_SPCD        = as.integer(FIA_SPCD),
      SCIENTIFIC_NAME = if ("scientific_name" %in% names(.)) scientific_name else NA_character_
    ) %>%
    select(SpeciesName, FIA_SPCD, SCIENTIFIC_NAME) %>%
    filter(!is.na(SpeciesName), !is.na(FIA_SPCD))

  df
}


# =============================================================================
# load_fia_trees
# Load FIA TREE table for the requested states.
#
# Strategy (in order):
#   1. For each state, look for {STATE}_TREE.csv in fia_dir (if provided).
#      Only use files that exist AND are non-empty (> 0 bytes).
#   2. For any states with missing/empty files, download via rFIA::getFIA()
#      into fia_dir (or a session temp folder if fia_dir is NULL).
#      getFIA() writes {STATE}_TREE.csv; subsequent calls skip re-download.
#   3. Read all TREE CSVs directly with read_csv() — avoids rFIA::readFIA()
#      which fails with "OneIndex" error on empty files.
#
# Returns a tibble filtered to live trees with bark & DBH measurements,
# with derived columns dia_cm and bark_cm.
# =============================================================================
load_fia_trees <- function(states, fia_dir = NULL) {

  if (!requireNamespace("rFIA", quietly = TRUE)) {
    stop("The 'rFIA' package is required. Install with: install.packages('rFIA')")
  }

  # Resolve the directory we will read from / download into
  if (!is.null(fia_dir) && nchar(trimws(fia_dir)) > 0) {
    if (!dir.exists(fia_dir)) dir.create(fia_dir, recursive = TRUE)
    work_dir <- fia_dir
  } else {
    work_dir <- file.path(tempdir(), "fia_download")
    dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Check which states already have non-empty TREE CSVs locally
  tree_paths     <- file.path(work_dir, paste0(states, "_TREE.csv"))
  names(tree_paths) <- states
  needs_download <- states[
    !file.exists(tree_paths) | file.info(tree_paths)$size == 0
  ]

  # Download missing / empty states via rFIA::getFIA()
  if (length(needs_download) > 0) {
    message("rFIA: downloading TREE table for: ",
            paste(needs_download, collapse = ", "), " → ", work_dir)
    tryCatch(
      rFIA::getFIA(
        states = needs_download,
        tables = "TREE",
        dir    = work_dir,
        nCores = 1
      ),
      error = function(e) {
        stop("rFIA::getFIA() failed for states [",
             paste(needs_download, collapse = ", "), "]: ",
             conditionMessage(e))
      }
    )
  }

  # Verify files now exist and are non-empty after download attempt
  still_missing <- states[
    !file.exists(tree_paths) | file.info(tree_paths)$size == 0
  ]
  if (length(still_missing) > 0) {
    stop("FIA TREE CSV is still empty or missing after download for: ",
         paste(still_missing, collapse = ", "),
         "\nExpected files:\n  ",
         paste(tree_paths[still_missing], collapse = "\n  "))
  }

  # Read CSVs directly — much more robust than rFIA::readFIA()
  message("Reading FIA TREE CSVs for: ", paste(states, collapse = ", "))
  trees_raw <- map_dfr(tree_paths, function(p) {
    read_csv(p, show_col_types = FALSE,
             col_types = cols(.default = col_guess()))
  })

  # Normalise column names to upper-case (FIA CSVs can vary in case)
  names(trees_raw) <- toupper(names(trees_raw))

  # Validate required columns
  req_cols <- c("SPCD", "DIA", "BARK_THICK", "STATUSCD")
  missing_cols <- setdiff(req_cols, names(trees_raw))
  if (length(missing_cols) > 0) {
    stop("FIA TREE CSV is missing expected columns: ",
         paste(missing_cols, collapse = ", "))
  }

  # Filter and convert units
  trees <- trees_raw %>%
    filter(
      STATUSCD == 1L,       # live trees only
      !is.na(DIA),
      !is.na(BARK_THICK),
      DIA > 0,
      BARK_THICK > 0
    ) %>%
    select(SPCD, DIA, BARK_THICK, any_of("TOTAGE")) %>%
    mutate(
      SPCD     = as.integer(SPCD),
      # DIA        : inches → cm
      dia_cm   = DIA * 2.54,
      # BARK_THICK : tenths-of-an-inch (one side) → cm
      #              value 6 = 0.6 in × 2.54 = 1.524 cm
      bark_cm  = (BARK_THICK / 10) * 2.54
    )

  if (nrow(trees) == 0) {
    stop("No live trees with bark/DBH measurements found for states: ",
         paste(states, collapse = ", "))
  }

  message(sprintf("FIA TREE data loaded: %s live trees across %s states.",
                  format(nrow(trees), big.mark = ","), length(states)))
  trees
}


# =============================================================================
# estimate_bark_params_from_fia
# Map LANDIS species codes → FIA SPCDs via the user lookup, then compute
# MaximumBarkThickness and AgeDBH from the FIA TREE data.
#
# lookup    : tibble from load_species_lookup() — cols SpeciesName, FIA_SPCD,
#             SCIENTIFIC_NAME
# Returns a tibble: SpeciesCode, AgeDBH, MaximumBarkThickness,
#                   ScientificName, FIA_SPCD, Note
# =============================================================================
estimate_bark_params_from_fia <- function(species_codes, fia_trees, lookup) {

  # Map LANDIS codes to FIA SPCD via lookup
  mapped <- tibble(SpeciesCode = species_codes) %>%
    left_join(lookup, by = c("SpeciesCode" = "SpeciesName"))

  # Unique SPCDs with data
  spcd_vec  <- unique(na.omit(mapped$FIA_SPCD))
  trees_sub <- fia_trees %>% filter(SPCD %in% spcd_vec)

  # --- Per-species bark summary ---------------------------------------------
  bark_summary <- trees_sub %>%
    group_by(SPCD) %>%
    summarise(
      n_trees = n(),
      # 95th-percentile bark thickness across all trees
      p95_bark_cm = quantile(bark_cm, 0.95, na.rm = TRUE),
      # Mean bark thickness for the top-5% largest trees by DBH
      bark_at_large_dia = {
        q95_dia <- quantile(dia_cm, 0.95, na.rm = TRUE)
        mean(bark_cm[dia_cm >= q95_dia], na.rm = TRUE)
      },
      max_dia_cm = max(dia_cm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # MaximumBarkThickness: take the larger of the two estimates (cm, 2 dp)
      MaximumBarkThickness = round(
        pmax(p95_bark_cm, bark_at_large_dia, na.rm = TRUE), 2
      )
    )

  # --- AgeDBH via Michaelis-Menten fit on TOTAGE ----------------------------
  # BarkThickness = MaxBT × Age / (Age + AgeDBH)
  # Only attempted when TOTAGE is present and ≥ 15 aged trees per species
  age_fits <- if ("TOTAGE" %in% names(fia_trees)) {
    trees_sub %>%
      filter(!is.na(TOTAGE), TOTAGE > 5) %>%
      group_by(SPCD) %>%
      filter(n() >= 15) %>%
      summarise(
        n_aged  = n(),
        age_dbh = {
          max_bt <- quantile(bark_cm, 0.95, na.rm = TRUE)
          fit <- tryCatch(
            nls(
              bark_cm ~ max_bt * TOTAGE / (TOTAGE + theta),
              data    = cur_data(),
              start   = list(theta = 60),
              control = nls.control(maxiter = 200, warnOnly = TRUE)
            ),
            error = function(e) NULL
          )
          if (!is.null(fit)) max(1, round(coef(fit)[["theta"]])) else NA_real_
        },
        .groups = "drop"
      )
  } else {
    tibble(SPCD = integer(0), n_aged = integer(0), age_dbh = numeric(0))
  }

  # --- Combine, apply bark-class heuristic where nls failed -----------------
  results <- bark_summary %>%
    left_join(age_fits %>% select(SPCD, n_aged, age_dbh), by = "SPCD") %>%
    mutate(
      AgeDBH = case_when(
        !is.na(age_dbh) & age_dbh > 0  ~ as.integer(round(age_dbh)),
        MaximumBarkThickness >= 30      ~ 80L,   # thick-barked → slow
        MaximumBarkThickness >= 15      ~ 60L,   # moderate
        TRUE                            ~ 40L    # thin-barked → fast
      ),
      fia_note = paste0(
        "n=", n_trees, " trees | MaxBT from 95th pct + large-DBH avg",
        if_else(!is.na(n_aged),
                paste0(" | AgeDBH nls fit (n_aged=", n_aged, ")"),
                " | AgeDBH from bark-thickness class heuristic")
      )
    )

  # --- Build final output tibble --------------------------------------------
  out <- mapped %>%
    left_join(
      results %>% select(SPCD, AgeDBH, MaximumBarkThickness, fia_note),
      by = c("FIA_SPCD" = "SPCD")
    ) %>%
    mutate(
      MaximumBarkThickness = if_else(!is.na(MaximumBarkThickness),
                                     MaximumBarkThickness, 10.0),
      AgeDBH = if_else(!is.na(AgeDBH), AgeDBH, 100L),
      Note = case_when(
        is.na(FIA_SPCD)  ~ "Not in lookup CSV — defaults used (AgeDBH=100, MaxBT=10 cm)",
        is.na(fia_note)  ~ paste0("SPCD ", FIA_SPCD,
                                   " — no FIA records in selected states; defaults used"),
        TRUE             ~ paste0("SPCD ", FIA_SPCD, " | ", fia_note)
      ),
      ScientificName = coalesce(SCIENTIFIC_NAME, "")
    ) %>%
    select(SpeciesCode, AgeDBH, MaximumBarkThickness,
           ScientificName, FIA_SPCD, Note)

  out
}


# =============================================================================
# build_default_species_table
# Populate a stub table (default AgeDBH=100, MaxBT=10) for the given species
# codes, joining scientific names from the lookup where available.
# Used to pre-fill the editable DT before FIA data are fetched.
# =============================================================================
build_default_species_table <- function(species_codes, lookup = NULL) {
  tbl <- tibble(SpeciesCode = species_codes) %>%
    mutate(AgeDBH = 100L, MaximumBarkThickness = 10.0)

  if (!is.null(lookup) && nrow(lookup) > 0) {
    tbl <- tbl %>%
      left_join(lookup %>% select(SpeciesName, FIA_SPCD, SCIENTIFIC_NAME),
                by = c("SpeciesCode" = "SpeciesName")) %>%
      mutate(
        ScientificName = coalesce(SCIENTIFIC_NAME, ""),
        Note = if_else(
          !is.na(FIA_SPCD),
          paste0("Default — SPCD ", FIA_SPCD,
                 "; click 'Fetch FIA Parameters' to estimate from data"),
          "Not in lookup CSV — enter values manually"
        )
      ) %>%
      select(SpeciesCode, AgeDBH, MaximumBarkThickness,
             ScientificName, FIA_SPCD, Note)
  } else {
    tbl <- tbl %>%
      mutate(ScientificName = "", FIA_SPCD = NA_integer_,
             Note = "Default — load lookup CSV and fetch FIA data")
  }

  tbl
}


# =============================================================================
# bark_curve_data
# Generate SCF Eq. 10 curve points for all species in spp_table.
# Returns a long tibble for ggplot: SpeciesCode, Age, BarkThickness_cm.
# =============================================================================
bark_curve_data <- function(spp_table, max_age = 500) {
  ages <- seq(0, max_age, by = 5)
  map_dfr(seq_len(nrow(spp_table)), function(i) {
    row <- spp_table[i, ]
    tibble(
      SpeciesCode          = row$SpeciesCode,
      Age                  = ages,
      BarkThickness_cm     = (row$MaximumBarkThickness * ages) /
                               (ages + row$AgeDBH),
      MaximumBarkThickness = row$MaximumBarkThickness
    )
  })
}
