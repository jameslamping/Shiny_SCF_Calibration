# =============================================================================
# species_helpers.R
# FIA-based bark parameter estimation for SCF Fire_Spp_Table.csv
#
# SCF Eq. 10:  BarkThickness = (MaximumBarkThickness * Age) / (Age + AgeDBH)
#   MaximumBarkThickness : double-sided bark thickness (cm) at maximum observed DBH
#   AgeDBH               : age (yr) at which bark thickness = MaximumBarkThickness / 2
#                          Derived here as the half-saturation constant from a
#                          Michaelis-Menten DBH~Age growth curve:
#                          DBH_cm = Dmax * Age / (Age + AgeDBH)
#
# FIA TREE fields used:
#   SPCD      : numeric species code
#   DIA       : DBH outside bark (inches) — converted to cm
#   TOTAGE    : total tree age (preferred; sparse — cored trees only)
#   BHAGE     : age at breast height (more common; missing years-to-BH)
#   STATUSCD  : 1 = live only
#
# MaximumBarkThickness is estimated as:
#   2 × bark_ratio × P95_DBH_cm
# where bark_ratio (one-side BT / DBH, dimensionless) comes from the built-in
# Cansler / Ryan & Reinhardt literature table, and P95_DBH_cm is the 95th-
# percentile DBH for the species across FIA plots in the selected states.
# =============================================================================

pacman::p_load(dplyr, tibble, readr, stringr, purrr)

if (!requireNamespace("rFIA", quietly = TRUE)) {
  message("Installing rFIA from CRAN...")
  install.packages("rFIA")
}

# =============================================================================
# Built-in bark thickness ratio table
# bark_ratio = one-side bark thickness (cm) / DBH (cm)
# MaximumBarkThickness (double-sided, cm) = 2 × bark_ratio × DBH_cm
#
# Primary sources:
#   Ryan & Reinhardt (1988) Predicting postfire mortality of seven western conifers.
#   Can. J. For. Res. 18:1291–1297.  [R&R88]
#   Cansler et al. (2020) Fire Tree Mortality Database.
#   USFS RDS-2020-0002.  [C20]
#   Silvics of North America (Burns & Honkala 1990).  [SNA]
# =============================================================================
bark_ratio_table <- tribble(
  ~fia_spcd, ~landis_code,               ~bark_ratio,  ~source,
  # --- Firs (Abies) — thin-barked -------------------------------------------
  11L,       "AbiesAmabilis",             0.017,        "SNA / C20",
  17L,       "AbiesGrandis",              0.022,        "R&R88 / C20",
  19L,       "AbiesLasiocarpa",           0.019,        "R&R88 / C20",
  22L,       "AbiesProcera",              0.022,        "SNA / C20",
  # --- Hardwoods ---------------------------------------------------------------
  312L,      "AcerMacrophyllum",          0.018,        "SNA",
  351L,      "AlnusRubra",                0.012,        "SNA",
  361L,      "ArbutusMenziesii",          0.025,        "SNA",
  375L,      "BetulaPapyrifera",          0.015,        "SNA",
  # --- Yellow-cedar -----------------------------------------------------------
  42L,       "ChamaecyparisNootkatensis", 0.025,        "SNA / C20",
  # --- Larch (Larix) — moderate to thick bark ---------------------------------
  72L,       "LarixLyallii",              0.030,        "SNA",
  73L,       "LarixOccidentalis",         0.035,        "R&R88 / C20",
  # --- Spruce (Picea) ---------------------------------------------------------
  93L,       "PiceaEngelmannii",          0.025,        "R&R88 / C20",
  98L,       "PiceaSitchensis",           0.020,        "SNA / C20",
  # --- Pine (Pinus) -----------------------------------------------------------
  101L,      "PinusAlbicaulis",           0.028,        "SNA / C20",
  108L,      "PinusContorta",             0.021,        "R&R88 / C20",
  119L,      "PinusMonticola",            0.038,        "SNA / C20",
  122L,      "PinusPonderosa",            0.074,        "R&R88 / C20",  # thick-barked
  # --- Poplar / Aspen — thin bark ----------------------------------------------
  741L,      "PopulusBalsamifera",        0.015,        "SNA",
  747L,      "PopulusBalsamifera",        0.015,        "SNA",  # trichocarpa (WA)
  746L,      "PopulusTremuloides",        0.012,        "SNA",
  # --- Douglas-fir — thick bark ------------------------------------------------
  202L,      "PseudotsugaMenziesii",      0.041,        "R&R88 / C20",
  # --- Oak, Yew, Cedars, Hemlocks ----------------------------------------------
  815L,      "QuercusGarryana",           0.040,        "SNA",
  231L,      "TaxusBrevifolia",           0.018,        "SNA",
  242L,      "ThujaPlicata",              0.023,        "SNA / C20",
  263L,      "TsugaHeterophylla",         0.021,        "SNA / C20",
  264L,      "TsugaMertensiana",          0.018,        "R&R88 / C20"
)

# Years-to-breast-height correction for BHAGE → approximate TOTAGE
# Used only when TOTAGE is absent. Sources: Silvics of North America.
years_to_bh <- tribble(
  ~fia_spcd, ~ytbh,
  11L,  5L,   # Abies amabilis — slow
  17L,  4L,   # Abies grandis
  19L,  5L,   # Abies lasiocarpa — slow, high elevation
  22L,  5L,   # Abies procera
  42L,  8L,   # Chamaecyparis nootkatensis — very slow
  72L, 10L,   # Larix lyallii — extremely slow, treeline
  73L,  4L,   # Larix occidentalis
  93L,  5L,   # Picea engelmannii
  98L,  4L,   # Picea sitchensis
  101L,10L,   # Pinus albicaulis — very slow, treeline
  108L, 4L,   # Pinus contorta
  119L, 4L,   # Pinus monticola
  122L, 3L,   # Pinus ponderosa
  202L, 3L,   # Pseudotsuga menziesii
  231L, 5L,   # Taxus brevifolia — slow
  242L, 4L,   # Thuja plicata
  263L, 4L,   # Tsuga heterophylla
  264L, 5L,   # Tsuga mertensiana
  312L, 2L,   # Acer macrophyllum — fast
  351L, 2L,   # Alnus rubra — fast
  361L, 3L,   # Arbutus menziesii
  375L, 2L,   # Betula papyrifera — fast
  741L, 2L,   # Populus balsamifera — fast
  747L, 2L,   # Populus trichocarpa — fast
  746L, 2L    # Populus tremuloides — fast
)


# =============================================================================
# load_species_lookup
# Read user-supplied LANDIS → FIA SPCD lookup CSV.
# Expected columns: SpeciesName, FIA_SPCD, SCIENTIFIC_NAME (case-insensitive).
# Returns standardised tibble and warns about SPCDs that differ from the
# built-in bark ratio table (common FIA reference codes).
# =============================================================================
load_species_lookup <- function(csv_path) {
  if (!file.exists(csv_path)) stop("Lookup CSV not found: ", csv_path)

  df <- read_csv(csv_path, show_col_types = FALSE, name_repair = "minimal")
  names(df) <- tolower(names(df))

  req <- c("speciesname", "fia_spcd")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) {
    stop("Lookup CSV missing columns: ", paste(miss, collapse = ", "),
         "\nExpected: SpeciesName, FIA_SPCD, SCIENTIFIC_NAME")
  }

  df <- df %>%
    rename(SpeciesName = speciesname, FIA_SPCD = fia_spcd) %>%
    mutate(
      FIA_SPCD        = as.integer(FIA_SPCD),
      SCIENTIFIC_NAME = if ("scientific_name" %in% names(.))
                          scientific_name else NA_character_
    ) %>%
    select(SpeciesName, FIA_SPCD, SCIENTIFIC_NAME) %>%
    filter(!is.na(SpeciesName), !is.na(FIA_SPCD))

  # Check for discrepancies against bark ratio reference table
  ref <- bark_ratio_table %>%
    distinct(fia_spcd, landis_code)

  mismatches <- df %>%
    left_join(ref, by = c("SpeciesName" = "landis_code")) %>%
    filter(!is.na(fia_spcd), FIA_SPCD != fia_spcd) %>%
    mutate(warning_msg = paste0(SpeciesName, ": lookup has SPCD ", FIA_SPCD,
                                " but reference table expects ", fia_spcd))

  if (nrow(mismatches) > 0) {
    warning(
      "Possible FIA SPCD mismatches in lookup CSV (check species assignments):\n",
      paste("  •", mismatches$warning_msg, collapse = "\n"), "\n",
      "FIA SPCD reference: https://apps.fs.usda.gov/fia/datamart/CSV/REF_SPECIES.csv"
    )
  }

  df
}


# =============================================================================
# load_fia_trees
# Load FIA TREE table for the requested states.
#
# Strategy:
#   1. For each state, check {STATE}_TREE.csv in fia_dir (non-empty only).
#   2. Download missing/empty states via rFIA::getFIA() into fia_dir.
#   3. Read CSVs directly with read_csv() — avoids rFIA::readFIA() which
#      crashes with "OneIndex" on empty files.
#   4. Selects SPCD, DIA, TOTAGE, BHAGE, STATUSCD; converts units.
#
# Returns filtered tibble (live trees only) with columns:
#   SPCD, dia_cm, TOTAGE (may be NA), BHAGE (may be NA)
# =============================================================================
load_fia_trees <- function(states, fia_dir = NULL) {

  if (!requireNamespace("rFIA", quietly = TRUE)) {
    stop("rFIA is required. Install with: install.packages('rFIA')")
  }

  # Resolve working directory
  work_dir <- if (!is.null(fia_dir) && nchar(trimws(fia_dir)) > 0) {
    dir.create(fia_dir, recursive = TRUE, showWarnings = FALSE)
    fia_dir
  } else {
    d <- file.path(tempdir(), "fia_download")
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    d
  }

  # Determine which states need downloading (missing or 0-byte files)
  tree_paths    <- setNames(file.path(work_dir, paste0(states, "_TREE.csv")),
                            states)
  needs_dl      <- states[
    !file.exists(tree_paths) | file.info(tree_paths[states])$size == 0
  ]

  if (length(needs_dl) > 0) {
    message("rFIA::getFIA() downloading TREE table for: ",
            paste(needs_dl, collapse = ", "), " → ", work_dir)
    tryCatch(
      rFIA::getFIA(states = needs_dl, tables = "TREE",
                   dir = work_dir, nCores = 1),
      error = function(e) stop("rFIA::getFIA() failed: ", conditionMessage(e))
    )
  }

  # Re-check after download
  still_bad <- states[
    !file.exists(tree_paths) | file.info(tree_paths[states])$size == 0
  ]
  if (length(still_bad) > 0) {
    stop("FIA TREE CSV still empty/missing after download for: ",
         paste(still_bad, collapse = ", "))
  }

  # Read directly — more robust than readFIA()
  message("Reading FIA TREE CSVs for: ", paste(states, collapse = ", "))
  raw <- map_dfr(tree_paths, \(p)
    read_csv(p, show_col_types = FALSE,
             col_types = cols(.default = col_guess()))
  )
  names(raw) <- toupper(names(raw))

  # Validate
  req <- c("SPCD", "DIA", "STATUSCD")
  miss <- setdiff(req, names(raw))
  if (length(miss) > 0) stop("FIA TREE CSV missing columns: ",
                              paste(miss, collapse = ", "))

  # Filter to live trees; keep age fields if present
  trees <- raw %>%
    filter(STATUSCD == 1L, !is.na(DIA), DIA > 0) %>%
    mutate(
      SPCD   = as.integer(SPCD),
      dia_cm = DIA * 2.54,
      # TOTAGE / BHAGE: 0 means not measured — treat as NA
      TOTAGE = if ("TOTAGE" %in% names(.))
                 if_else(TOTAGE > 0, as.numeric(TOTAGE), NA_real_) else NA_real_,
      BHAGE  = if ("BHAGE" %in% names(.))
                 if_else(BHAGE > 0, as.numeric(BHAGE), NA_real_) else NA_real_
    ) %>%
    # Combine Populus balsamifera (741) + black cottonwood (747) — WA uses 747
    mutate(SPCD = if_else(SPCD == 747L, 741L, SPCD)) %>%
    select(SPCD, dia_cm, TOTAGE, BHAGE)

  message(sprintf("Loaded %s live trees across %d state(s).",
                  format(nrow(trees), big.mark = ","), length(states)))
  trees
}


# =============================================================================
# estimate_bark_params_from_fia
# Estimate MaximumBarkThickness and AgeDBH per species from FIA TREE data.
#
# MaximumBarkThickness (cm, double-sided):
#   = 2 × bark_ratio × P95_DBH_cm
#   bark_ratio from built-in table (literature coefficients)
#
# AgeDBH (years):
#   Half-saturation constant from Michaelis-Menten DBH~age growth curve:
#   DBH_cm = Dmax × Age / (Age + AgeDBH)
#   Fitted with nls(); Age = TOTAGE if available, else BHAGE + years-to-BH.
#   Falls back to bark-class heuristic when < 15 aged records.
#
# lookup  : tibble from load_species_lookup() — SpeciesName / FIA_SPCD / SCIENTIFIC_NAME
# Returns : tibble with SpeciesCode, AgeDBH, MaximumBarkThickness,
#           ScientificName, FIA_SPCD, Note
# =============================================================================
estimate_bark_params_from_fia <- function(species_codes, fia_trees, lookup) {

  # Map LANDIS codes → FIA SPCD via user lookup
  mapped <- tibble(SpeciesCode = species_codes) %>%
    left_join(lookup, by = c("SpeciesCode" = "SpeciesName"))

  spcd_vec  <- unique(na.omit(mapped$FIA_SPCD))
  trees_sub <- fia_trees %>% filter(SPCD %in% spcd_vec)

  # ---------------------------------------------------------------------------
  # MaximumBarkThickness via bark ratio × P95 DBH
  # ---------------------------------------------------------------------------
  dbh_summary <- trees_sub %>%
    group_by(SPCD) %>%
    summarise(
      n_trees    = n(),
      p95_dia_cm = quantile(dia_cm, 0.95, na.rm = TRUE),
      max_dia_cm = max(dia_cm, na.rm = TRUE),
      .groups    = "drop"
    )

  bark_params <- dbh_summary %>%
    left_join(
      bark_ratio_table %>%
        group_by(fia_spcd) %>%
        slice(1) %>%              # one row per SPCD (741/747 share ratio)
        select(fia_spcd, bark_ratio, source),
      by = c("SPCD" = "fia_spcd")
    ) %>%
    mutate(
      # Double-sided (× 2) bark thickness at 95th-pct DBH
      MaximumBarkThickness = round(2 * bark_ratio * p95_dia_cm, 2),
      bt_source = if_else(!is.na(bark_ratio),
                          paste0("bark_ratio=", bark_ratio, " [", source, "]"),
                          "No bark ratio — default used")
    )

  # ---------------------------------------------------------------------------
  # AgeDBH via MM DBH~Age curve
  # ---------------------------------------------------------------------------
  # Build age column: prefer TOTAGE, else BHAGE + correction
  trees_age <- trees_sub %>%
    left_join(years_to_bh, by = "SPCD") %>%
    mutate(
      ytbh = coalesce(ytbh, 4L),
      Age  = case_when(
        !is.na(TOTAGE)                   ~ TOTAGE,
        !is.na(BHAGE)                    ~ BHAGE + ytbh,
        TRUE                             ~ NA_real_
      )
    ) %>%
    filter(!is.na(Age), Age > 0, dia_cm > 0)

  age_fits <- trees_age %>%
    group_by(SPCD) %>%
    filter(n() >= 15) %>%           # minimum records for a reliable fit
    summarise(
      n_aged = n(),
      age_dbh_fit = {
        dmax_est <- quantile(dia_cm, 0.95, na.rm = TRUE)
        fit <- tryCatch(
          nls(
            dia_cm ~ dmax_est * Age / (Age + theta),
            data    = cur_data(),
            start   = list(theta = 60),
            control = nls.control(maxiter = 300, warnOnly = TRUE)
          ),
          error = function(e) NULL
        )
        if (!is.null(fit)) max(1, round(coef(fit)[["theta"]])) else NA_real_
      },
      age_data_source = if_else(
        any(!is.na(TOTAGE), na.rm = TRUE), "TOTAGE", "BHAGE+correction"
      ),
      .groups = "drop"
    )

  # ---------------------------------------------------------------------------
  # Combine; apply heuristic fallback for AgeDBH where nls failed
  # ---------------------------------------------------------------------------
  results <- bark_params %>%
    left_join(age_fits %>% select(SPCD, n_aged, age_dbh_fit, age_data_source),
              by = "SPCD") %>%
    mutate(
      AgeDBH = case_when(
        !is.na(age_dbh_fit) & age_dbh_fit > 0 ~ as.integer(round(age_dbh_fit)),
        # Bark-class heuristic (thicker bark ≈ slower diameter growth)
        !is.na(MaximumBarkThickness) & MaximumBarkThickness >= 25 ~ 120L,
        !is.na(MaximumBarkThickness) & MaximumBarkThickness >= 12 ~  80L,
        TRUE                                                       ~  50L
      ),
      age_method = case_when(
        !is.na(age_dbh_fit) & age_dbh_fit > 0 ~
          paste0("MM DBH~Age nls (n=", n_aged, ", ", age_data_source, ")"),
        TRUE ~ "Bark-class heuristic (< 15 aged trees)"
      )
    )

  # ---------------------------------------------------------------------------
  # Build output tibble
  # ---------------------------------------------------------------------------
  out <- mapped %>%
    left_join(
      results %>% select(SPCD, AgeDBH, MaximumBarkThickness, bt_source, age_method, n_trees),
      by = c("FIA_SPCD" = "SPCD")
    ) %>%
    mutate(
      MaximumBarkThickness = if_else(!is.na(MaximumBarkThickness),
                                     MaximumBarkThickness, 10.0),
      AgeDBH = if_else(!is.na(AgeDBH), AgeDBH, 100L),
      Note = case_when(
        is.na(FIA_SPCD)  ~ "Not in lookup — defaults used",
        is.na(n_trees)   ~ paste0("SPCD ", FIA_SPCD,
                                   " — no FIA records in selected states; defaults used"),
        TRUE             ~ paste0("n=", n_trees, " | MaxBT: ", bt_source,
                                   " | AgeDBH: ", age_method)
      ),
      ScientificName = coalesce(SCIENTIFIC_NAME, "")
    ) %>%
    select(SpeciesCode, AgeDBH, MaximumBarkThickness, ScientificName, FIA_SPCD, Note)

  out
}


# =============================================================================
# build_default_species_table
# Stub table (AgeDBH=100, MaxBT=10) — fills in before FIA data are fetched.
# =============================================================================
build_default_species_table <- function(species_codes, lookup = NULL) {
  tbl <- tibble(SpeciesCode = species_codes,
                AgeDBH = 100L, MaximumBarkThickness = 10.0)

  if (!is.null(lookup) && nrow(lookup) > 0) {
    tbl <- tbl %>%
      left_join(lookup %>% select(SpeciesName, FIA_SPCD, SCIENTIFIC_NAME),
                by = c("SpeciesCode" = "SpeciesName")) %>%
      mutate(
        ScientificName = coalesce(SCIENTIFIC_NAME, ""),
        Note = if_else(!is.na(FIA_SPCD),
                       paste0("Default — SPCD ", FIA_SPCD,
                              "; click 'Fetch FIA Parameters' to estimate"),
                       "Not in lookup CSV — enter values manually")
      ) %>%
      select(SpeciesCode, AgeDBH, MaximumBarkThickness, ScientificName, FIA_SPCD, Note)
  } else {
    tbl %>%
      mutate(ScientificName = "", FIA_SPCD = NA_integer_,
             Note = "Default — load lookup CSV and fetch FIA data")
  }
}


# =============================================================================
# bark_curve_data
# SCF Eq. 10 curve data for preview plot.
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
