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
#   BARK_THICK : bark thickness in tenths of an inch, one-side measurement
#   TOTAGE     : total tree age (rings at breast height + growth years to BH)
#   STATUSCD   : 1 = live, 2 = dead (filter to 1)
# =============================================================================

pacman::p_load(dplyr, tibble, readr, stringr, purrr)

# rFIA is available on CRAN; install if missing
if (!requireNamespace("rFIA", quietly = TRUE)) {
  message("Installing rFIA from CRAN...")
  install.packages("rFIA")
}

# -----------------------------------------------------------------------------
# LANDIS species code → FIA SPCD lookup table
# Covers full PascalCase binomials (NECN style) and 8-char abbreviations.
# Add rows here for project-specific codes not yet listed.
# -----------------------------------------------------------------------------
landis_to_fia_lookup <- tribble(
  ~landis_code,               ~fia_spcd,  ~common_name,
  # ---- Firs (Abies) -------------------------------------------------------
  "AbiesAmabilis",             11L,        "Pacific silver fir",
  "AbieAmab",                  11L,        "Pacific silver fir",
  "AbiesConcolor",             15L,        "white fir",
  "AbieConc",                  15L,        "white fir",
  "AbiesGrandis",              17L,        "grand fir",
  "AbieGran",                  17L,        "grand fir",
  "AbiesLasiocarpa",           19L,        "subalpine fir",
  "AbieLasi",                  19L,        "subalpine fir",
  "AbiesMagnifica",            20L,        "California red fir",
  "AbieMagn",                  20L,        "California red fir",
  "AbiesProcera",              22L,        "noble fir",
  "AbieProc",                  22L,        "noble fir",
  # ---- Incense-cedar -------------------------------------------------------
  "CalocedrusDecurrens",       81L,        "incense-cedar",
  "CaloDecu",                  81L,        "incense-cedar",
  # ---- Yellow-cedar --------------------------------------------------------
  "ChamaecyparisNootkatensis", 42L,        "Alaska yellow-cedar",
  "ChamNoot",                  42L,        "Alaska yellow-cedar",
  # ---- Larch (Larix) -------------------------------------------------------
  "LarixLyallii",              72L,        "subalpine larch",
  "LariLyal",                  72L,        "subalpine larch",
  "LarixOccidentalis",         73L,        "western larch",
  "LariOcci",                  73L,        "western larch",
  # ---- Maples, Alders, Madrone, Birch (hardwoods) -------------------------
  "AcerMacrophyllum",          312L,       "bigleaf maple",
  "AcerMacr",                  312L,       "bigleaf maple",
  "AlnusRubra",                351L,       "red alder",
  "AlnuRubr",                  351L,       "red alder",
  "ArbutusMenziesii",          361L,       "Pacific madrone",
  "ArbuMenz",                  361L,       "Pacific madrone",
  "BetulaPapyrifera",          375L,       "paper birch",
  "BetuPapi",                  375L,       "paper birch",
  # ---- Spruce (Picea) ------------------------------------------------------
  "PiceaEngelmannii",          93L,        "Engelmann spruce",
  "PiceEnge",                  93L,        "Engelmann spruce",
  "PiceaSitchensis",           98L,        "Sitka spruce",
  "PiceSitc",                  98L,        "Sitka spruce",
  # ---- Pine (Pinus) --------------------------------------------------------
  "PinusAlbicaulis",           101L,       "whitebark pine",
  "PinuAlbi",                  101L,       "whitebark pine",
  "PinusContorta",             108L,       "lodgepole pine",
  "PinuCont",                  108L,       "lodgepole pine",
  "PinusJeffreyi",             116L,       "Jeffrey pine",
  "PinuJeff",                  116L,       "Jeffrey pine",
  "PinusLambertiana",          117L,       "sugar pine",
  "PinuLamb",                  117L,       "sugar pine",
  "PinusMonticola",            119L,       "western white pine",
  "PinuMont",                  119L,       "western white pine",
  "PinusPonderosa",            122L,       "ponderosa pine",
  "PinuPond",                  122L,       "ponderosa pine",
  # ---- Poplar / Aspen (Populus) --------------------------------------------
  "PopulusBalsamifera",        743L,       "balsam poplar",
  "PopuBals",                  743L,       "balsam poplar",
  "PopulusTremuloides",        746L,       "quaking aspen",
  "PopuTrem",                  746L,       "quaking aspen",
  # ---- Douglas-fir ---------------------------------------------------------
  "PseudotsugaMenziesii",      202L,       "Douglas-fir",
  "PseuMenz",                  202L,       "Douglas-fir",
  # ---- Oak -----------------------------------------------------------------
  "QuercusGarryana",           815L,       "Oregon white oak",
  "QuerGarr",                  815L,       "Oregon white oak",
  # ---- Yew, Cedars, Hemlocks -----------------------------------------------
  "TaxusBrevifolia",           231L,       "Pacific yew",
  "TaxuBrev",                  231L,       "Pacific yew",
  "ThujaPlicata",              242L,       "western redcedar",
  "ThujPlic",                  242L,       "western redcedar",
  "TsugaHeterophylla",         263L,       "western hemlock",
  "TsugHete",                  263L,       "western hemlock",
  "TsugaMertensiana",          264L,       "mountain hemlock",
  "TsugMert",                  264L,       "mountain hemlock"
)


# -----------------------------------------------------------------------------
# parse_necn_species
# Read SpeciesCode column from a NECN species parameter CSV.
# Returns a character vector of species codes.
# -----------------------------------------------------------------------------
parse_necn_species <- function(csv_path) {
  if (!file.exists(csv_path)) stop("NECN CSV not found: ", csv_path)
  df <- read_csv(csv_path, show_col_types = FALSE)
  if (!"SpeciesCode" %in% names(df)) {
    stop("No 'SpeciesCode' column found in: ", basename(csv_path))
  }
  df$SpeciesCode
}


# -----------------------------------------------------------------------------
# load_fia_trees
# Download or read FIA TREE table for the given states.
# Returns a data frame with columns: SPCD, DIA, BARK_THICK, TOTAGE.
#
# fia_dir  : path to a directory containing pre-downloaded FIA state CSVs
#            (e.g. OR_TREE.csv).  If NULL or empty, data are downloaded fresh
#            using rFIA::getFIA() into a session-temp folder.
# states   : two-letter state abbreviation vector (default PNW).
# -----------------------------------------------------------------------------
load_fia_trees <- function(states = c("OR", "WA", "ID", "MT"),
                           fia_dir = NULL) {

  if (!requireNamespace("rFIA", quietly = TRUE)) {
    stop("The 'rFIA' package is required. Install it with:\n",
         "  install.packages('rFIA')")
  }

  # Decide whether to read local or download
  if (!is.null(fia_dir) && nchar(fia_dir) > 0 && dir.exists(fia_dir)) {
    message("Reading FIA data from: ", fia_dir)
    db <- rFIA::readFIA(fia_dir, states = states,
                        tables  = "TREE",
                        nCores  = 1,
                        inMemory = TRUE)
  } else {
    dl_dir <- file.path(tempdir(), "fia_download")
    dir.create(dl_dir, showWarnings = FALSE, recursive = TRUE)
    message("Downloading FIA TREE table for: ", paste(states, collapse = ", "))
    db <- rFIA::getFIA(states  = states,
                       tables  = "TREE",
                       dir     = dl_dir,
                       nCores  = 1)
  }

  # Extract TREE table (rFIA stores it as db$TREE)
  trees <- db$TREE
  if (is.null(trees) || nrow(trees) == 0) {
    stop("No TREE records returned from FIA for states: ", paste(states, collapse = ", "))
  }

  trees <- trees %>%
    filter(
      STATUSCD == 1L,         # live trees only
      !is.na(DIA),
      !is.na(BARK_THICK),
      DIA > 0,
      BARK_THICK > 0
    ) %>%
    select(SPCD, DIA, BARK_THICK, any_of("TOTAGE")) %>%
    mutate(
      # Unit conversions
      # DIA        : inches → cm
      # BARK_THICK : tenths-of-an-inch (one side) → cm
      dia_cm       = DIA * 2.54,
      bark_cm      = (BARK_THICK / 10) * 2.54
    )

  trees
}


# -----------------------------------------------------------------------------
# estimate_bark_params_from_fia
# Given a data frame of FIA TREE records (from load_fia_trees) and a vector of
# LANDIS species codes, estimate MaximumBarkThickness and AgeDBH per species.
#
# Returns a tibble with columns:
#   SpeciesCode, AgeDBH, MaximumBarkThickness, CommonName, FIA_SPCD, Note
# -----------------------------------------------------------------------------
estimate_bark_params_from_fia <- function(species_codes,
                                          fia_trees,
                                          lookup = landis_to_fia_lookup) {
  # Map LANDIS codes → FIA SPCD
  mapped <- tibble(SpeciesCode = species_codes) %>%
    left_join(lookup %>% select(landis_code, fia_spcd, common_name),
              by = c("SpeciesCode" = "landis_code"))

  # Unique SPCDs to process
  spcd_vec <- unique(na.omit(mapped$fia_spcd))

  # Filter FIA trees to relevant species
  trees_sub <- fia_trees %>% filter(SPCD %in% spcd_vec)

  # --- Per-species summaries ------------------------------------------------
  bark_summary <- trees_sub %>%
    group_by(SPCD) %>%
    summarise(
      n_trees           = n(),
      # MaximumBarkThickness: 95th percentile of bark thickness (in cm)
      # This approximates the asymptotic bark thickness at very large DBH
      max_bark_cm       = quantile(bark_cm, 0.95, na.rm = TRUE),
      # Bark thickness at top-5% trees by DBH (more direct measurement)
      bark_at_large_dia = {
        q95_dia <- quantile(dia_cm, 0.95, na.rm = TRUE)
        mean(bark_cm[dia_cm >= q95_dia], na.rm = TRUE)
      },
      max_dia_cm        = max(dia_cm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # Use the larger of the two estimates as MaximumBarkThickness
      MaximumBarkThickness = round(pmax(max_bark_cm, bark_at_large_dia, na.rm = TRUE), 2)
    )

  # --- AgeDBH fitting -------------------------------------------------------
  # Only attempt if TOTAGE column is present
  age_fits <- if ("TOTAGE" %in% names(fia_trees)) {
    trees_sub %>%
      filter(!is.na(TOTAGE), TOTAGE > 5) %>%
      group_by(SPCD) %>%
      filter(n() >= 15) %>%                    # need at least 15 aged trees
      summarise(
        n_aged    = n(),
        age_dbh   = {
          this_max_bt <- quantile(bark_cm, 0.95, na.rm = TRUE)
          fit <- tryCatch(
            nls(
              bark_cm ~ this_max_bt * TOTAGE / (TOTAGE + theta),
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
    tibble(SPCD = integer(), n_aged = integer(), age_dbh = numeric())
  }

  # --- Combine and fall back to defaults ------------------------------------
  results <- bark_summary %>%
    left_join(age_fits %>% select(SPCD, n_aged, age_dbh),
              by = "SPCD") %>%
    mutate(
      AgeDBH = case_when(
        !is.na(age_dbh) & age_dbh > 0 ~ as.integer(round(age_dbh)),
        # Rough proxy: larger-barked = slower growth = higher AgeDBH
        MaximumBarkThickness >= 30     ~ 80L,
        MaximumBarkThickness >= 15     ~ 60L,
        TRUE                           ~ 40L
      ),
      fia_note = paste0(
        "FIA SPCD ", SPCD, " | n=", n_trees, " trees",
        if_else(!is.na(n_aged), paste0(" | AgeDBH from nls (n=", n_aged, ")"),
                " | AgeDBH estimated from bark thickness class")
      )
    )

  # --- Build output tibble --------------------------------------------------
  out <- mapped %>%
    left_join(results %>% select(SPCD, AgeDBH, MaximumBarkThickness, fia_note),
              by = c("fia_spcd" = "SPCD")) %>%
    mutate(
      MaximumBarkThickness = if_else(!is.na(MaximumBarkThickness),
                                     MaximumBarkThickness, 10.0),
      AgeDBH               = if_else(!is.na(AgeDBH), AgeDBH, 100L),
      Note = case_when(
        is.na(fia_spcd)              ~ "No FIA mapping — default values (AgeDBH=100, MaxBT=10 cm)",
        is.na(fia_note)              ~ paste0("FIA SPCD ", fia_spcd,
                                               " — no records in selected states; defaults used"),
        TRUE                         ~ fia_note
      ),
      CommonName = coalesce(common_name, "Unknown")
    ) %>%
    select(SpeciesCode, AgeDBH, MaximumBarkThickness, CommonName, FIA_SPCD = fia_spcd, Note)

  out
}


# -----------------------------------------------------------------------------
# build_default_species_table
# Quick stub table for species codes with no FIA data yet — used to populate
# the editable DT so the user can enter values manually.
# -----------------------------------------------------------------------------
build_default_species_table <- function(species_codes,
                                        lookup = landis_to_fia_lookup) {
  tibble(SpeciesCode = species_codes) %>%
    left_join(lookup %>% distinct(landis_code, .keep_all = TRUE) %>%
                select(landis_code, fia_spcd, common_name),
              by = c("SpeciesCode" = "landis_code")) %>%
    mutate(
      AgeDBH               = 100L,
      MaximumBarkThickness = 10.0,
      CommonName           = coalesce(common_name, "Unknown"),
      FIA_SPCD             = fia_spcd,
      Note                 = "Default — click 'Fetch FIA Parameters' to estimate from data"
    ) %>%
    select(SpeciesCode, AgeDBH, MaximumBarkThickness, CommonName, FIA_SPCD, Note)
}


# -----------------------------------------------------------------------------
# Bark thickness curve helper for the preview plot
# Returns a data frame of Age vs BarkThickness for each species in spp_table
# -----------------------------------------------------------------------------
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
