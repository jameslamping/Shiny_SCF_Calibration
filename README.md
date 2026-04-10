<p align="center">
  <img src="scf_fire_calibration_icon.svg" width="350" alt="SCF Fire Calibration App icon"/>
</p>

# SCF Fire Calibration Shiny App

An interactive R Shiny application for calibrating the
[LANDIS-II Social-Climate-Fire (SCF) v4](https://github.com/LANDIS-II-Foundation/Extension-Social-Climate-Fire)
extension. The app derives ignition coefficients, spread probability
coefficients, maximum daily spread area coefficients, bark thickness
parameters, and all required spatial maps directly from publicly available
fire and climate data.

---

## Contents

1. [Quick Start](#1-quick-start)
2. [Required Data](#2-required-data)
3. [App Structure](#3-app-structure)
4. [Sidebar — Landscape Setup](#4-sidebar--landscape-setup)
5. [Tab 1 — Ignition Calibration](#5-tab-1--ignition-calibration)
6. [Tab 2 — Species Table](#6-tab-2--species-table)
7. [Tab 3 — Spread Calibration](#7-tab-3--spread-calibration)
8. [SCF Parameter Mappings](#8-scf-parameter-mappings)
9. [Output Files Reference](#9-output-files-reference)
10. [Troubleshooting](#10-troubleshooting)
11. [Data Sources & Citations](#11-data-sources--citations)

---

## 1. Quick Start

### Prerequisites

- R ≥ 4.4.0 and RStudio (recommended)
- The `renv` package for reproducible package management

### Installation

```r
# 1. Clone the repository
#    git clone https://github.com/jameslamping/Shiny_SCF_Calibration

# 2. Open app.R in RStudio — renv activates automatically via .Rprofile
#    OR from the R console in the project directory:

install.packages("renv")          # if not already installed
renv::restore()                   # installs all packages at exact locked versions
                                  # (~5–10 min on first run)

# 3. Launch the app
shiny::runApp()
```

`renv::restore()` reads `renv.lock` and installs every package at the exact
version used during development — including GitHub-sourced packages (`rFIA`,
`climateR`, `ncmeta`). This ensures the app runs identically regardless of
future package updates.

To update the lockfile after adding or upgrading packages:

```r
renv::snapshot()
# then commit renv.lock
```

---

## 2. Required Data

The app relies on several large external datasets that are too big for
version control. Download them once and store locally; the app accepts file
paths entered in the sidebar or tab inputs.

### 2.1 ERA-Interim FWI Indices *(~28 GB total, four files)*

**Source:** Vitolo et al. (2019) — <https://doi.org/10.5281/zenodo.1065400>

Download four NetCDF files:

| File | Variable |
|------|----------|
| `fwi.nc`  | Fire Weather Index |
| `ffmc.nc` | Fine Fuel Moisture Code |
| `dmc.nc`  | Duff Moisture Code |
| `dc.nc`   | Drought Code |

Global daily reanalysis at ~80 km resolution, 1980–2018. Longitudes are
stored 0→360; `terra` handles this transparently.

### 2.2 GRIDMET Wind *(fetched automatically — no download required)*

Daily 10 m wind speed and direction are fetched directly from
[GRIDMET](https://www.climatologylab.org/gridmet.html) via OPeNDAP using
the `climateR` package. CONUS coverage at ~4 km, 1979–present. Simply click
**Fetch GRIDMET Wind** in the sidebar after loading the calibration boundary.

### 2.3 FPA Fire Occurrence Database *(~500 MB)*

**Source:** Short (2022) — <https://www.fs.usda.gov/rds/archive/Catalog/RDS-2013-0009.6>

File: `FPA_FOD_20221014.gdb`

All reported US wildfires 1992–2020 with ignition date, location (latitude/
longitude), cause (lightning vs human), and final fire size. Used for
ignition calibration.

### 2.4 GeoMAC Historic Fire Perimeters

**Source:** NIFC / USGS GeoMAC

File: `Historic_Geomac_Perimeters_All_Years_2000_2018.gdb`

Daily fire perimeter progressions for major US wildfires 2000–2018. Used
for spread calibration. Available from the
[NIFC data portal](https://www.nifc.gov/fire-information/maps-imagery-data-products/historical-fire-data).

### 2.5 LANDFIRE FBFM40 *(optional; fetched via REST API)*

Fine fuel load index (Scott & Burgan 40 fire behavior fuel models) used as
the B2 predictor in the spread probability equation. Select
**Fetch LANDFIRE FBFM40** in the Spread Calibration inputs; the app
downloads the raster for your calibration boundary via the LFPS REST API.
Alternatively, supply any 0–1 scaled GeoTIFF via the **Local file** option.

### 2.6 FIA Tree Data *(for species table only)*

Downloaded automatically by state via the `rFIA` package when you click
**Fetch FIA Bark Parameters** in the Species Table tab. To avoid
re-downloading (~50–500 MB per state), provide the path to a directory of
pre-downloaded FIA state CSV files (e.g. `OR_TREE.csv`). Files from the
[USFS FIA DataMart](https://apps.fs.usda.gov/fia/datamart/) are compatible.

---

## 3. App Structure

```
scf_calibration_app/
├── app.R                        # UI + server (all tabs)
├── renv.lock                    # Pinned package versions
├── .Rprofile                    # Auto-activates renv on project open
├── renv/
│   ├── activate.R               # renv bootstrap (do not edit)
│   └── settings.json
└── modules/
    ├── era_helpers.R            # ERA-Interim FWI extraction; GRIDMET wind fetch
    ├── ignition_helpers.R       # FPA FOD loading, ZIP/Poisson fitting, kernel surfaces
    ├── spread_helpers.R         # GeoMAC loading, perimeter pairs, spread fitting,
    │                            #   effective wind, LANDFIRE REST API, terrain
    ├── species_helpers.R        # FIA bark thickness estimation, species lookup
    └── ui_helpers.R             # Minor UI utilities
```

### Spatial design

The app uses two distinct spatial extents:

| Object | Purpose |
|--------|---------|
| **Calibration boundary** (`cal_vect`) | Large regional polygon (e.g. WA Cascades + Coast Range). Clips ERA-Interim FWI, GRIDMET wind, FPA FOD, and GeoMAC perimeters. Should be large enough to capture the full fire climate driving your landscape. |
| **LANDIS template raster** (`template_r`) | Your park-specific LANDIS grid raster. Defines the working CRS, cell size, and spatial extent for all output maps (ignition surfaces, terrain maps, suppression maps). |

All data are reprojected to the CRS of `template_r` during processing.
Interactive maps display in EPSG:4326 (WGS84).

---

## 4. Sidebar — Landscape Setup

Complete these steps before running any calibration tab.

| Step | Input | Notes |
|------|-------|-------|
| 1 | **Calibration Boundary (.shp)** | Regional boundary shapefile. The `.dbf`, `.shx`, and `.prj` files must be in the same directory. |
| 2 | **LANDIS Template Raster (.tif)** | Any LANDIS input raster (e.g. climate regions, initial communities). Defines cell size and working CRS. |
| 3 | **ERA-Interim FWI Files** | Enter paths to `fwi.nc`, `ffmc.nc`, `dmc.nc`, `dc.nc`. |
| 4 | **GRIDMET Wind** | Click **Fetch GRIDMET Wind** — no local files needed. Requires internet and that the calibration boundary is loaded first. |
| 5 | **Calibration Period** | Slider to set the year range used across all tabs (default 1992–2018). |

After loading the template raster, a landscape preview map appears at the
bottom of the sidebar.

---

## 5. Tab 1 — Ignition Calibration

Calibrates the fire ignition parameters in SCF using historical ignition
records from the FPA FOD database.

### Workflow

1. Enter the path to `FPA_FOD_20221014.gdb` in the Ignition Inputs panel.
2. Select the **distribution** for count modelling:
   - **Zero-Inflated Poisson (ZIP):** Recommended. Handles the large number
     of zero-ignition days common in most landscapes.
   - **Poisson:** Simpler; use if ZIP fails to converge.
3. Set **kernel bandwidth** and **max distance** for the ignition surface:
   - *Bandwidth* (σ, metres): standard deviation of the Gaussian kernel
     applied around each historical ignition point. Typical range 5–15 km.
     Larger values spread ignition probability more broadly.
   - *Max distance* (metres): hard cutoff beyond which kernel weight drops
     to zero. Set to ~3–4× the bandwidth to retain the full Gaussian tail.
4. Click **Run Ignition Calibration**.

### Output tabs

| Sub-tab | Content |
|---------|---------|
| **Overview Map** | Spatial distribution of FPA FOD ignitions coloured by type |
| **Coefficients** | Fitted b0, b1, bz0, bz1 per ignition type; download CSV |
| **Diagnostics** | Observed vs fitted ignition counts over time |
| **ERA Startup Codes** | Click **Compute Startup Values** to derive recommended FFMC/DMC/DC/FWI startup values from ERA climatology over the calibration boundary |
| **Landscape Maps** | Download all spatial outputs (see [Output Files](#9-output-files-reference)) |

### Landscape Maps sub-tab

Three sections under the **Landscape Maps** sub-tab generate all spatial
outputs required by SCF:

**Ignition Allocation Surfaces**
Gaussian kernel-smoothed probability surfaces on the LANDIS template grid.
Download as double-precision float rasters (FLT8S):
- `Lightning_Ignition_Map.tif`
- `Accidental_Ignition_Map.tif`
- `Rx_Ignition_Map.tif`

**Terrain Maps**
Click **Fetch DEM & Generate Terrain Maps** to download a SRTM/3DEP DEM
via `elevatr` (zoom level 8, ~600 m), compute slope and aspect, and
resample to the LANDIS template grid. Preview maps are displayed in the app.
Download as 16-bit signed integer rasters (INT2S):
- `GroundSlope.tif` — slope in degrees
- `UphillSlope.tif` — uphill azimuth 0–359° from North

**Suppression Maps**
Blank suppression maps (value = 0.0 everywhere) for initial calibration.
Replace with zone-based maps once fire behaviour is validated.
Download as single-precision float rasters (FLT4S):
- `Suppression_Accidental_3Zones.tif`
- `Suppression_Lightning_3Zones.tif`
- `Suppression_Rx_3Zones.tif`

---

## 6. Tab 2 — Species Table

Generates `Fire_Spp_Table.csv`, which SCF uses to compute bark thickness
and fire-induced tree mortality via:

**SCF Equation 10:**
```
BarkThickness = (MaximumBarkThickness × Age) / (Age + AgeDBH)
```

where `MaximumBarkThickness` (cm, double-sided) is the asymptotic bark
thickness at maximum observed DBH, and `AgeDBH` (years) is the age at which
bark thickness equals half the maximum.

Parameters are estimated from USFS FIA tree measurement data.

### Workflow

**Step 1 — Load Species Lookup CSV**

Provide the path to a CSV with three columns:

| Column | Description |
|--------|-------------|
| `SpeciesName` | LANDIS species code (must match exactly what is used in your NECN and other extension files) |
| `FIA_SPCD` | Integer FIA species code — verify against [REF_SPECIES.csv](https://apps.fs.usda.gov/fia/datamart/CSV/REF_SPECIES.csv) |
| `SCIENTIFIC_NAME` | Scientific name (display only) |

The app will flag any SPCD values that differ from its internal reference
table. Mismatched SPCDs will pull data for the wrong species from FIA — 
check the warning banner carefully before proceeding.

**Step 2 — Filter or Add Species**

- *All species in lookup:* Uses every row in the lookup CSV.
- *Filter to NECN CSV:* Loads the `SpeciesCode` column from a NECN species
  parameter CSV and restricts the table to those species.
- *Manual entry:* Type species codes one per line.

Click **Build Species Table** to create a stub table with default values
(AgeDBH = 100, MaxBarkThickness = 10 cm).

**Step 3 — Fetch FIA Bark Parameters**

1. Select the US states to query (multi-select; default OR, WA, ID, MT).
2. Optionally enter a **local FIA directory** containing pre-downloaded
   `{STATE}_TREE.csv` files. Any empty or missing files for a state are
   automatically downloaded via `rFIA::getFIA()` and saved to this folder
   for future use. Leave blank to download to a session temp folder.
3. Click **Fetch FIA Bark Parameters**.

For each species the app:
- **MaximumBarkThickness:** Computes `2 × bark_ratio × P95_DBH_cm`, where
  `bark_ratio` (one-side bark thickness / DBH) comes from a built-in
  literature table (Ryan & Reinhardt 1988; Cansler et al. 2020) and
  `P95_DBH_cm` is the 95th-percentile DBH across all live FIA trees of
  that species in the selected states.
- **AgeDBH:** Fits a Michaelis-Menten DBH growth curve
  (`DBH = Dmax × Age / (Age + AgeDBH)`) to FIA records with age
  measurements (`TOTAGE` preferred; `BHAGE` + species-specific
  years-to-breast-height correction as fallback). Falls back to a
  bark-thickness-class heuristic when fewer than 15 aged records exist.

### Output tabs

| Sub-tab | Content |
|---------|---------|
| **Species Parameters** | Editable table — double-click AgeDBH or MaxBarkThickness cells to override values |
| **Bark Thickness Curves** | SCF Eq. 10 preview curves for all species at once |
| **Species Lookup Table** | Contents of the loaded lookup CSV |

Download `Fire_Spp_Table.csv` (3-column SCF format) or the annotated
version with scientific names, FIA SPCDs, and estimation notes.

---

## 7. Tab 3 — Spread Calibration

Calibrates fire spread parameters using historical daily fire perimeter
progressions from GeoMAC.

### Workflow

1. Enter the path to `Historic_Geomac_Perimeters_All_Years_2000_2018.gdb`.
2. Adjust **cleaning thresholds**:

| Parameter | Purpose |
|-----------|---------|
| Max Gap Days — Spread Probability | Only perimeter pairs ≤ N days apart feed the logistic spread probability model. Shorter gaps give a cleaner signal. |
| Max Gap Days — Daily Spread Area | Maximum day gap for the daily area model. Growth is divided by gap length, so wider gaps are tolerable. |
| Negative Growth Tolerance (ha) | Drops pairs where a later perimeter is more than N ha smaller than the earlier (GeoMAC remapping noise vs real data problem). |

3. Set **cell size** for perimeter rasterization (defaults to the LANDIS
   template resolution).
4. Adjust **sampling controls** for the logistic regression:
   - *Failure : Success ratio:* For each newly-burned cell (success), sample
     N adjacent unburned cells (failure). Higher ratios give more balanced
     logistic regression.
   - *Max cells per pair:* Caps cell count per perimeter pair to prevent
     large fires from dominating.
5. Select **fine fuels source** (None / LANDFIRE / local GeoTIFF).
6. Click **Run Spread Fitting**.

The fitting runs in two stages:

**Stage A — Maximum Daily Spread Area**
Linear regression on `log1p(ha/day)` vs FWI and effective wind speed.
Maps to `MaximumSpreadAreaB0`, `B1`, `B2` in the SCF parameter file.

**Stage B — Spread Probability**
Logistic regression on cell-level burn/no-burn outcomes vs FWI, fine fuels
index, and effective wind speed. Maps to `SpreadProbabilityB0`–`B3`.

### Effective Wind Speed

Rather than using raw wind speed, the app computes effective wind speed at
each fire centroid incorporating topographic enhancement (Rothermel 1972):

```
EffectiveWind = max(0, WindSpeed × cos(WindDir − UphillAzimuth))
              + tan(slope_deg × π/180) × 100 × 0.3
```

Slope and aspect are extracted from a DEM fetched via `elevatr` at the
centroid of each perimeter pair. The wind × terrain interaction ensures that
the calibrated B3 coefficient matches the same effective wind metric SCF
computes at runtime.

### Output tabs

| Sub-tab | Content |
|---------|---------|
| **Fit Coefficients** | Coefficient tables for both models; download CSV |
| **Maps** | Fine fuels index map; GeoMAC calibration fires coloured by year |
| **Diagnostics** | Spread probability vs FWI and wind; ROC curve; observed vs predicted daily area scatter |
| **Terrain** | Effective vs raw wind comparison; slope distribution at fire centroids; wind-to-upslope alignment; slope and aspect maps |
| **Perimeter Pairs** | Cleaned perimeter progression table; download CSV |
| **Candidate Grid** | *(Optional)* Define a parameter search grid, export SCF snippet files, run LANDIS externally, return to score candidates |
| **Score Candidates** | Import FireEventLog outputs from external LANDIS runs to score and rank candidate parameter sets |

---

## 8. SCF Parameter Mappings

### Ignition (SCF Eqs. 1–3)

```
AccidentalIgnitionB0    ← b0  (Accidental row in coefficient table)
AccidentalIgnitionB1    ← b1
AccidentalIgnitionBZ0   ← bz0  (ZIP only)
AccidentalIgnitionBZ1   ← bz1  (ZIP only)

LightningIgnitionB0     ← b0  (Lightning row)
LightningIgnitionB1     ← b1
LightningIgnitionBZ0    ← bz0  (ZIP only)
LightningIgnitionBZ1    ← bz1  (ZIP only)
```

### Spread Probability (SCF Eq. 5)

```
SpreadProbabilityB0     ← B0_Intercept
SpreadProbabilityB1     ← B1_FWI
SpreadProbabilityB2     ← B2_FineFuels    (0.0 if fine fuels = None)
SpreadProbabilityB3     ← B3_EffectiveWind
```

### Maximum Daily Spread Area (SCF Eq. 6)

```
MaximumSpreadAreaB0     ← B0_Intercept
MaximumSpreadAreaB1     ← B1_FWI
MaximumSpreadAreaB2     ← B2_EffectiveWind
```

> **Note:** The MaximumSpreadArea model is fit on a `log1p(ha/day)` scale.
> Verify that your SCF version expects coefficients on the log scale (most do).
> Back-transform the intercept with `expm1()` only if SCF expects raw ha.

### Bark Thickness (SCF Eq. 10)

```
Species_CSV_File  ← Fire_Spp_Table.csv   (SpeciesCode, AgeDBH, MaximumBarkThickness)
```

---

## 9. Output Files Reference

| File | Tab | Format | Description |
|------|-----|--------|-------------|
| `ignition_coefficients.csv` | Ignition → Coefficients | CSV | b0, b1, bz0, bz1 per type |
| `startup_values.csv` | Ignition → ERA Startup | CSV | FFMC/DMC/DC/FWI startup values |
| `Lightning_Ignition_Map.tif` | Ignition → Landscape Maps | FLT8S | Lightning ignition probability surface |
| `Accidental_Ignition_Map.tif` | Ignition → Landscape Maps | FLT8S | Accidental ignition probability surface |
| `Rx_Ignition_Map.tif` | Ignition → Landscape Maps | FLT8S | Prescribed fire ignition surface |
| `GroundSlope.tif` | Ignition → Landscape Maps | INT2S | Ground slope in degrees |
| `UphillSlope.tif` | Ignition → Landscape Maps | INT2S | Uphill azimuth 0–359° from North |
| `Suppression_Accidental_3Zones.tif` | Ignition → Landscape Maps | FLT4S | Blank suppression map (0.0) |
| `Suppression_Lightning_3Zones.tif` | Ignition → Landscape Maps | FLT4S | Blank suppression map (0.0) |
| `Suppression_Rx_3Zones.tif` | Ignition → Landscape Maps | FLT4S | Blank suppression map (0.0) |
| `Fire_Spp_Table.csv` | Species Table | CSV | SpeciesCode, AgeDBH, MaximumBarkThickness |
| `Fire_Spp_Table_annotated.csv` | Species Table | CSV | Full table with FIA SPCDs and notes |
| `spread_coefficients.csv` | Spread → Fit Coefficients | CSV | All spread model coefficients |
| `geomac_clean_pairs.csv` | Spread → Perimeter Pairs | CSV | Cleaned perimeter pair dataset |

---

## 10. Troubleshooting

### FIA data fetch fails: "OneIndex" error or empty files

Your `_TREE.csv` files may be 0-byte placeholders (common if the directory
was created by a cloud sync service before files downloaded). The app detects
empty files and re-downloads them via `rFIA::getFIA()` automatically. If you
see this error, check that R has write access to the FIA directory and that
you have an internet connection.

### GRIDMET wind fetch fails

GRIDMET data is fetched via OPeNDAP and requires an active internet
connection. Ensure the calibration boundary is loaded first (the fetch uses
the boundary extent). Retrying usually succeeds if the server is temporarily
unavailable.

### Spread fitting is very slow

The spread probability stage rasterizes each perimeter pair and samples
cells across the calibration boundary. Processing time scales with the number
of fires, the rasterization cell size, and the failure-to-success ratio.
Increasing the cell size, reducing the calibration period, or reducing the
max cells per pair will speed up fitting at some cost to precision.

### SPCD mismatch warnings when loading species lookup CSV

The app compares your FIA SPCD values against an internal reference table
derived from standard FIA codes. If a species code in your lookup CSV maps
to a different integer than expected, the app shows a warning. Incorrect
SPCDs will cause the wrong species to be matched in FIA data. Verify against
[REF_SPECIES.csv](https://apps.fs.usda.gov/fia/datamart/CSV/REF_SPECIES.csv).

### renv::restore() is slow or fails

On first use, renv installs all packages from source or CRAN binaries. This
can take 5–15 minutes. If a GitHub-sourced package fails (climateR, rFIA),
ensure you have a working internet connection and that `pak` or `remotes` is
available. Run `renv::restore(prompt = FALSE)` to skip confirmation prompts.

---

## 11. Data Sources & Citations

**ERA-Interim FWI:**
Vitolo, C., Di Giuseppe, F., Barnard, C., Coughlan, R., San-Miguel-Ayanz, J.,
Libertá, G., & Krzeminski, B. (2019). ERA-interim based global fire weather
index series. *Scientific Data*, 6(1), 190032.
<https://doi.org/10.5281/zenodo.1065400>

**GRIDMET:**
Abatzoglou, J. T. (2013). Development of gridded surface meteorological data
for ecological applications and modelling. *International Journal of
Climatology*, 33(1), 121–131. <https://doi.org/10.1002/joc.3413>

**FPA Fire Occurrence Database:**
Short, K. C. (2022). Spatial wildfire occurrence data for the United States,
1992–2020 (FPA FOD). Forest Service Research Data Archive.
<https://doi.org/10.2737/RDS-2013-0009.6>

**GeoMAC Historic Perimeters:**
National Interagency Fire Center (NIFC). Historic Geomac Perimeters 2000–2018.
<https://www.nifc.gov/fire-information/maps-imagery-data-products/historical-fire-data>

**LANDFIRE FBFM40:**
LANDFIRE Program (2024). LF2024 Fire Behavior Fuel Model 40 (FBFM40), CONUS.
U.S. Department of the Interior. <https://landfire.gov>

**FIA Tree Data:**
U.S. Forest Service Forest Inventory and Analysis Program.
FIADB DataMart. <https://apps.fs.usda.gov/fia/datamart/>

**Bark Thickness Coefficients:**
Ryan, K. C., & Reinhardt, E. D. (1988). Predicting postfire mortality of
seven western conifers. *Canadian Journal of Forest Research*, 18(10),
1291–1297. <https://doi.org/10.1139/x88-199>

Cansler, C. A., Hood, S. M., Varner, J. M., et al. (2020). The Fire and
Tree Mortality Database, for empirical modeling of individual tree mortality
after fire. *Scientific Data*, 7(1), 194.
<https://doi.org/10.1038/s41597-020-0522-7>

**SCF Extension:**
Shinneman, D. J., et al. Social-Climate-Fire (SCF) v4 LANDIS-II Extension.
<https://github.com/LANDIS-II-Foundation/Extension-Social-Climate-Fire>
