# SCF Fire Calibration Shiny App

Shiny application for calibrating ignition and spread parameters for the
LANDIS-II Social-Climate-Fire (SCF) v4 extension.

---

## App Structure

```
scf_calibration_app/
├── app.R                       # Main Shiny app (UI + server)
└── modules/
    ├── era_helpers.R           # ERA-Interim FWI + ERA5 wind extraction
    ├── ignition_helpers.R      # FPA FOD loading, model fitting, surfaces
    ├── spread_helpers.R        # GeoMAC loading, spread fitting, scoring
    └── ui_helpers.R            # Minor UI utilities
```

---

## Required R Packages

Install all dependencies before running:

```r
install.packages(c(
  "shiny",
  "bslib",
  "terra",
  "dplyr",
  "lubridate",
  "readr",
  "tidyr",
  "ggplot2",
  "glue",
  "DT",
  "leaflet",
  "purrr",
  "stringr",
  "broom",
  "sf",         # needed for leaflet polygon display
  "raster",     # needed for leaflet raster display
  "scales",     # needed for axis formatting
  "pscl"        # needed for Zero-Inflated Poisson ignition model
))
```

---

## Data Downloads

### 1. ERA-Interim FWI Indices
Source: Vitolo et al. (2019) — https://doi.org/10.5281/zenodo.1065400

Download the following NetCDF files (~7 GB each):
- `fwi.nc`  — Fire Weather Index
- `ffmc.nc` — Fine Fuel Moisture Code
- `dmc.nc`  — Duff Moisture Code
- `dc.nc`   — Drought Code

These are global daily reanalysis data at ~0.7° (~80 km) resolution, 1980–2018.
Longitudes are stored 0→360 and are handled transparently by terra.

### 2. ERA5 Wind Components (u10 / v10)
Source: Copernicus Climate Data Store
URL: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels

Variables to download:
- `10m_u_component_of_wind`
- `10m_v_component_of_wind`

Settings:
- Product type: Reanalysis
- Years: 1980–2018
- Time: 12:00 UTC (noon, to match ERA-Interim FWI timing)
- Format: NetCDF

You can download as two separate files (u10.nc, v10.nc) or as a single
combined file containing both variables. The app handles both.

Note: ERA5 (~31 km) and ERA-Interim (~80 km) are used together here.
Both datasets are derived from ECMWF reanalyses representing noon
conditions, making them compatible for this calibration workflow.

### 3. FPA FOD Fire Occurrence Database
Source: Karen Short / USDA Forest Service
URL: https://www.fs.usda.gov/rds/archive/Catalog/RDS-2013-0009.6
File: `FPA_FOD_20221014.gdb` (~500 MB)

Used for ignition calibration. Contains all reported US wildfires 1992–2020
with ignition date, location, cause, and final fire size.

### 4. GeoMAC Historic Perimeters
Source: NIFC / USGS GeoMAC
File: `Historic_Geomac_Perimeters_All_Years_2000_2018.gdb`

Used for spread calibration. Contains daily fire perimeter progressions
for major US wildfires 2000–2018.
Download from: https://www.geomac.gov or NIFC data portal.

---

## Running the App

```r
shiny::runApp("/path/to/scf_calibration_app")
```

Or from RStudio: open `app.R` and click "Run App".

---

## Workflow

### Ignition Tab

1. Enter paths for park shapefile and template raster in the sidebar
2. Load ERA-Interim FWI files
3. Load ERA5 wind files (used in spread tab, but loaded once in sidebar)
4. In the Ignition tab: enter FPA FOD GDB path, choose distribution,
   set kernel parameters, click "Run Ignition Calibration"
5. Review diagnostics, download coefficient CSV and ignition surface TIFs
6. Click "Compute Startup Values" in the ERA Startup Codes sub-tab
   to get recommended FFMC/DMC/DC/FWI startup values

### Spread Tab

1. Enter GeoMAC GDB path, adjust cleaning thresholds
2. Click "Run Spread Fitting"
   - Stage A: Fits MaximumDailySpreadArea parameters (B0–B2) using linear
     regression on daily area growth vs FWI + wind
   - Stage B: Fits SpreadProbability parameters (B0–B3) using logistic
     regression on rasterized perimeter progression
     NOTE: This rasterization loop can be slow for parks with many fires.
     A progress bar tracks pair-by-pair progress.
3. Download coefficient tables as CSV
4. (Optional) In the Candidate Grid tab: define a parameter search grid,
   export SCF parameter snippets, run LANDIS externally for each candidate,
   then return to Score Candidates to rank them.

---

## SCF Parameter File

The coefficient tables export directly as CSVs. To enter them into your
SCF parameter text file, use these mappings:

### Ignition (Eq. 1–3)
```
AccidentalIgnitionB0   <b0 from Accidental row>
AccidentalIgnitionB1   <b1 from Accidental row>
AccidentalIgnitionBZ0  <bz0 from Accidental row>  # ZIP only
AccidentalIgnitionBZ1  <bz1 from Accidental row>  # ZIP only

LightningIgnitionB0    <b0 from Lightning row>
LightningIgnitionB1    <b1 from Lightning row>
LightningIgnitionBZ0   <bz0 from Lightning row>   # ZIP only
LightningIgnitionBZ1   <bz1 from Lightning row>   # ZIP only
```

### Spread Probability (Eq. 5)
```
SpreadProbabilityB0    <B0_Intercept>
SpreadProbabilityB1    <B1_FWI>
SpreadProbabilityB2    <B2_FineFuels>    # placeholder = 0 until fuels raster added
SpreadProbabilityB3    <B3_EffectiveWind>
```

### Maximum Daily Spread Area (Eq. 6)
```
MaximumSpreadAreaB0    <B0_Intercept>
MaximumSpreadAreaB1    <B1_FWI>
MaximumSpreadAreaB2    <B2_EffectiveWind>
```

NOTE: The MaximumDailySpreadArea coefficients are fit on log1p(ha/day) scale.
Back-transform the intercept with expm1() if your SCF version expects raw ha.
Verify against typical observed daily spread rates for your landscape before use.

---

## Known Limitations

- **Fine fuels:** The spread probability B2 coefficient (fine fuels) uses a
  placeholder value of 1.0 everywhere. Add a fine fuels raster path when
  available (LANDFIRE FBFMxx or NECN litter output).
- **Effective wind speed:** The SCF formula incorporates slope angle and aspect
  relative to wind direction (Nelson 2002). The current implementation uses
  plain wind speed as a proxy. Add slope/aspect rasters and wire in the full
  formula to improve the B3 estimate.
- **ERA-Interim resolution:** ~80 km spatial resolution. Park-mean FWI smooths
  over internal climate gradients. Acceptable for most NPS landscapes; flag
  for very large or topographically complex parks.
- **Candidate scoring:** Requires external LANDIS runs. The app generates the
  parameter snippets and ingests FireEventLog outputs; the actual LANDIS
  simulations must be run outside the app.
