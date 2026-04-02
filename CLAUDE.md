# SCF Fire Calibration Shiny App — Project Context

## What this is
A Shiny app for calibrating LANDIS-II Social-Climate-Fire (SCF v4) extension 
parameters. Produces ignition coefficients, ignition allocation maps, spread 
probability coefficients, and maximum daily spread area coefficients.

## App structure
app.R                        # Main Shiny UI + server
modules/era_helpers.R        # ERA-Interim FWI + ERA5 wind extraction
modules/ignition_helpers.R   # FPA FOD loading, ZIP/Poisson fitting, kernel surfaces
modules/spread_helpers.R     # GeoMAC loading, perimeter pairs, spread fitting, scoring
modules/ui_helpers.R         # Minor UI utilities

## Spatial design (critical)
- cal_vect      : Large regional calibration boundary (e.g. WA Cascades)
                  Used to clip ERA data, FPA FOD, and GeoMAC perimeters
- template_r    : Small LANDIS park grid raster
                  Used only for ignition surface generation and perimeter rasterization
- working_crs   : Always taken from template_r; everything reprojects to match it
- All leaflet/tmap display reprojects to EPSG:4326 on the fly

## Key packages
pacman::p_load(shiny, bslib, terra, dplyr, lubridate, readr, tidyr,
               ggplot2, glue, DT, tmap, purrr, stringr, broom, sf, scales)
pscl  -- for Zero-Inflated Poisson ignition model (loaded conditionally)

## Data sources
- ERA-Interim FWI (ffmc/dmc/dc/fwi .nc): Vitolo et al. 2019, Zenodo DOI 10.5281/zenodo.1065400
- ERA5 wind (u10/v10 .nc): Copernicus CDS, daily 12:00 UTC, 1980-2018
- FPA FOD ignitions: Short 2022, FPA_FOD_20221014.gdb
- GeoMAC perimeters: Historic_Geomac_Perimeters_All_Years_2000_2018.gdb

## Ignition tab — COMPLETE AND WORKING
Workflow:
1. load_fpa_fod()              FPA FOD -> filter to cal_vect, project to working_crs
2. build_ignition_model_data() Daily count x FWI join, zero-fill
3. fit_ignition_models()       Poisson or ZIP per ignition type (Lightning/Accidental)
4. make_ignition_surfaces()    terra focal() Gaussian kernel on template grid (NOT a loop)
5. compute_startup_values()    ERA climatology -> FFMC/DMC/DC/FWI startup values

Key fixes applied:
- FPA FOD DISCOVERY_DATE is M/D/YYYY format; parse with "%m/%d/%Y" first
- Always use FIRE_YEAR (integer column) for year filter, not derived from date
- Kernel surface uses terra::rasterize() + focal() for speed, not an R loop
- tmap used for all maps (not leaflet) due to rendering reliability with terra objects
- output$landscape_preview_ui renderUI wraps tmapOutput("landscape_map") — 
  this container must exist before output$landscape_map renderTmap fires

## Spread tab — IN PROGRESS
Current state:
- UI is fully built in app.R (all panels, inputs, download buttons)
- spread_helpers.R has skeleton functions but needs testing/debugging:
  - load_geomac_perimeters()    needs same CRS fix pattern as FPA FOD
  - bind_climate_to_pairs()     builds sequential perimeter pairs + joins ERA climate
  - fit_max_daily_area()        lm: log1p(daily_area_ha) ~ FWI + EffectiveWind
  - fit_spread_probability()    logistic glm from rasterized perimeter pairs (slow loop)
  - score_candidates()          scores post-LANDIS runs against FPA FOD targets
  - write_scf_snippet()         formats parameters for SCF text file

Known issues to address on spread side:
- fit_spread_probability() uses an R loop over perimeter pairs + rasterization;
  may need same optimization treatment as make_kernel_surface() received
- EffectiveWind is currently a placeholder (= WindSpeed_kmh directly);
  full SCF formula incorporates slope/aspect but those rasters are not yet wired in
- Fine fuels predictor is placeholder (= 1.0 everywhere)
- score_candidates() depends on external LANDIS runs being present on disk;
  read_scf_events() parser may need adjustment for actual SCF v4 output format

## SCF parameter mappings
Ignition:
  b0, b1       -> AccidentalIgnitionB0/B1, LightningIgnitionB0/B1
  bz0, bz1     -> AccidentalIgnitionBZ0/BZ1 (ZIP only)

Spread probability (SCF Eq. 5):
  B0_Intercept    -> SpreadProbabilityB0
  B1_FWI          -> SpreadProbabilityB1
  B2_FineFuels    -> SpreadProbabilityB2 (placeholder, = 0 until fuels raster added)
  B3_EffectiveWind -> SpreadProbabilityB3

Max daily spread area (SCF Eq. 6):
  B0_Intercept    -> MaximumSpreadAreaB0
  B1_FWI          -> MaximumSpreadAreaB1
  B2_EffectiveWind -> MaximumSpreadAreaB2
  Note: fit on log1p(ha/day) scale; verify back-transform before entering in SCF

## Style conventions
- terra for all raster/vector ops (not sf/raster directly except for tmap interop)
- dplyr pipelines throughout
- Error messages include extent printouts to help diagnose CRS/spatial mismatches
- pacman::p_load() for all package loading
- Working CRS always from template_r, never hardcoded