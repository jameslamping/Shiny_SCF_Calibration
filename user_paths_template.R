# =============================================================================
# user_paths_template.R
#
# HOW TO USE
#   1. Copy this file to `user_paths.R` in the same folder.
#   2. Fill in the paths you have; leave the rest as "".
#   3. Restart the app.  Every path box starts with your value.
#
# `user_paths.R` is gitignored: your paths stay on your machine and are never
# committed.  The diagnostic and analysis scripts read the same file.
#
# Use FORWARD slashes, also on Windows:  "C:/Users/me/data/fwi.nc"
# (a single backslash in an R string is an escape character and will fail).
#
# On start-up the R console lists any entry the app doesn't recognise and any
# path that doesn't exist on this machine.
# =============================================================================

user_paths <- list(

  # ---- Sidebar: landscape setup ---------------------------------------------
  shp_path            = "",   # calibration boundary (.shp / .gpkg)
  tif_path            = "",   # LANDIS template raster (e.g. ecoregions .tif)
  era_fwi_path        = "",   # ERA-Interim fwi.nc   (only for the ERA FWI source)
  era_ffmc_path       = "",   # ERA-Interim ffmc.nc
  era_dmc_path        = "",   # ERA-Interim dmc.nc
  era_dc_path         = "",   # ERA-Interim dc.nc
  gridmet_cache_dir   = "",   # GRIDMET download cache (blank = "gridmet_cache" in the app folder)
  landis_climate_path = "",   # LANDIS daily climate CSV (runtime FWI check)

  # ---- Ignition calibration -------------------------------------------------
  fpa_gdb_path        = "",   # FPA_FOD_20221014.gdb
  events_log_path     = "",   # socialclimatefire-events-log.csv (Ignition > LANDIS Validation)

  # ---- Species table --------------------------------------------------------
  lookup_csv_path     = "",   # LANDIS <-> FIA species lookup CSV
  necn_csv_path       = "",   # NECN species parameter CSV (to filter species)
  fia_local_dir       = "",   # folder of {STATE}_TREE.csv files

  # ---- Spread calibration ---------------------------------------------------
  nirops_shp_path        = "", # NIROPS perimeters shapefile
  geomac_gdb_path        = "", # Historic_Geomac_Perimeters_All_Years_2000_2018.gdb
  fired_gpkg_path        = "", # FIRED CONUS DAILY GeoPackage
  fired_events_gpkg_path = "", # FIRED CONUS EVENTS GeoPackage
  fine_fuels_path        = "", # local 0-1 fine fuels GeoTIFF (optional)
  scf_runs_root          = "", # folder holding candidate_0001/, candidate_0002/, ...

  # ---- Mortality calibration ------------------------------------------------
  solus_clay_dir      = "",   # SOLUS100 clay folder
  sg_clay_dir         = "",   # SoilGrids clay folder
  mtbs_out_dir        = "",   # where downloaded MTBS rasters are saved
  mtbs_dnbr_dir       = "",   # folder of local MTBS dNBR TIFFs
  tc_cache_dir        = "",   # TerraClimate cache (blank = session temp dir)
  et_local_path       = "",   # local annual PET raster (mm/yr)
  cwd_local_path      = "",   # local annual CWD raster (mm/yr)
  ladder_fuels_path   = "",   # local ladder fuel biomass raster (g m-2)
  ftm_dir             = "",   # Fire and Tree Mortality database, RDS-2020-0001-2/Data

  # ---- Model validation -----------------------------------------------------
  val_events_path     = "",   # socialclimatefire-events-log.csv (or its folder)
  mtbs_local_path     = "",   # MTBS fire occurrence points shapefile
  fri_maps_dir        = "",   # LANDIS run folder with SCF fire maps
  fri_lf_tif_path     = "",   # LANDFIRE mean fire return interval GeoTIFF

  # ---- Session ---------------------------------------------------------------
  rdata_load_path     = ""    # a saved scf_cache_*.RData
)
