# SCF Fire Calibration Methods

## Overview

Fire behavior in the Social-Climate-Fire (SCF v4) extension for LANDIS-II is governed
by three sets of empirically calibrated parameters: (1) ignition probability coefficients
that determine how many fires start each day as a function of fire weather, (2) spread
probability coefficients that determine whether fire advances from one cell to an adjacent
cell, and (3) maximum daily spread area coefficients that set a ceiling on how much area
a single fire can burn in one day. Each parameter set is calibrated independently against
historical fire records and gridded climate reanalysis data.

---

## Data Sources

### Fire Weather Index (FWI)
Daily FWI system components — Fine Fuel Moisture Code (FFMC), Duff Moisture Code (DMC),
Drought Code (DC), and the composite Fire Weather Index (FWI) — were obtained from the
ERA-Interim reanalysis (Vitolo et al. 2019; Zenodo DOI 10.5281/zenodo.1065400). These
fields are available at approximately 0.75-degree spatial resolution for 1980–2018.

### Wind
Daily 10-m wind components (u10, v10) at 12:00 UTC were obtained from the ERA5 reanalysis
(Copernicus Climate Change Service, 1980–2018) at approximately 0.25-degree resolution.
Wind speed magnitude was computed from vector components and converted from m/s to km/h
prior to all calculations, consistent with the SCF climate library convention (SCF User
Guide, §2.43).

### Fire Ignitions
Lightning and human-caused (accidental) fire ignitions were obtained from the Fire Program
Analysis Fire-Occurrence Database (FPA FOD; Short 2022), which records all fires reported
to US federal, state, and local agencies from 1992–2020. Discovery date and cause class
were used to assign each ignition to a calendar day and ignition type.

### Fire Perimeters
Historical fire perimeter data were sourced from one of two products depending on the
calibration objective:

- **GeoMAC Historic Perimeters (2000–2018):** Daily or near-daily cumulative fire
  perimeters mapped by incident management teams using aerial infrared and GPS surveys.
  Sequential perimeter pairs represent genuine cumulative growth and provide accurate
  fire-edge geometry for spread probability calibration.

- **FIRED (Fire Events Delineation; Balch et al. 2020):** Daily burn detection patches
  derived from the MODIS MCD64A1 500-m burned area product. Because MODIS detections
  are incremental (each record represents newly detected burn on a given day rather than
  cumulative area), a cumulative reconstruction step was applied before use (see below).

### Terrain
A digital elevation model (DEM) covering the calibration domain was obtained from the
SRTM/3DEP dataset via the elevatr R package. Slope (degrees) and aspect (degrees from
north, clockwise) were derived from the DEM using standard first-difference gradient
methods and resampled to the LANDIS grid resolution.

### Fine Fuels
Fine fuel load was represented using the LANDFIRE Fire Behavior Fuel Model 40-class
(FBFM40) dataset, normalized to a 0–1 index. Where the LANDFIRE product was unavailable,
a spatially uniform placeholder value of 1.0 was used, effectively setting the fine fuels
coefficient (B2) to zero in the fitted model.

---

## Calibration Domain and Spatial Framework

All data were clipped to a large regional calibration boundary (e.g., the WA Cascades,
Pacific Northwest states) chosen to encompass the range of fire weather conditions and
fuel types present in the target LANDIS landscape. ERA data were spatially extracted to
this boundary; perimeter and ignition records were filtered to intersect it. The working
coordinate reference system (CRS) was taken from the LANDIS template raster to ensure
all outputs align with the simulation grid.

---

## 1. Ignition Calibration

### Objective
Calibrate daily ignition rates for lightning and accidental fire types as a function of
the Fire Weather Index.

### Data Preparation
FPA FOD ignition records were filtered to the calibration boundary and classified by
ignition cause (NWCG cause code: lightning vs. human). Records were aggregated to daily
ignition counts by type and joined to the ERA-Interim FWI time series. Days with no
recorded ignitions were retained as zero counts, producing a complete daily time series
suitable for count regression.

### Model Fitting
Separate count-regression models were fit for lightning and accidental ignition types.
Two model families were evaluated:

**Poisson regression:**
```
lambda = exp(b0 + b1 * FWI)
```
The expected daily count follows a Poisson distribution with mean lambda. Daily draws
during simulation use this rate parameter directly.

**Zero-Inflated Poisson (ZIP, optional):**
When the proportion of zero-count days exceeded what a standard Poisson model predicts,
a zero-inflated Poisson model was fit using the pscl package. The ZIP model augments the
Poisson rate equation with a logistic component (BZ0, BZ1) controlling the probability
of a structural zero (no fire possible regardless of FWI).

### SCF Parameter Mapping
| Calibrated coefficient | SCF parameter |
|---|---|
| b0 (lightning) | LightningIgnitionB0 |
| b1 (lightning) | LightningIgnitionB1 |
| b0 (accidental) | AccidentalIgnitionB0 |
| b1 (accidental) | AccidentalIgnitionB1 |
| bz0, bz1 (ZIP only) | LightningIgnitionBZ0/BZ1, AccidentalIgnitionBZ0/BZ1 |

### Ignition Surface
The spatial allocation of ignitions across the LANDIS grid was represented as a kernel
density surface computed from the historical FPA FOD point locations. Each ignition point
was rasterized onto the LANDIS template grid and smoothed using an isotropic Gaussian
kernel via iterative focal-mean operations. The resulting surface represents the relative
probability that an ignition, when it occurs, will start in each cell.

### Climate Startup Values
ERA-Interim climatology was used to compute representative startup values for the FFMC,
DMC, and DC moisture codes at the beginning of each simulated fire season, ensuring that
the FWI trajectory during simulation begins from a physically realistic moisture state.

---

## 2. Spread Probability Calibration

### Objective
Calibrate the probability that fire advances from a currently burning cell to an adjacent
unburned cell as a function of FWI, fine fuel load, and effective wind speed.

### SCF Formulation
SCF evaluates the spread decision (CanSpread function, FireEvent.cs) as:

```
P(spread) = 1 / (1 + exp(-(B0 + B1*FWI + B2*FineFuels + B3*EWS)))
```

Calibration produces logistic regression coefficients that enter the SCF parameter file
directly without transformation.

### Perimeter Pair Construction

**GeoMAC data:** Because GeoMAC perimeters record cumulative burned area, sequential pairs
(perimeter at time t, perimeter at time t+1) directly represent the state of a fire at
the start and end of each time step. Pairs were filtered to a maximum gap of 1 day for
spread probability fitting to ensure that each pair captures a single day of fire
behavior.

**FIRED cumulative reconstruction:** Because FIRED records individual daily burn detection
patches rather than cumulative area, a reconstruction step was applied before building
pairs. For each fire event, daily patches were sorted by detection date and a running
spatial union was computed, producing a synthetic cumulative perimeter series equivalent
in structure to GeoMAC. Quality filters applied prior to reconstruction included a minimum
total fire size threshold (default 100 ha), a minimum number of distinct detection days
(default 3 days), and a maximum daily growth factor to exclude pairs where the area
increment was implausibly large relative to the existing fire size (a signature of
multi-day MODIS detection gaps compressed into a single record).

Pairs with more than 1 day between perimeter dates, with negative area growth exceeding
a tolerance threshold, or flagged by the growth factor filter were excluded from spread
probability fitting.

### Cell Classification
Each perimeter pair was rasterized to the spread calibration grid. Individual grid cells
were classified into four categories:

- **Already burned (t):** Within the cumulative perimeter at time t.
- **Newly burned — proximal:** Within the t+1 perimeter but not t, AND within a
  circular dilation (user-specified radius, default 5 cells) of the t perimeter edge.
  These cells are the primary "success" cases for the logistic model.
- **Newly burned — distant:** Within the t+1 perimeter but not t, AND outside the
  dilation ring. These represent potential spotting or data artifacts and are included
  in model fitting but tracked separately in diagnostics.
- **Frontier fail:** Within the dilation ring of the t perimeter but not in the t+1
  perimeter. These are the "failure" cases — cells that were adjacent to the active fire
  edge but did not burn.

A cumulative burn accumulator ensured that the "already burned" category grew
monotonically across all frames within a fire event, filling detection gaps that would
otherwise cause the burned footprint to shrink between pairs.

### Effective Wind Speed
Effective Wind Speed (EWS) was computed following Nelson (2002, Int. J. Wildland Fire,
Eq. 5), which combines ambient wind speed and terrain slope into a single aerodynamic
forcing term:

```
EWS = Cb * sqrt((Ua/Cb)^2 + 2*(Ua/Cb)*sin(theta)*cos(phi) + sin^2(theta))
```

where Ua is ambient wind speed (km/h), theta is slope angle (radians), phi is the angle
between wind direction and the uphill azimuth (radians), and Cb is a combustion buoyancy
parameter (10 for low-intensity, 25 for moderate, 50 for high-intensity fire). The uphill
azimuth was derived as (terrain aspect + 180) mod 360. ERA5 wind direction was computed
from the u10/v10 vector components.

### Model Fitting
A binary logistic regression (binomial GLM, logit link) was fit to the rasterized
cell-level success/failure samples pooled across all eligible perimeter pairs. Predictors
were FWI (at the centroid of each pair's t perimeter), FineFuels (LANDFIRE FBFM40
normalized 0–1, spatially sampled at each cell), and EWS (from the Nelson formula above).
FineFuels was dropped from the model and B2 set to zero when the fine fuels raster was
unavailable or when the variance of the fuel index across samples was below a
detectability threshold (< 0.001). A configurable failure-to-success sampling ratio
(default 3:1) was used to balance the logistic regression without requiring that all
frontier-fail cells be retained.

A VIF (Variance Inflation Factor) diagnostic was computed for all predictor pairs.
Multiple model variants (FWI only; FWI + EWS; FWI + FineFuels + EWS) were evaluated
using AUC and AIC to select the best-supported predictor combination.

### SCF Parameter Mapping
| Calibrated coefficient | SCF parameter |
|---|---|
| B0 (intercept) | SpreadProbabilityB0 |
| B1 (FWI) | SpreadProbabilityB1 |
| B2 (FineFuels) | SpreadProbabilityB2 |
| B3 (EWS) | SpreadProbabilityB3 |

Note: SCF applies a 0.71× weight to diagonal-direction spread internally; this does not
affect the calibrated coefficient values.

---

## 3. Maximum Daily Spread Area Calibration

### Objective
Calibrate the upper bound on area (ha) that a single fire can burn in one day as a
function of FWI and effective wind speed.

### SCF Formulation
SCF's Spread function (FireEvent.cs) advances the simulation clock to the next day when
cumulative daily area exceeds:

```
MaxSpreadArea = (int)(B0 + B1 * FWI + B2 * EWS)   [hectares]
```

Because the coefficients enter SCF directly without transformation, the model must be fit
on the raw hectare scale with a strictly linear (no log-transform) formulation.

### Data Preparation
Daily area increments were computed from either sequential GeoMAC perimeter pairs or the
FIRED events product. For GeoMAC pairs, the daily increment equals the difference in
cumulative area between consecutive perimeters, divided by the gap in days. For the FIRED
events product, the maximum single-day burned area within each fire event was used
directly, paired with the FWI and EWS on the date of maximum growth.

Pairs and events were filtered to a maximum gap of 3 days, positive area increments only,
and the optional maximum growth factor filter described above.

### Calibration Target
MaxSpreadArea is a ceiling, not a mean: calibrating to the mean daily increment would
cap roughly half of all observed fire-days, producing unrealistically slow simulated
spread. Models were therefore fit at three quantile targets — the median (50th
percentile), the 75th percentile, and the 90th percentile — and the 75th-percentile fit
is recommended as the default starting value for SCF.

### Model Fitting
Quantile regression was used as the primary fitting method, targeting each quantile (tau)
separately:

**Unconstrained (default):** Standard quantile regression via quantreg::rq() without
inequality constraints. Coefficients are free to take any sign; the user should inspect
signs and apply constraints if collinearity produces wrong-sign coefficients.

**Constrained (optional):** Quantile regression with explicit inequality constraints
enforcing B0 ≥ 50 ha, B1 ≥ 0, B2 ≥ 0. This prevents SCF "MaxSpreadArea < 0" warnings
at low FWI/EWS values. Fallback methods in priority order: constrained quantile regression
(quantreg, br method with R/r matrix), non-negative least squares with reparameterization
(nnls), and weighted OLS with post-hoc coefficient clamping.

Collinearity between FWI and EWS was assessed using the condition number of the design
matrix and VIF for each predictor. EWS can be excluded from the model when FWI/EWS
correlation is strongly negative.

### SCF Parameter Mapping
| Calibrated coefficient | SCF parameter |
|---|---|
| B0 (intercept) | MaximumSpreadAreaB0 |
| B1 (FWI) | MaximumSpreadAreaB1 |
| B2 (EWS) | MaximumSpreadAreaB2 |

---

## References

Balch, J.K., et al. (2020). FIRED (Fire Events Delineation): An open, flexible algorithm
and database of US fire events derived from the MODIS burned area product (2001–2019).
Remote Sensing, 12(21), 3498.

Nelson, R.M. (2002). An effective wind speed for models of fire spread. International
Journal of Wildland Fire, 11(2), 153–161.

Short, K.C. (2022). Spatial wildfire occurrence data for the United States, 1992–2020
(FPA FOD). Fort Collins, CO: Forest Service Research Data Archive.
https://doi.org/10.2737/RDS-2013-0009.6

Vitolo, C., et al. (2019). ERA-Interim based global fire weather indices [Data set].
Zenodo. https://doi.org/10.5281/zenodo.1065400
