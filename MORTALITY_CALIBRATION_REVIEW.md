# SCF Mortality Sub-model — Source Audit and Calibration Review

Source reviewed: `LANDIS-II-Foundation/Extension-Social-Climate-Fire` @ master (SCF v4.3)
Files: `src/FireEvent.cs`, `src/InputParameterParser.cs`, `src/SiteVars.cs`,
`src/TopographySoils.cs`, `src/SpeciesData.cs`, `docs/LANDIS-II Social-Climate-Fire v4 User Guide.pdf`,
`Supporting R Code/Mortality/Calculating_Mortality.Rmd` (Robbins 2021 reference implementation).

App code reviewed: `modules/mortality_helpers.R`, `modules/species_helpers.R`,
`modules/spread_helpers.R` (fine fuels), `app.R` (mortality tab server block).

> **STATUS (2026-08-21): implemented.** All recommendations below have been applied
> except the continuous-RdNBR data switch (#3), which is now *supported* by the code
> but still needs you to point the app at per-fire MTBS rasters. See
> "What changed" at the end of this document for the as-built summary.

---

## Part 1 — What SCF actually does

Mortality is a strict two-stage cascade. Stage 1 produces a *site* severity value on the
dNBR scale; Stage 2 uses that value, plus species bark thickness, to kill cohorts.

### Stage 1 — Site mortality (`FireEvent.CalculateDNBR`, FireEvent.cs:312–412)

Called once per burned cell, immediately after `CalculateEffectiveWindSpeed`, from the
spread loop at FireEvent.cs:186.

```
eta = B0
    + B1 * Clay                  // SiteVars.Clay, from ClayMap raster
    + B2 * Previous_Year_PET     // SiteVars.PotentialEvapotranspiration = "Succession.PET"
    + B3 * siteEffectiveWindSpeed
    + B4 * WaterDeficit          // SiteVars.ClimaticWaterDeficit = "Succession.CWD"
    + B6 * ladderFuelBiomass
    + B5 * fineFuelPercent

siteMortality = max(eta, 0.0005) ^ -1      // Gamma inverse link; hard ceiling 2000
siteMortality = max(siteMortality, 0.0)
this.SiteMortality = (int) siteMortality
SiteVars.DNBR[site] = this.SiteMortality
```

Details that matter for calibration:

| Term | What SCF feeds it | Units / range |
|---|---|---|
| `Clay` | `ClayMap` raster, read verbatim | User Guide §2.18: **0–1 fraction** |
| `Previous_Year_PET` | NECN `Succession.PET`; falls back to `TimeZeroPET` at year 0 | mm/yr, NECN's own PET formulation |
| `WaterDeficit` | NECN `Succession.CWD`; falls back to `TimeZeroCWD` | mm/yr |
| `siteEffectiveWindSpeed` | Nelson (2002) Eq. 5, day-of-burn wind, **dynamic Cb** | km/h |
| `fineFuelPercent` | `min(SiteVars.FineFuels / MaximumFineFuels, 1.0)` | **0–1**, rescaled by the parameter |
| `ladderFuelBiomass` | `Σ cohort.Data.Biomass` for species in `LadderFuelSpeciesList` with `Age <= LadderFuelMaxAge` | **raw g m⁻², unbounded (0 – several thousand)** |

Two consequences of the inverse link:
- **Coefficient signs are inverted relative to intuition.** A variable that *increases*
  severity must have a *negative* coefficient. B0 must be positive and small
  (B0 = 3e-3 → baseline dNBR ≈ 333). The reference implementation confirms this:
  B0 = +2.9e-3, B3(EWS) = −3.0e-5, B5(fuels) = −1.6e-6.
- **The clamp at `eta = 0.0005` is a cliff, not a soft bound.** Any cell whose linear
  predictor goes to zero or negative is assigned dNBR = 2000 — maximum severity.
  Miscalibrated coefficients don't degrade gracefully; they saturate the landscape.

Side effects of `siteMortality` beyond cohort killing:

- `SiteVars.Intensity = min(max(siteMortality/100, 1), 10)`
- That intensity feeds back into **combustion buoyancy** on subsequent spread
  (FireEvent.cs:606–610): `Cb = 10`, → `25` if source-site Intensity > 3,
  → `50` if > 6. So severity raises EWS, and (with negative B3) EWS raises severity.
  There is a positive feedback loop between Stage 1 and the spread engine.
- Rx fires: intensity clamped by `RxMaxFireIntensity`, and `SiteMortality` is *snapped*
  to `standardSeverityIndex * 100` — coarse 100-unit bins.

### Stage 2 — Cohort mortality (`ReduceOrKillMarkedCohort`, FireEvent.cs:416–480)

```
BarkThickness = MaximumBarkThickness[spp] * Age / (Age + AgeDBH[spp])    // Eq. 10
P  = exp(B0 + B1*BarkThickness + B2*SiteMortality)
p  = P / (1 + P)                                                         // logit
if (p > runif()) -> kill                                                 // returns 1.0
```

- `MaximumBarkThickness` and `AgeDBH` come from `Species_CSV_File`.
- `SiteMortality` here is the **integer dNBR value from Stage 1** (0–2000), not a 0–1 index.
- Expected signs: **B1 negative** (thicker bark → survives), **B2 positive**.
- Returning `1.0` removes the **entire cohort's biomass**. There is no partial mortality —
  it is one Bernoulli draw per cohort per fire.
- User Guide §Eq. 10 and the code comment both say `MaximumBarkThickness` is
  "the bark thickness of the largest DBH observed within the species", derived by
  "linear scaling from the Fire Tree Mortality Database (Cansler et al. 2020)" —
  i.e. `BT_coef × maxDBH`. **Cansler's BT_coef is single-sided bark thickness.**

### Dead code

`src/FireDamages.cs` (`IFireDamage` / `FireDamage`, species × age-range × fixed mortality
probability) is **not referenced anywhere** in the v4 source. It is a v1/v2 leftover.
Ignore it — the severity path is entirely Equations 7–10.

---

## Part 2 — Review of our calibration

### ✅ What is already correct

- Gamma GLM with inverse link, coefficients transferred raw — matches SCF exactly, no
  back-transformation needed. (`fit_site_mortality`)
- Clay as a 0–1 fraction — matches User Guide §2.18. (`process_clay_raster`)
- **PET** (TerraClimate `pet`), not AET, for B2 — matches `Succession.PET`. Correctly
  called out in the code comment; a common mistake avoided.
- CWD from TerraClimate `def`, previous year — matches `Succession.CWD`.
- Fixed B0–B6 slot mapping with zero-padding for unfitted terms — correct, since SCF
  reads positional keywords and renumbering would silently mis-assign coefficients.
- Wind in km/h throughout, Nelson Eq. 5 shared with the spread module.
- Binomial logit for cohort mortality with `dead ~ BarkThickness + dNBR` — structurally
  identical to the reference implementation.
- Predictor standardisation with correct algebraic unscaling of both slopes and intercept.

### 🔴 Blocking defects

**1. Ladder fuels are on the wrong scale by ~3 orders of magnitude.**
`fetch_landfire_cbh()` builds a **0–1 index** (`1/(CBH_m+0.5)`, normalised by its own max).
SCF feeds B6 **raw summed cohort biomass in g m⁻²**, typically 0–3000.

Illustration: suppose the 0–1 index yields a plausible B6 ≈ −1e-3. At runtime with
`ladderFuelBiomass = 800 g m⁻²`, the B6 term contributes −0.8 to a linear predictor
whose intercept is ~+3e-3. `eta` goes hard negative → clamped to 0.0005 →
**dNBR = 2000 in every cell with any ladder fuel**. The landscape saturates at maximum
severity and every cohort dies.

This is the single most consequential finding. Options, in order of preference:
- (a) Build the ladder-fuel predictor as *biomass*, not an index: convert LANDFIRE CBH
  (or better, canopy bulk density × the sub-`LadderFuelMaxAge` canopy layer) into
  g m⁻², or reconstruct it from a LANDIS spin-up biomass map filtered to the species/age
  set you intend to put in `LadderFuelSpeciesList`. The latter is more faithful, because
  SCF computes it from its own cohort list, not from remote sensing.
- (b) If (a) isn't feasible, **leave B6 = 0** and say so explicitly in the snippet. A
  zeroed B6 is far safer than a mis-scaled one, and Eq. 7 in the User Guide doesn't even
  include the ladder term.
- The app currently defaults to a **0.5 placeholder** when no ladder raster is loaded.
  That is worse than NA: a constant predictor with zero variance still consumes an
  intercept-shifting coefficient. It should be dropped from the model, not filled.

**2. `MaximumBarkThickness` is double-sided in our species table, single-sided in the calibration.**
`species_helpers.R:309` computes `MaximumBarkThickness = 2 * bark_ratio * p95_dia_cm`
("double-sided"). But `load_ftm_cohort_data()` fits B1 against
`BarkThickness = BT_coef * DBH_cm` — Cansler's single-sided convention — and the User
Guide/`FireEvent.cs` comment specify single-sided as well.

So SCF applies a B1 fitted on single-sided bark to a bark thickness that is **2× too large**.
For Douglas-fir (BT_coef ≈ 0.048, P95 DBH ≈ 60 cm): calibration BT ≈ 2.9 cm, runtime
BT ≈ 5.8 cm. At the reference B1 = −0.79 that is a logit shift of −2.3, an **~10× reduction
in mortality odds** for mature cohorts. Fire will barely kill anything.

Fix: drop the `× 2`, or equivalently double the FTM bark thickness in the calibration.
Whichever you choose, it must be the same convention in both places. My recommendation is
single-sided in both, since that is what the User Guide specifies.

Also worth reconciling: the reference implementation sets `MaxBark = Dmax × BT_coef` where
`Dmax` is the **fitted asymptote** of the DBH~Age curve, and `AgeDBH` is the
half-saturation constant from the *same* fit. Our code fixes `Dmax = P95 DBH` and fits only
`theta`, which is internally consistent — fine — but the P95 choice should be stated in the
snippet since it sets the bark ceiling.

**3. The site-mortality response is a 3-level categorical, not a continuous dNBR.**
`fetch_mtbs_from_imageserver()` pulls the MTBS **thematic 6-class** CONUS mosaic and
`subst()`s it to three midpoints:

```
class 1 (unburned/low)   -> NA      <- low-severity tail deleted
class 2 (low)            -> 185
class 3 (moderate)       -> 355
class 4 (high)           -> 600
class 5 (incr. greenness)-> NA
class 6 (non-mapping)    -> NA
```

Fitting `Gamma(link="inverse")` to a variable that takes exactly three values means:
- the Gamma variance structure is fiction;
- coefficient magnitudes are driven by the *class boundaries*, not by real severity gradients;
- the predicted dNBR range at runtime will be badly compressed (see #4);
- dropping class 1 truncates the left tail, so the fitted mean severity is biased high.

The reference implementation uses **continuous RdNBR** (their example spans ~200–1400,
mean 548). The User Guide (§Eq. 7) says "gridded observation of fire severity (such as
those available in the dNBR style maps)" and §2.56–2.58 say "DNBR, or CBI".

Recommendation: switch the primary source to **continuous per-fire MTBS rasters**
(`*_rdnbr.tif`, or `*_dnbr.tif` from the MTBS fire bundles) rather than the thematic
CONUS mosaic. Prefer **RdNBR** — it's what Robbins used, and it's less sensitive to
pre-fire biomass, which matters when you're regressing against fuels. Keep the
class-midpoint path only as an explicitly-labelled fallback, and warn loudly when it's active.

Two related issues in the same function:
- `min_dnbr = 100` is a hard-coded left truncation that is never exposed in the UI. With
  continuous data this systematically inflates the intercept and biases every slope. It
  should be a UI control defaulting to something much lower (e.g. 1, just enough to keep
  the Gamma support positive), and any truncation should be reported.
- MTBS continuous products carry a per-fire **dNBR offset** in their metadata that should
  be subtracted before use. Worth at least checking against a couple of fires.

**4. Cohort-model dNBR is one value per fire, assigned to thousands of trees.**
`load_ftm_cohort_data()` takes the **median dNBR in a 1 km buffer around the fire's
lat/lon** and joins it to every tree in that fire. Combined with #3, most fires collapse to
a single value (usually 355), so B2 — the entire link between site severity and cohort
death — is estimated from between-fire variation across a handful of distinct numbers.

Two separate problems:
- **Identifiability.** With continuous RdNBR this improves, but FTM only publishes
  fire-level coordinates, so within-fire severity variation is genuinely unavailable.
  Realistically B2 will always be a between-fire estimate. That should be stated in the
  snippet rather than implied to be pixel-level.
- **Pseudoreplication.** Thousands of trees share one dNBR value, so the reported standard
  errors and p-values on B2 are wrong by a large factor. Fit with cluster-robust SEs by
  `YrFireName` (`sandwich::vcovCL`) or a fire-level random intercept
  (`lme4::glmer(dead ~ BarkThickness + dNBR + (1|YrFireName))`, taking the fixed effects).
  The point estimate matters more than the SE for SCF, but the app currently presents a
  p-value that will read as wildly significant regardless of the evidence.

Minor, same function: the 1 km buffer is centred on the fire's reported coordinate, which
for many FTM fires is a centroid nowhere near the sampled plots. Consider widening the
buffer and reporting the within-buffer dNBR IQR as a diagnostic, so you can see how much
the median is actually pinning down.

### 🟡 Consistency issues worth fixing

**5. `MaximumFineFuels` rescales B5 too, not just `SpreadProbabilityB2`.**
User Guide §2.28: "The value rescales the fine fuel parameter in equations 5 **and 7**".
Our `CLAUDE.md` documents this linkage for spread but the mortality snippet does not
mention it at all. Since our fine fuels are FBFM40 normalised by 9.0 t/ac (≈1000 g C m⁻²),
`MaximumFineFuels` must be set to a comparable load or B5 is on the wrong scale — exactly
the same trap already documented for B2. The mortality snippet should carry the same
warning and echo the currently-set `MaximumFineFuels`.

**6. Calibration EWS uses an annual fire-season mean; SCF uses day-of-burn wind.**
`extract_site_predictors()` calls `compute_annual_fire_season_wind()` (Jun–Sep mean) and
joins one wind speed/direction per *year*. The reference implementation says explicitly
"the effective wind speed **on the day of**", and our own validation tab (`validation_helpers.R:1310`)
already computes EWS from wind on the MTBS ignition date.

Effect: annual means have a fraction of the variance of daily values, so B3 is estimated
over a narrow x-range and then extrapolated to the full daily distribution at runtime.
Given the Stage-1 ↔ buoyancy feedback loop, an over-steep B3 is how you get a landscape
that flips to dNBR = 2000. The fix is small: use the MTBS per-fire ignition date (already
available in the validation path) instead of the annual summary.

**7. Combustion buoyancy is fixed at Cb = 10 in calibration, dynamic in SCF.**
SCF escalates Cb to 25 and 50 once site intensity exceeds 3 and 6. Cb = 10 is a defensible
calibration choice (it's the low-intensity default and matches our spread module), but B3
will then be applied to EWS values that are systematically larger in exactly the severe
fires that matter. At minimum this should be documented in the snippet; a sensitivity
readout of B3 refit at Cb = 25 would be a cheap and informative addition.

**8. PET/CWD formulation mismatch between TerraClimate and NECN.**
TerraClimate `pet` is Penman-Monteith reference ET (~1000–1400 mm/yr in the Cascades);
NECN's `Succession.PET` uses its own formulation and typically runs substantially lower.
B2 and B4 are calibrated in TerraClimate units and applied in NECN units. Since these
enter a linear predictor whose intercept is ~3e-3, even a 30% offset is a meaningful
shift in baseline severity.

Recommended check: after a short LANDIS test run, read SCF's own PET and CWD output maps
(`PlugIn.cs:671, 696` write them) and compare the landscape mean against the TerraClimate
mean over the same footprint. If they differ materially, either rescale before fitting or
expose a scale factor. This is a *verification step*, not a code change — worth doing once.

**9. `TimeZeroPET` / `TimeZeroCWD` are never produced by the app.**
SCF uses these at year 0 before NECN has produced PET/CWD (FireEvent.cs:344–360). If they
default to 0.0, the first simulated fire year gets a linear predictor missing two terms —
which, with negative B2/B4, means an artificially *low* eta and therefore artificially
*high* severity. The app already has the landscape climatology; emitting the landscape-mean
PET and CWD into the snippet is a few lines and closes a real gap.

**10. `LadderFuelMaxAge` and `LadderFuelSpeciesList` are not produced at all.**
Both are required for B6 to have any meaning. If we pursue the biomass-based ladder fuel
fix (#1a), the species list and max age used to *build* the predictor must be the same
ones written to the parameter file.

### 🟢 Diagnostics I'd add before trusting any of these numbers

**A. Linear-predictor positivity check.** After fitting, evaluate `eta` across the observed
predictor ranges (and ideally across the 1st–99th percentile of each predictor
independently) and report the fraction of cells where `eta <= 0.0005`. Anything above ~0%
means the clamp is firing and severity is saturating. This single check would have caught
the ladder-fuel scale bug immediately.

**B. Predicted-vs-observed dNBR *range*, not just fit.** A GLM's fitted values compress
variance. If the site model only spans, say, 300–450 dNBR while MTBS spans 100–1200, then
the cohort model effectively sees a constant `SiteMortality` at runtime and B2 does nothing
— severity becomes species-driven only. Plot both distributions side by side. Given issue
#3, I'd expect the current model's range to be very narrow.

**C. Runtime mortality curve preview.** Combine the two stages and plot P(mortality) vs
cohort age for each species at the 10th/50th/90th percentile of predicted site dNBR, using
the *actual* `AgeDBH`/`MaximumBarkThickness` values that will be written to
`Species_CSV_File`. This is the reference implementation's key figure and it would surface
the bark-thickness convention bug (#2) visually.

**D. Sign guard.** Warn when B0 ≤ 0, when any of B1–B6 has a sign implying the wrong
direction of effect, or when B1(cohort) > 0 / B2(cohort) < 0.

### Minor / cosmetic

- `fit_site_mortality()` filters `dNBR > 0`; with a proper continuous product this drops
  genuine low-severity pixels. Combined with `min_dnbr = 100` the truncation is double.
- Stale comment in `plot_yr1status_distribution` (`mortality_helpers.R` ~line 1727) refers
  to "yr1status=2 (confirmed dead)", contradicting the loader's documented 0/1 coding.
- `write_mortality_snippet()` labels B6 as "ladder fuel index 0-1" — accurate to what we
  compute, but that *is* the bug; once fixed the note should read "g m⁻² summed cohort
  biomass".
- The snippet header writes Eq. 7 as `dNBR = 1 / (B0 + ...)` without noting the
  `max(eta, 0.0005)` clamp or the 2000 ceiling. Worth adding.

---

## Suggested order of work

1. **Bark thickness convention** (#2) — one-line change, large effect, no new data needed.
2. **Ladder fuels** (#1) — decide between biomass-based predictor and B6 = 0. Until then,
   stop filling the 0.5 placeholder.
3. **Positivity + range diagnostics** (A, B) — cheap, and they tell you whether anything
   else you do is working.
4. **Continuous RdNBR source** (#3) — the biggest data-side change, gates real improvement
   in B0/B3/B5 and in cohort B2.
5. **Day-of-burn EWS** (#6) — reuse the validation-tab path.
6. **Cluster-robust SEs on the cohort fit** (#4).
7. **Snippet completeness**: `MaximumFineFuels` warning, `TimeZeroPET`/`TimeZeroCWD`,
   `LadderFuelMaxAge`/`LadderFuelSpeciesList`, clamp/ceiling documentation (#5, #9, #10).
8. **PET/CWD scale verification against a LANDIS test run** (#8) — do once, document result.

---

## Appendix — Source provenance

There is a second SCF clone on this machine at
`~/Documents/GitHub/Extension-Social-Climate-Fire` (HEAD `4b089af`, "Updating csproj with
correct hint path"). This review was written against a **fresh clone of GitHub master**
(`0babccb`), so I diffed the two.

**The mortality logic is identical.** Diffing every line of `CalculateDNBR` and
`ReduceOrKillMarkedCohort` that touches the equations — `siteMortality`, `Pm`,
`probabilityMortality`, `BarkThickness`, `fineFuelPercent`, `ladderFuelBiomass`,
`standardSeverityIndex` — yields a single difference: one commented-out debug `WriteLine`.
The mortality parameter keywords in `InputParameterParser.cs` are identical as well.
Every finding above applies to both copies.

**But the local clone is stale, and it matters for spread, not mortality.** Its
`origin/master` ref still points at its own HEAD, and the current GitHub master commit
doesn't exist in it. Concretely, the local copy has **no `BurningSitesThreshold`** —
zero occurrences in either `FireEvent.cs` or `InputParameterParser.cs`. Our `CLAUDE.md`
documents that feature at `FireEvent.cs:260-272, 292-310` and warns it "substantially
raises effective spread and invalidates a B0 tuned without it."

So: if the LANDIS build being used for spread calibration came from that local clone,
`BurningSitesThreshold` is not available regardless of what the parameter file says, and
the `CLAUDE.md` line numbers for the spread notes won't match the local source. Worth a
`git pull` there before the next build, and worth confirming which source tree produced
the binary currently in use.


---

## What changed (2026-08-21)

### `modules/species_helpers.R`
- **Removed the `× 2`.** `MaximumBarkThickness = bt_coef × P95_DBH_cm`, single-sided.
- **Replaced the hard-coded bark table with FOFEM/Cansler BT_coef values.** The old
  table disagreed with Cansler's published coefficients by 0.72×–3.67× with no
  consistent factor, so removing the ×2 alone would not have made the species table
  and the cohort GLM agree. Congener substitutions for the four species FOFEM does
  not list (*Abies procera*, *Betula papyrifera*, *Larix lyallii*,
  *Populus balsamifera*) are recorded in a new `bt_note` column.
- **Added `load_ftm_bt_coef()` and `resolve_ftm_data_dir()`.** When the FTM directory
  is set, `estimate_bark_params_from_fia()` reads BT_coef straight from the user's own
  `Species_BarkThickness.csv`, so the species table and `fit_cohort_mortality()` use
  literally the same numbers. The built-in table is now a fallback.

### `modules/mortality_helpers.R`
- **`fetch_landfire_cbh()` → `fetch_landfire_ladder_fuels()`.** Now fetches CBH, CBD
  and CH and returns **g m⁻²**:
  `CBD_kg_m3 × max(0, min(CH, ladder_ceiling_m) − CBH) × 1000`.
  Verified on synthetic LANDFIRE-native input: 0–1625 g m⁻², mean ~700, and 0 where
  CBH sits above the ladder ceiling.
- **`load_mtbs_dnbr_pixels()`**: `min_dnbr` default 100 → 1 (with a warning above 50);
  MTBS class 1 now maps to dNBR 60 instead of being deleted; carries a per-fire
  `fire_date`; added a `dnbr_offset` argument; warns when a thematic product is
  detected or when fewer than 25 distinct dNBR values survive.
- **`extract_site_predictors()`**: uses **day-of-burn wind** when per-fire dates are
  available (detects and rejects the 1-Jan annual-mosaic bound), falling back to the
  annual mean with a warning; no more constant placeholders — absent predictors stay
  NA so they are dropped and written as 0.0; warns when ladder values look like an
  index; records `ews_wind_basis` and `combustion_buoyancy` as attributes.
- **`fit_site_mortality()`**: drops zero-variance predictors; adds clamp-exposure,
  predicted-vs-observed range, and sign-guard diagnostics; returns them for the snippet.
- **`cluster_robust_vcov()`** (new): CR0 / Liang–Zeger with the standard finite-sample
  correction, implemented directly so no `sandwich` dependency is needed. Validated to
  5e-17 against an independent reference; on data with fire-level heterogeneity it
  inflates the SE on the severity coefficient ~5×.
- **`fit_cohort_mortality()`**: cluster-robust SEs by fire, sign guard, and a warning
  when fewer than 5 distinct dNBR values are available to B2.
- **`write_mortality_snippet()`**: rewritten. Documents the clamp and inverted signs,
  emits `MaximumFineFuels`, `TimeZeroPET`, `TimeZeroCWD`, `LadderFuelMaxAge` and
  `LadderFuelSpeciesList`, and prints the diagnostics.

### `app.R`
- Ladder fuel source relabelled to biomass; added `ladder_ceiling_m`,
  `ladder_max_age` and `ladder_species` inputs; "None" now means B6 = 0 rather than a
  0.5 placeholder; the status badge turns red if the raster looks like an index.
- Exposed `mtbs_min_dnbr`; passes `wind_daily` into predictor extraction; passes the
  companion parameters into the snippet; new inputs registered for session save/restore.
- Species tab text corrected: FIA supplies no bark thickness, and MaximumBarkThickness
  is single-sided.
- Fixed the stale `yr1status=2` comment (FTM codes 0 = alive, 1 = dead).

### Verification
- Site GLM recovers all seven known coefficients on synthetic data at SCF scales
  (B0 3.0e-3 → 2.98e-3; B3 −2.0e-5 → −1.86e-5; B6 −8.0e-7 → −7.93e-7), confirming the
  unscaling algebra.
- Guard rails exercised: constant/missing predictors dropped and written as 0.0;
  wrong-sign coefficients flagged; index-scale ladder fuels flagged; annual-mean wind
  fallback warns.
- All seven modules source cleanly and `app.R` parses.

### One bug found while testing
`extract_site_predictors()` initially used dot-prefixed temporary column names
(`.d`, `.ws_yr`). In dplyr 1.1.4 a `mutate()` that creates a dot-prefixed column
breaks the data mask for bare column references in the same call, so `as.Date(date)`
silently resolved to `base::date` and errored at runtime. Renamed to `burn_day` /
`ws_yr` / `wd_yr`. A scan found no other instances in the codebase.

### Still open — your call
- **Continuous RdNBR (#3).** The loader now handles continuous products, exposes the
  truncation, and warns loudly when it sees a thematic one, but the app still defaults
  to the MTBS thematic CONUS mosaic. Pointing it at per-fire MTBS rasters is a data
  step, and RdNBR is what Robbins used.
- **PET/CWD scale check (#8).** Needs one LANDIS test run to compare SCF's own PET/CWD
  output maps against the TerraClimate means. Verification, not a code change.
