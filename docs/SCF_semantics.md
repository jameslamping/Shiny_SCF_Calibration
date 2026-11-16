# How SCF Actually Uses These Parameters

What the Social-Climate-Fire (SCF v4) extension does with each calibrated
parameter at runtime, and the traps that follow from it. Read this before
trusting any coefficient the app produces.

Source citations are to `LANDIS-II-Foundation/Extension-Social-Climate-Fire`
master at commit `4b089af` (2026-04-15): `src/FireEvent.cs`, `src/PlugIn.cs`,
`src/InputParameterParser.cs`. Line numbers move between versions; search for
the quoted code if they don't match yours.

---

## 1. Parameter file syntax and order

SCF's parameter file is read by `InputParameterParser.cs` with the standard
LANDIS-II text parser.

- **Comments:** whole-line comments start with `>>`. Trailing comments start
  with `<<`. A `#` is **not** a comment.
- **Values:** each parameter is `Keyword  value`, with no `=` (it would be read
  as the value).
- **Keywords are plural:** `LightningIgnitionsB0`, not `LightningIgnitionB0`.
  The C# properties are singular, but the file keywords are not; the singular
  form fails to parse.
- **Keywords are read in a fixed order.** A block pasted out of order fails.
  The relevant part of the order:

```
Timestep
TimeZeroPET                       << optional
TimeZeroCWD                       << optional
... maps (Species_CSV_File ... ClayMap) ...
LightningIgnitionsB0
LightningIgnitionsB1
AccidentalIgnitionsB0
AccidentalIgnitionsB1
IgnitionDistribution  ZeroInflatedPoisson   << optional; REQUIRED before the next 4
LightningIgnitionsBinomialB0
LightningIgnitionsBinomialB1
AccidentalIgnitionsBinomialB0
AccidentalIgnitionsBinomialB1
MaximumFineFuels
... Rx parameters (MaximumRxWindSpeed ... TargetRxSize, RxZonesMap) ...
MaximumSpreadAreaB0 .. B2
SpreadProbabilityB0 .. B3
SiteMortalityB0 .. B6             << keywords are singular "SiteMortality"
CohortMortalityB0 .. B2
LadderFuelMaxAge                  << integer
LadderFuelSpeciesList             << alone on its line; species on the NEXT line(s)
SuppressionMaxWindSpeed
Suppression_CSV_File
...
```

Every snippet the app exports is labelled with the section of the file it goes
in. Each snippet has been checked against this order by pasting it into SCF's
own v4 test parameter file.

---

## 2. Fire weather (FWI)

- **SCF consumes the LANDIS Climate Library's FWI.**
  - That FWI is computed from daily **mean** temperature
    (0.5 × (MinTemp + MaxTemp)) and **mean** relative humidity, not
    afternoon max temperature / min RH.
  - Calibrate on the same convention (`FWI convention: landis` in the
    sidebar). Otherwise every coefficient is fitted to a systematically higher
    FWI than SCF will feed it.
- **Known Climate Library issues (as of mid-2026):**
  - The FFMC rain phase is missing.
  - Precipitation in cm is fed into equations written for mm.
  - The DC day-length sign is wrong for Jan–Mar.

  These were being fixed upstream; check which library version your runs use
  and set **LANDIS Climate Library version** to match. The app cannot
  replicate the missing FFMC rain phase, so "pre-patch" FWI here still runs
  somewhat lower than an unpatched library's.
- **Check runtime against calibration FWI.** Load your LANDIS climate file
  under **Runtime FWI check** in the sidebar to compare the FWI SCF will see
  with the FWI you calibrated on.
- **Ignitions** use the **landscape-average** FWI for the day
  (PlugIn.cs:238). **Spread** uses the FWI of the ecoregion the fire is
  burning in (FireEvent.cs:177).

---

## 3. Ignitions

- **Daily ignition count per type:** a Poisson draw with
  λ = exp(b0 + b1 × FWI).
- **Zero-inflated Poisson (ZIP):**
  - SCF first computes zerosprob = plogis(bz0 + bz1 × FWI), the probability of
    a structural zero (PlugIn.cs:818).
  - It draws from the Poisson only when a uniform draw is ≥ zerosprob.
  - This matches `pscl::zeroinfl`, so E[n] = (1 − π) × λ transfers directly.
- **Area offset.** Coefficients are fitted over the calibration boundary but
  SCF applies them to the landscape.
  - Landscape b0 = regional b0 + log(landscape area / boundary area).
  - b1, bz0 and bz1 transfer unchanged.
  - The Tuning export applies this offset automatically. Without the areas
    loaded, it writes the ignition block commented out.

---

## 4. Spread

### Spread probability (FireEvent.cs `CanSpread`)

```
P = 1 / (1 + exp(-(B0 + B1*FWI + B2*FineFuels + B3*EWS)))
P = P * 0.71           for diagonal neighbours (FireEvent.cs:538)
P = P * suppressEffect from the suppression table
```

- **FineFuels** = min(SiteVars.FineFuels / `MaximumFineFuels`, 1)
  (FireEvent.cs:462).
  - The divisor is the **parameter**, not the raster maximum.
  - The calibration scales LANDFIRE FBFM40 so 1.0 = 9 t/ac (GR9), about
    **900 g C m⁻²**.
  - `MaximumFineFuels` must denote that load, or B2 (and SiteMortalityB5) is
    on the wrong scale.
  - If you raise it to avoid saturation, multiply B2 and B5 by
    (new value / 900).

### Spread is single-trial percolation

- **One trial per cell.** `CanSpread` marks a cell `Disturbed` **whether or
  not it ignites** (FireEvent.cs:456), and neighbour lookups skip `Disturbed`
  cells (:679, :707, :735). The ignition cell is marked at :105.
- **Consequence:** fire only keeps spreading when P is above the percolation
  threshold of this lattice, roughly 0.45–0.5 with the 0.71 diagonal weight.
  - A P of 0.05–0.2 at typical fire weather extinguishes nearly every
    ignition, whatever MaximumSpreadArea is.
  - Fire size jumps by orders of magnitude over a narrow range of B0.
- **B0 is not identifiable from the regression.** Case-control sampling biases
  the intercept, and the perimeter pairs measure a daily spread rate, not
  SCF's one-shot trial. **Set B0 by simulation**, using Parameter Tuning >
  **B0 Sweep**, then confirm with LANDIS runs on the Score Candidates tab.

### Maximum daily spread area (FireEvent.cs `Spread`)

```
MaxSpreadArea = (int)(B0 + B1*FWI + B2*EWS)        hectares
```

- **It is pacing, not a ceiling.**
  - The fire is a FIFO queue. Each burned cell adds one cell's area to the
    day's total.
  - When the total exceeds MaxSpreadArea, the day advances and the fire **keeps
    burning** (FireEvent.cs:204–216).
  - Each day therefore burns floor(MSA / cell area) + 1 cells. Fire size is set
    by when the queue empties.
- **The year-end cutoff.** A fire stops at day 364 (`DaysPerYear`,
  PlugIn.cs:47). So with easy spread, fire size can level off at about
  (cells per day) × (days left), and MaxSpreadArea then does limit size.
- **MSA ≤ 0** (B0 negative at low FWI):
  - SCF logs a warning and burns one cell per day.
  - Fires ignited on low-FWI days then creep until a hotter day lets them run.
  - Keep B0 > 0. The spread tab can constrain it to ≥ 50 ha.

### Suppression double-counting

- SCF multiplies P by the suppression effect from `Suppression_CSV_File`.
- The historical perimeters used for calibration were **already suppressed**.
- A non-zero suppression table therefore suppresses twice. Leave it at 0, or
  re-tune B0 upward.

### BurningSitesThreshold

Earlier notes described a `BurningSitesThreshold` option that forces spread
into any cell with N burning neighbours. **It is not in SCF master (4b089af)**
or in any copy we have found. If your SCF build has it, set it to 0 or re-tune
B0, because it bypasses P(spread) entirely.

---

## 5. Effective wind speed (Nelson 2002, Eq. 5)

```
EWS = Cb * sqrt((Ua/Cb)^2 + 2*(Ua/Cb)*sin(θ)*cos(φ) + sin(θ)^2)
```

- **θ** is the slope angle, and **φ** is the angle between the wind direction
  and the uphill direction. With terra's aspect (which is downhill), uphill
  azimuth = (aspect + 180) mod 360.
- **Cb (combustion buoyancy)** comes from the **source cell's** fire intensity
  (FireEvent.cs:572–578):
  - 10 by default
  - 25 if intensity > 3
  - 50 if intensity > 6

  Intensity = dNBR / 100, clamped to 1–10 (:347). Severity therefore feeds
  back into EWS, and EWS back into severity. The app calibrates at Cb = 10.
- **Wind units are km/h.** The Climate Library converts wind to km/h. ERA5 and
  GRIDMET wind (m/s) must be multiplied by 3.6 before calibration. Coefficients
  fitted on m/s are wrong by a factor of 3.6.

---

## 6. Mortality

### Site severity (FireEvent.cs:329)

```
eta  = B0 + B1*Clay + B2*PET + B3*EWS + B4*CWD + B5*FineFuels + B6*LadderFuels
dNBR = max(eta, 0.0005) ^ -1          (Gamma GLM, inverse link)
```

- **Signs are inverted.**
  - B0 must be **positive** (1 / B0 is the baseline dNBR).
  - Anything that **increases** severity has a **negative** coefficient.
- **The 0.0005 floor is a cliff.** Any cell whose eta reaches it gets
  dNBR = 2000. Treat any calibration pixels at the floor as a scale error.
- **Predictor units:**

| Term | Unit |
|------|------|
| Clay | fraction 0–1 (not percent) |
| PET | previous-year **potential** ET, mm/yr (TerraClimate `pet`, not `aet`) |
| EWS | day-of-burn effective wind, km/h |
| CWD | previous-year climatic water deficit, mm/yr |
| FineFuels | min(FineFuels / MaximumFineFuels, 1), 0–1 |
| LadderFuels | **raw summed cohort biomass, g m⁻²** (0 to ~3000), never a 0–1 index |

- **Ladder fuels** are the biomass of cohorts of `LadderFuelSpeciesList`
  species aged ≤ `LadderFuelMaxAge`. B6 must be calibrated on the same stratum.
  A 0–1 index gives a B6 roughly 1000 times too large, which drives eta
  negative and pins dNBR at 2000 everywhere.
- **`TimeZeroPET` / `TimeZeroCWD`** are used in year 0, before NECN has
  produced PET/CWD. Leaving them out sets those terms to 0 in the first year,
  which with negative B2/B4 **raises** first-year severity.

### Cohort mortality

```
P = exp(B0 + B1*BarkThickness + B2*SiteMortality);   p = P / (1 + P)
```

- **Expected signs:** B1 negative (thick bark survives), B2 positive.
- SiteMortality is SCF's own predicted dNBR (0–2000).
- **All or nothing:** each cohort gets one Bernoulli draw, and a kill removes
  the entire cohort's biomass.
- **Bark thickness is single-sided, in cm, everywhere.**
  - MaximumBarkThickness = FOFEM `BT_coef` × P95 DBH, with **no factor of 2**.
    Doubling it roughly halves simulated mortality.
  - FIA supplies no bark thickness. Cansler's FTM `BT_coef` is the FOFEM
    coefficient.

### Response data

- SCF severity is meant to track dNBR/RdNBR maps.
- Prefer **continuous per-fire MTBS rasters** (RdNBR, as in Robbins' reference
  calibration) over the 6-class thematic mosaic, which collapses the response
  to a few class midpoints.
- Keep `min_dnbr` at 1: higher truncation inflates B0 and flattens every slope.

Reference implementation: `Supporting R Code/Mortality/Calculating_Mortality.Rmd`
in the SCF repo (Robbins 2021). `FireDamages.cs` is dead code in v4.

---

## 7. Testing parameters without LANDIS

Parameter Tuning > **Fire Simulation** and **Fire Size Distribution** use
`simulate_fire_simple()`, a port of the SCF spread loop. It covers:

- the FIFO queue in SCF's neighbour order
- one trial per cell
- the 0.71 diagonal weight
- MaxSpreadArea pacing with C# `(int)` truncation
- the day-364 cutoff
- optional daily FWI and wind

It reproduces a literal per-cell translation of the C# loop cell for cell.

It does **not** model:

- suppression
- the Cb escalation of EWS
- spatially varying fuels or terrain
- the irregular landscape edge
- cells already disturbed by earlier fires in the same year

Treat its fire sizes as a fast screen and LANDIS runs as the test.

---

## 8. R gotcha (dplyr 1.1+)

Never create a dot-prefixed column inside `mutate()`, such as
`mutate(.d = as.Date(date))`. It breaks the data mask for bare column
references in the same call, so `date` resolves to `base::date` and errors.
