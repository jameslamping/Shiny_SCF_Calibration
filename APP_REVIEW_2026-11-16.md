# SCF Calibration App — Review, 2026-11-16

Reviewed at commit `5d0dae8` (now `main` on GitHub). Covers `app.R` and all
`modules/`. Every R file parses. The app was not launched end-to-end.

Findings are ranked by how badly they could produce wrong SCF parameters.

---

## Status at end of 2026-11-16 session

**Fixed and on `main`:**
- #1 candidate scoring
- #2 SCF-faithful simulator, plus a real-daily-weather ensemble
- #3 renv.lock, plus Windows/macOS fixes
- #5, #6, #7 parameter exports
- #11 single events-log reader
- most of #10

**New this session:**
- Parameter Tuning > **B0 Sweep** tab sets SpreadProbabilityB0 by simulation.
  This covers the B0 part of #9.
- Fixed the Observed Targets plot: its P75/P90 labels showed the P50 value.

**Still open:** #4, #8, the B1/B3 part of #9, #10 in the Validation tab, #12–#16,
and the MORA findings below.

### MORA findings from this session (data, not code)
- **`SiteMortalityB6` is NA in the saved MORA site-mortality fit.** It is marked
  fitted, but the estimate is missing. The export now writes it as `<SET THIS>`.
  Refit, or set it to 0.
- **`MaximumSpreadAreaB0` = −134 ha.** MSA is negative below about FWI 4, so SCF
  burns one cell per day there. Fires ignited on low-FWI days creep along until a
  hot day lets them run.
- **The saved fit still has `SpreadProbabilityB2` = −1.77.** You had already
  decided B2 = 0 for MORA (fuel load runs opposite to flammability in maritime
  forest), so set it to 0 after "Load from Calibration".
- **B0 sweep with the settled MORA choices** (B2 = 0, fire size only, since
  annual area doesn't transfer from region to park; 200 fires, runtime weather):
  - best B0 ≈ −1.56
  - median fire 5.7 ha vs 4.9 observed
  - **P90 only 17 ha vs 292 observed: the large-fire tail is far too thin**
  - B0 alone cannot match both typical and large fires; a steeper B1 (fires grow
    faster on hot days) is the obvious lever, hence the joint sweep in #9.
  - Set the sweep weights to size = 1, area = 0 for MORA. An earlier run
    (B2 = −1.77, 50/50 weights) gave B0 ≈ −0.12; disregard it.
- **CLAUDE.md is out of date on `BurningSitesThreshold`.** It describes this
  setting as part of FireEvent.cs, but it is absent from SCF master (4b089af) and
  from every local SCF copy.

---

## High — wrong answers or a broken setup for a new user

### 1. Score Candidates ranks nothing (spread_helpers.R:2661-2745) — FIXED 2026-11-16
> Rewritten: shared events-log reader, log-ratio size + fire-rotation scores
> on one scale, size truncation to [1 cell, landscape] (also covers #10),
> empty runs rank last. Snippet ZIP now writes each candidate's B1 too.
- `read_scf_events()` looks for `FireEventLog.csv` / `fire-event-log.csv` /
  `SCF-event-log.csv`, with a column named `sites_burned`. SCF v4 writes
  `socialclimatefire-events-log.csv` with `TotalSitesBurned`. The rest of the app
  already uses those names (`app.R:691`, `validation_helpers.R:79`). So every
  candidate reads as an empty event set.
- `score_annual_area_rmse()` returns `sqrt(mean(sim^2))`. That is the RMSE
  against **zero**, not against observed annual area burned, so it rewards
  candidates that burn less.
- The size score is on a log1p scale (about 0–5). The area score is in raw
  hectares (thousands), so the area score swamps the total even after weighting.
- Net effect: all candidates come out at about 0.5, or the ones that burn least
  rank first. **Don't trust this tab until it's fixed.**

### 2. The in-app fire simulator doesn't behave like SCF (tuning_helpers.R:52) — FIXED 2026-11-16
> Rewritten to FireEvent.cs rules (single trial, 0.71 diagonals, MSA pacing,
> SCF queue order, optional daily weather). Matches a literal per-cell port of
> the C# loop cell-for-cell. The Tuning ensemble now also runs on real daily
> weather (default mode).
`simulate_fire_simple()` drives the Tuning ensemble. It differs from FireEvent.cs
in three ways that bias exactly the parameters it's used to tune:
- **Retries:** a cell that fails its spread trial stays unburned, and a
  neighbour that ignites the next day can try it again. SCF marks it
  `Disturbed` after one trial. Near the percolation threshold (P ≈ 0.4) this
  overstates spread, so B0 tuned here will be too low for SCF.
- **MaxSpreadArea:** excess cells are **discarded** (`sample(new_idx, msa_cells)`).
  SCF carries them into the next day and the fire keeps burning. The simulator
  treats MSA as a ceiling on fire size, which CLAUDE.md says is the wrong
  reading.
- **Diagonals:** no 0.71 weight on diagonal spread.

### 3. renv.lock is missing packages the app needs — FIXED 2026-11-16
The lockfile lacks `cffdrs`, `quantreg`, `nnls`, `patchwork` and `tigris`. A
labmate running `renv::restore()` won't get them.
- Without `cffdrs`, GRIDMET FWI and runtime FWI stop with an error. These are
  the core of the recent FWI-alignment work.
- Without `quantreg`, `fit_max_daily_area()` **silently** falls back to nnls or
  clamped OLS. Two machines can then produce different MaxSpreadArea
  coefficients from the same data.
- Fix: `renv::install(c("cffdrs","quantreg","nnls","patchwork","tigris"))`, then
  `renv::snapshot()`.

### 4. Hardcoded personal paths as input defaults
About 13 `textInput` defaults in `app.R` point at `/Users/jlamping/...`
(lines 89, 130, 276–282, 336, 850, 876, 938, 2016, 2500). The same goes for
`spread_probability_diagnostic.R:63–246` and `analysis/*.R`.
- A new user gets file-not-found errors on first run.
- The defaults also mix parks: the template is NOCA, the NECN table is OLYM, and
  the diagnostic script is MORA.
- Fix: use empty defaults with placeholders, or a `config.yml` per user (gitignored).

---

## Medium — export and method gaps

### 5. The mortality snippet may not paste into SCF (mortality_helpers.R:1644) — FIXED 2026-11-16
> All three snippet writers now use LANDIS syntax and SCF's parse order (checked
> against an emulation of InputParameterParser.cs on SCF's own v4 test file).
> The spread snippet also had MaximumSpreadArea and SpreadProbability in the
> wrong order, and TimeZeroPET/CWD belong right after Timestep.
`write_mortality_snippet()` writes `# comments` and `SiteMortalityB0 = 0.003`.
The spread and tuning writers use LANDIS syntax: `>>` comments and
`Keyword  value` with no `=`. If you paste the mortality snippet directly, the
parser will most likely reject it. Make it match the other writers, then confirm
with one SCF run.

### 6. Tuning export loses precision and has a stale default (tuning_helpers.R:1274) — FIXED 2026-11-16
> Correction: the sidebar default is 18000, not 500 (500 was only the export
> fallback, now a <SET THIS> line). Default changed to 900 on 2026-11-16.
- `SiteMortalityB*` are printed with `%14.8f`. A B6 of about −7.9e-7 prints as
  `-0.00000079`, which keeps only 2 significant digits. B3 has the same problem.
  Use `%.6e` or `%.8g`.
- The `MaximumFineFuels` default is 500 (`app.R:2682`, `7640`). CLAUDE.md says
  about 900 is needed to match the FBFM40 9 t/ac normalisation that B2 and B5
  were fitted on.

### 7. The tuning export writes regional ignition b0 without the area offset — FIXED 2026-11-16
The ignition block in `build_scf_tuning_snippet()` writes the
calibration-boundary b0 and leaves only a comment saying to add
log(landscape/boundary). It looks ready to paste but is wrong by that offset.
Apply the offset automatically, or refuse to export until the Scale Adjustment
step has been done.

### 8. Runtime FWI copies only some of the Climate Library's behaviour (era_helpers.R:673)
- `compute_landis_climate_fwi()` copies mean temp/RH (Issue 4) and the cm/mm
  precip mismatch (Issue 2).
- It does **not** copy the missing FFMC rain phase (Issue 1; the code comments
  this but leaves it in) or the DC day-length sign error (Issue 3).
- Until the patched library ships and is the version you run, the
  "Calibration vs Runtime FWI" plot understates the gap.
- Option: add an `unpatched` mode that reimplements FFMC and DC the way the
  library does.
- Minor: `lat = 46.85` (MORA) is hardcoded and never passed from the app
  (`app.R:5750`). This doesn't matter in the PNW because cffdrs uses latitude
  bands, but it would elsewhere.

### 9. The spread-probability response is a daily rate, but SCF draws once per fire — B0 PART DONE 2026-11-16
> B0 now comes from the B0 Sweep tab. B1/B3 are still taken from the GLM; a
> joint sweep (B0 × B1) is the natural next step: for MORA no single B0 matches
> both the median and the P90 fire size.
- The perimeter-pair frontier ring counts a cell as a failure if it is unburned
  at t+1, even if it burns at t+2. SCF never gives that cell a second chance.
- This affects more than the intercept: it also distorts the B1 and B3 slopes,
  because daily and per-event probabilities aren't proportional.
- CLAUDE.md already says B0 has to come from a sweep. Once #2 is fixed, consider
  fitting B0, B1 and B3 together against observed sizes and growth in the
  simulator (a simulation-based fit) instead of reading the slopes off the GLM.

### 10. The fire-size comparison mixes scales — MOSTLY FIXED 2026-11-16
> Candidate scoring, the real-weather ensemble and the B0 sweep truncate both
> sides to (1 cell, landscape/grid]. The Validation tab size plots still apply
> only the lower cutoff (`plot_val_size_ecdf`, `plot_val_fire_size`).
`score_candidates()` and the validation size plots compare FPA FOD sizes from
the whole calibration region with fires simulated on the landscape.
- The observed set includes many fires smaller than one cell, and some larger
  than the whole landscape. SCF can produce neither.
- Truncate both sides to [1 cell, landscape area] before comparing quantiles.

---

## Low — maintenance and handoff

11. *(Fixed 2026-11-16: all three now go through `load_scf_events_full()`.)*
    **Three separate SCF event-log parsers.** They are
    `ignition_helpers.R:905`, `spread_helpers.R:2716` and
    `validation_helpers.R:59`, and they assume different filenames and columns.
    Merging them into one parser fixes #1 at the root.
12. **Stale wind-unit comments.** `era_helpers.R:78–82` and `:193` say wind is
    returned in m/s, but the code returns km/h.
13. **README out of date.** It was last updated 2026-04-29. It has no
    Mortality Calibration tab, GRIDMET/cffdrs FWI convention, FIRED/NIROPS
    perimeters, percolation diagnostics, `MaximumFineFuels` or
    `BurningSitesThreshold`.
14. **CLAUDE.md is gitignored.** It holds the most important SCF-semantics notes:
    plural keywords, single-trial percolation, suppression double-counting, wind
    units and mortality sign inversion. Your labmate never sees them. Move them
    into a tracked `docs/SCF_semantics.md`.
15. **No tests, and `app.R` is 8,800 lines.** A small `testthat` suite would catch
    regressions like #1 and #5. Candidate tests:
    - EWS formula
    - snippet round-trip parse
    - mortality coefficient recovery on synthetic data (this was already done
      once, by hand)
    - candidate scoring on a fixture log

    Splitting the tabs into Shiny modules would make review possible.
16. **Open items from the August mortality review:**
    - switch the response to continuous per-fire RdNBR
    - check the PET/CWD scale against one LANDIS run

---

## Suggested order (remaining)
1. **#14 + #13 handoff docs:** your labmate is using the app now.
   - Move the SCF semantics out of the gitignored CLAUDE.md into a tracked
     `docs/SCF_semantics.md`, and fix the BurningSitesThreshold note.
   - Bring the README up to date: Mortality tab, real-weather ensemble, B0 Sweep,
     scoring, exports.
2. **#4 default paths:** a gitignored per-user config instead of
   `/Users/jlamping/...` defaults.
3. **#9 joint B0 × B1 sweep:** fixes the shape of the size distribution, not just
   its level.
4. **MORA data fixes:** refit B6 (or set 0), positive MSA B0.
5. **#15 tests:** turn this session's scratch checks into a `testthat` suite:
   - simulator vs SCF reference loop
   - scoring known-answer
   - snippet parser emulation
6. #8, #10 (Validation tab), #12, #16
