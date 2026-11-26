# spatial_synop_maps

Station-based (SYNOP) verification maps that compare two or more IFS
experiments with observations, built on the VINO tools:

1. **upload**: `~moz/bin/vino_genfc` retrieves the experiment's forecast
   fields from MARS, picks the nearest grid-point values at the observation
   stations and stores them in the VINO database.
   **collocate** (experiments with `collocate: local`) does the same with
   `vino_genfc`, but writes the station values to geopoint files instead of
   uploading them to VINO.
2. **download**: `~moz/bin/vino_getgeo` fetches the forecasts (VINO table
   `model`) and the verifying observations (table `observation`) as Metview
   geopoint files.
3. **plot**: pairs each forecast with the observations at its valid time and
   draws station maps of bias, RMSE and MAE, maps of mean bias (MEF), RMSE and
   MAE differences against a reference experiment, and score-vs-step curves.
   It then writes an `index.html` page for browsing all of them.

Each stage runs as a Slurm job, and each one can be rerun safely: anything
already done is skipped, so you can add dates or parameters to the config
and resubmit.

All geopoint files go to a **store** shared by every config (default
`./store` next to `synop_maps.py`, i.e. on `/perm`). Forecast files are filed
under what defines them (class, stream, expver, type, ...), not under the
run's experiment label. Any run that needs data another run already fetched
or computed finds it there and skips that work.

## Files

| File | Purpose |
|------|---------|
| `synop_maps.py` | Main script with the stages `submit`, `upload`, `download`, `plot` and `status` |
| `vino_upper.py` | Helper for upper-air (TEMP) parameters, run by `synop_maps.py` |
| `browser_template.html` | Template for the HTML browser. `plot` fills in the list of images and scores |
| `config_example.yaml` | Annotated example configuration |
| `README.md` | This file |

## Quick start

```bash
module load python3                      # Python >= 3.8 with yaml, pandas, matplotlib, cartopy
cd /perm/pavn/spatial_synop_maps
cp config_example.yaml my_exp.yaml        # edit dates, steps, parameters, experiments
python3 synop_maps.py submit my_exp.yaml  # submits upload -> download -> plot
squeue -u $USER                           # monitor
python3 synop_maps.py status my_exp.yaml  # how many files are available locally
```

When the plot job finishes, open `<workdir>/index.html` in a browser. The
default workdir is `./runs/<name>`. From your laptop you can copy the whole
`<workdir>` (or just `index.html` and `plots/`), or view it through any web
server you have. The page only uses relative links.

## Configuration

See `config_example.yaml` for all options. The main ones:

```yaml
name: my_comparison
dates: 20250601/to/20250630/by/1   # forecast base dates (MARS-style list/range)
times: [0, 12]                     # forecast base times
steps: 24/to/240/by/24
parameters:
  - param: 2t
  - {param: tp, period: 24}        # accumulations need a period; steps < period are dropped
  # add `seeps: true` to tp/period 24 for SEEPS maps and scores: Quaver's station
  # climatology (STVL, 1980-2009, valid month; fetched once into store/clim),
  # stations with 0.1 <= p1 < 0.85 only, values rounded to 0.1 mm as in Quaver
experiments:                       # FIRST experiment = reference for difference maps
  - {label: oper, class: od, stream: oper,  expver: "0001", type: fc, upload: false}
  - {label: exp1, class: rd, stream: lwda,  expver: abcd,   type: fc, anoffset: 9}
plot:
  regions: {global: null, europe: [-25, 45, 33, 72]}   # [west, east, south, north]
```

* `database: fdb` (per experiment) retrieves the forecasts from an alternative
  MARS database, e.g. experiments whose data is still in the FDB only. Needs
  explicit end dates (MARS cannot list such a database).
* `database: mars+fdb` (per experiment) takes each date from MARS and, where
  MARS does not have it (yet), from the FDB, e.g. an experiment whose older
  dates are archived and whose recent ones are still in the FDB only. Open end
  dates (`-`) use the dates found in either.
* `store` (top level) sets the shared data directory. Environment variables
  are expanded, e.g. `store: $PERM/synop_store`.
* `collocate: vino | local` (top level default, can be overridden per
  experiment and per period): `vino` uploads to and downloads from VINO,
  `local` computes the station values with `vino_genfc` and keeps them only in
  the store.
* `upload: false` is for data that is already in VINO. Operational
  forecasts (`class od, expver 0001`) are normally there already.
* `expver_label` (optional, per experiment) is passed to `vino_genfc
  --expver_label`. The data is then stored and read back from VINO under that
  label.
* `genfc_args` (per experiment) and `genfc.extra_args` (global) pass any other
  `vino_genfc` options, for example `--grid`, `--orog_from_file` or
  `--plugins`.
* `interpolation` is passed to both `vino_genfc` and `vino_getgeo`. The value
  must match: VINO stores the method together with the forecast values.

## Stages in detail

### Upload (`vino_genfc`)

The script runs one job per experiment that has `upload: true`. For every
parameter and base time:

1. **Check VINO.** `vino_getgeo` asks for the first and the last step of every
   date (with `--table model`). A date counts as uploaded only if both come
   back, so dates that are already in VINO are skipped, whoever uploaded
   them.
2. **Upload the rest.** The missing dates are sent to `vino_genfc
   --driven_by_forecast_date` in chunks of `chunk_size` dates. The defaults
   are `--ignore_missing` and `--get_orog_once`. `--get_lsm_once` is added for
   the land/sea interpolation methods.
3. **Check again.** After each chunk the script checks VINO again and logs
   any dates that are still missing (usually because the data isn't in MARS).
   A later rerun tries those dates again.

`python3 synop_maps.py upload cfg.yaml --force` re-uploads everything, for
example after an experiment was rerun under the same expver. The store copies
of re-uploaded dates are deleted, so `download` fetches them again.

If a `vino_genfc` call fails (a date missing in MARS, or a VINO database
clash between concurrent calls), the chunk is redone one date at a time and
dates whose call still fails are reported as missing. Their old data in VINO
(if any) is not trusted.

### Collocate (`vino_genfc`, no VINO upload)

For experiments with `collocate: local`, one job per experiment runs
`vino_genfc --store_in_database 0 --store_geopoint_files 1`: the same MARS
retrieval, orography and interpolation to the VINO observation stations as
the upload, but the result is written to geopoint files. The files are
written to a temporary directory inside the store and moved into place only
when the call succeeded, so the store never holds partial output. Dates
already in the store are skipped (`--force` recomputes them). Failed chunks
are redone one date at a time, as for the upload. These experiments are not
uploaded or downloaded; observations still come from VINO through
`download`.

### Download (`vino_getgeo`)

* Observations are fetched for every valid time (base time + step) and
  stored in `<store>/obs/<sources>/YYYY/MM/DD/<param>_obs_YYYYMMDDHH.geo`.
* Forecasts are stored in
  `<store>/fc/<vino|local>/<class>/<stream>/<expver>/<type>[/anoffset<n>...]/<interpolation>_<sources>/YYYY/MM/DD/<param>_mod_YYYYMMDDHH_<step>_<stream>_<expver>_<type>.geo`.
  `status` prints the directory of every experiment.
* A file that already exists in the store is not requested again.
  `vino_getgeo` writes to a temporary directory in the store and complete
  files are then moved into place.
* Runs of different configs can work on the store at the same time: each
  (directory, parameter, time) is locked while it is fetched or computed, and
  a second run waits and then only does what is still missing. Files that VINO
  doesn't have are simply not created, and the next run tries them again.
* At the end, the job prints how many of the expected files exist (the same
  output as `status`).

### Plot

* For each parameter, step and base date, the observation and forecast files
  are matched on station (id + lat/lon). Only station/forecast pairs that
  exist for **every** experiment are used, so all experiments are scored on
  the same sample.
* Per-station statistics over all base dates/times: bias (fc − obs), RMSE and
  MAE. MEF, RMSE and MAE differences are each experiment minus the reference
  (the first experiment). For RMSE/MAE differences, **red = worse than the
  reference**. As in quaver's
  `error_map`, a station is used only when it has valid pairs for at least
  `plot.min_fraction` (default 0.5) of the forecasts.
* Regional scores (scores.png, browser table) are weighted means of the
  station scores: bias and MAE are station means, RMSE is
  sqrt(mean of station RMSE²). By default (`plot.spatial_mean_weights:
  station_density`, as quaver's `spatial_mean_weights = "station_density"`)
  the stations are weighted by 1 / Σ exp(−(d/scangle)²) over their neighbours
  within 4·scangle (`plot.scangle`, default 0.75°), so each station of a dense
  network counts less; `none` weights all stations equally, as in `error_map`.
* Observations outside plausible ranges are rejected as in `error_map`
  (2t/tmax/tmin −100…60 °C, tp 0…900 mm, 10ff 0…50 m/s, tcc 0…100 %), and
  station elevations above 9000 m are treated as missing.
* Each colour scale is fixed across all steps of a parameter/region (98th
  percentile of the values), so maps for different steps can be compared.
* Wind direction (`10dd`) errors are wrapped into [−180°, 180°].
* Only pairs where the observation and every experiment's forecast are valid
  are used, so all maps and scores share the same stations.
* 2t forecasts are corrected from model orography to station height with
  `plot.lapse_rate` (default 0.0065 K/m, i.e. 6.5 K/km; 0 switches it off).
  Pairs with a missing model or station elevation are dropped for these
  parameters. tmax/tmin are corrected the same way.
* For all parameters, stations where the model and station elevations differ
  by more than `plot.max_dz` are excluded (stations without an elevation are
  kept, as in `error_map --orography_filter_threshold`); `max_dz: null`
  switches this off.
* Stations are classed by their identifier: `wmo` (numeric WMO block/station
  number), `icao` (4 letters, METAR), `noid` (no identifier; the ODB gives
  these reports `_` + a hash, e.g. dense national networks such as most of the
  Italian and Balkan stations, which often record much weaker winds than the
  WMO stations) and `other` (ships, mobile stations). `plot.station_ids`
  (e.g. `[wmo, icao]`) keeps only those classes in all statistics; the default
  `null` keeps all. `plot.station_groups` adds extra periods with only some
  classes, alongside the full set, e.g.
  `station_groups: {wmo: [wmo], noid: [noid]}` gives `JJA_wmo` and `JJA_noid`
  (and `JJA_00UTC_wmo` etc. with `valid_hours`).
* `orog_<experiment>.png` maps show the model minus station elevation per
  station (brown = model surface above the station), over the stations
  verified at any step for that parameter, region and period. They do not
  depend on the step.
* `rmserel` maps show the relative RMSE change of each experiment against the
  reference per station, `100 * (RMSE_exp / RMSE_ref - 1)` in %, so negative
  (blue) is better. The colour scale is capped at ±100 %.
* `scores.png` shows the bias/RMSE/MAE of all experiments against forecast
  step, averaged over all pairs in the region, plus the mean bias difference
  (MEF) and relative RMSE change against the reference.
* With `plot.daily_maps: true` the script also draws the observed values and
  each experiment's error for every single forecast. That produces
  dates × times × steps × (experiments + 1) images per region.
* With `plot.valid_hours: [0, 12]` every period gets extra periods (e.g.
  `DJF_00UTC`, `DJF_12UTC`, selectable in the browser) using only the
  forecasts valid at that hour (UTC), as an indication of the diurnal cycle.
  Each step then needs base times that reach those valid hours (e.g. `times:
  [0, 12]` with steps 24/to/240/by/48).
* The images used for each map are fingerprinted in `plots/signatures.json`,
  and an image is redrawn only when its data or colour scale changes.
  Appending dates therefore redraws the summary maps, and only the new daily
  maps. `--force` redraws everything.

### HTML browser (`<workdir>/index.html`)

* Selectors for parameter, region, view (Bias / MEF / RMSE / MAE / RMSE
  difference / MAE difference / Relative RMSE change / Model − station height /
  Scores vs step / Single forecasts) and forecast step.
  With ensembles there are also Fair CRPS / FCRPS difference / Relative FCRPS
  change views, and FCRPS columns in the table.
  Use the slider, the ◀ ▶ buttons or the ← → keys to change the step.
* All experiments are shown side by side. Below the maps, a table gives the
  bias, RMSE, MAE and relative RMSE change for the current step: green is
  better than the reference, red is worse.
* The current selection is kept in the URL (`#param=2t&view=rmse&...`), so you
  can bookmark or share a particular view.

## Slurm

`submit` writes job scripts to `<workdir>/jobs/` and submits them:

```
upload_<exp1>     ─┐
upload_<exp2>     ─┴─(afterany)─> download ─┐
collocate_<exp3>  ──────────────────────────┴─(afterany)─> plot
```

* Upload and collocate jobs for different experiments run in parallel.
* The dependencies use `afterany`, so a partly failed upload doesn't stop
  the rest of the chain. Whatever is available still gets downloaded and
  plotted, and a resubmission fills in the gaps.
* Logs go to `<workdir>/logs/<jobname>-<jobid>.out`.
* Resources are set under `slurm:` in the config (default `qos: nf`). The
  plot job uses `--cpus-per-task` worker processes.
* `--stages download plot` submits only some of the stages, and `--dry-run`
  writes the job scripts without submitting them.
* Each stage takes a lock file in `<workdir>/locks/`, so two copies of the
  same stage can't run on the same workdir at once. The store has its own
  locks in `<store>/locks/`.
* The VINO tools retrieve their MARS fields into a temporary directory per
  call under `$SCRATCHDIR`, which is deleted when the call ends. Everything
  that is kept goes to the store.

## Rerunning and appending

* **Resubmit unchanged.** Nothing is uploaded or downloaded again. Images
  whose data hasn't changed are not redrawn.
* **Extend `dates`** (for example to a longer period). Only the new dates
  are uploaded and downloaded, and the statistics are recomputed over the
  full period.
* **Add a parameter or an experiment.** Only that parameter or experiment is
  processed. A new experiment changes the common sample, so the maps are
  redrawn.
* The stages can also be run interactively without Slurm, for example
  `python3 synop_maps.py plot cfg.yaml`.

## Notes and limitations

* Without `control: true` an ensemble is only its perturbed members (`type: pf`).
* The units are whatever VINO stores (for example tp in mm and tcc in %).
  Forecasts and observations come from VINO in the same units.
* The upload check reads only the first and the last step. A date where a
  middle step failed to upload is counted as present. Use `--force` if you
  suspect this.
* `METVIEW_PYTHON_START_TIMEOUT` is set to 120 s for the VINO tools. The
  default of 8 s sometimes timed out during testing.

## Testing done

The pipeline was tested with operational data that is already in VINO:
`od/oper/0001/fc` against `od/enfo/0001/cf`, 2t and tp24,
20250101–20250104 00z, steps 24–72.

* `download` fetched all files, and a second run did nothing.
* Extending the dates downloaded only the new date.
* The upload check correctly told the dates in VINO apart from an expver
  that isn't there.
* `plot` drew all the maps, and a rerun skipped the unchanged images.
* The `submit` job chain ran under Slurm.

`vino_genfc` itself (the actual upload) has not been run yet: that needs one
of your experiments.

## Ensembles (fair CRPS)

An experiment with `type: pf` is an ensemble. Give its members with `numbers`
(passed to `vino_genfc --numbers`):

```yaml
  - {label: ctl, class: rd, stream: enfo, expver: j5me, type: pf, numbers: 1/to/8}
```

Both `collocate: vino` and `local` work: `vino_genfc` and `vino_getgeo` put all
members of one forecast in one geopoint file (columns `value_0 ... value_<n-1>`).
`numbers` is part of the store directory name, so different member sets don't
mix.

`control: true` adds the control forecast (`type: cf`, same class/stream/expver)
as one more member. It is uploaded/collocated and downloaded as a data set of its
own, and a forecast is used only when both its pf and its cf files are there.

* Per pair, the fair CRPS (Ferro 2014) is
  `mean_i |x_i - y| - sum_{i<j} |x_i - x_j| / (m (m - 1))` over the m members.
  Unlike the plain CRPS it does not favour ensembles with more members.
  Station FCRPS is the mean over the forecasts, and the regional score the mean
  over the stations.
* The maps `fcrps_<exp>`, `fcrpsdiff_<exp>` (minus the reference) and `fcrpsrel_<exp>`
  (`100 * (FCRPS_exp / FCRPS_ref - 1)`, i.e. minus the CRPSS against the reference)
  are drawn, and `scores.png` gets FCRPS panels.
  These views appear only when the config has at least one ensemble.
* Bias, RMSE and MAE of an ensemble are of the ensemble mean.
* A deterministic experiment can be compared with ensembles: its CRPS is its
  absolute error.
* The height correction is applied to every member. A station is used only
  when all of its members are valid.
* For `10dd` the differences are angular and the ensemble mean is the circular mean.

## Publishing on sites.ecmwf.int

Set `publish.site` (e.g. `https://sites.ecmwf.int/<user>/<site>/`) and put the
site's authentication token (Settings tab of the site) in `~/.spatial_maps_sites`
with `chmod 600`. The `plot` stage then uploads `index.html` and the plots to
`<site>/<run name>/` (change the folder with `publish.path`). Only new or
changed files are uploaded; `python3 synop_maps.py publish config.yaml [--force]`
uploads by hand.

## Several periods

Instead of `dates`, a config can define named periods:

```yaml
periods:
  DJF: 20231201/to/20240229/by/1
  JJA: 20240601/to/20240831/by/1
```

A period (or `dates`) can be open-ended: `20241201/to/-/by/1` runs up to the
last base date that MARS has for **every** experiment of that period. "Has"
means 2t at the longest configured step for every base time, and the date is
capped so that the observations for that step already exist. `submit` checks
MARS each time it is run (all experiments in parallel, ~30 s) and saves the
answer in `<workdir>/open_periods.json`. Its jobs, and any stage run by hand
afterwards, reuse that answer, so all stages work on the same dates; a stage
run without a saved answer checks MARS itself. Resubmitting later picks up
dates the experiments have reached since.

Upload and download cover all dates of all periods. Statistics, maps and
scores are computed per period and stored under `plots/<period>/`; the HTML
page gets a Period selector. Colour scales are shared between periods so the
maps can be compared directly. A config with plain `dates` is a single period
called `all`.

Experiments can use different settings in different periods but keep one
label (and so one line in the plots): give any setting as `{period: value}`,
e.g. `expver: {DJF: abcd, JJA: efgh}`. Every period must be listed.

## Observations from MARS feedback

By default the surface observations come from VINO (`sources`, default
`[synop]`). They can instead come from the feedback (ODB) of an assimilation in
MARS, using all of its stations:

```yaml
observations:
  source: mars                 # default: vino
  mars: {class: od, expver: 1, stream: lwda, type: ofb, obsgroup: conv}
```

* The observations are retrieved for every assimilation window the valid times
  need (one MARS request per analysis time and chunk of dates, about a minute
  per window), and written as one geopoint file per parameter and valid time to
  `<store>/obs/ofb_.../`. The directory name holds all settings, so different
  selections don't mix.
* Reports within `tolerance` minutes (default 30) of the valid hour are used,
  the closest one per station (active ones first if equally close). The window
  of analysis time T is (T + `window[0]`, T + `window[1]`] hours, by default
  that of the 12-hour 4D-Var of `lwda`. Change `window` and `window_times`
  for other streams.
* `obstype` (default 1, SYNOP incl. SHIP and METAR), `codetype` (e.g.
  `[11, 14]` for land SYNOP only) and `active_only` select the reports.
* The forecasts of every experiment (also `collocate: vino` ones) are picked at
  the stations of the observation file valid at each step by `vino_upper.py`
  in the `collocate` stage, which first retrieves the observations it needs.
  They are taken at the nearest grid point of the native grid (or the
  experiment's `grid`) and go to `<store>/fc/ofb/.../`. Nothing is uploaded to
  or downloaded from VINO for these parameters.
* The model height at the stations (for the height correction and `max_dz`) is
  z at step 0 of the forecast, or the file given with `--orog_from_file` in
  `genfc.extra_args` / `genfc_args`. The station height is the one in the
  feedback.
* Available are 2t, 2d, 10ff, 10dd, 10u and 10v (ODB varno 39, 40, 112, 111, 41,
  42) with `interpolation: nearest_point`. Accumulated parameters are not.
* Upper-air parameters (see below) then also come from the feedback: radiosonde
  reports (`upper_air.obstype`, default 5 = TEMP) at the standard pressure
  levels, ODB varno 2 (t), 1 (z), 112 (ff), 111 (dd) and 59 (td), in
  `<store>/obs/ofbtemp_.../`. Their time is the launch time (about an hour
  before the nominal hour), so the closest report within
  `upper_air.tolerance` minutes (default 90) is used. Surface and upper-air
  reports are retrieved in the same MARS request, and they are told apart by
  obstype, so the two obstype lists must not overlap:

  ```yaml
  observations:
    source: mars
    upper_air: {obstype: [5], codetype: null, tolerance: 90}
  ```
* Reading the feedback needs `pyodc` (in the `python3` module).
* A parameter can override the observation source with `obs: vino` (or
  `obs: mars`), e.g. tp, which is not in the feedback, from VINO SYNOP
  observations next to 2t from MARS:

  ```yaml
  parameters:
    - param: 2t
    - {param: tp, period: 24, obs: vino}
  ```

## Upper-air (TEMP) parameters

Radiosonde observations from the VINO `temp` data set (or, with
`observations.source: mars`, from the MARS feedback, see above) can be
verified like the SYNOP parameters. Give the pressure levels (hPa) of a parameter:

```yaml
parameters:
  - {param: t, levels: [850, 500, 250]}
  - {param: ff, levels: [250]}
upper_air:
  grid: O1280            # required: MARS grid for the forecast fields
  # interpolation: bilinear   # default nearest_point
  # sources: [temp]
```

Each level becomes a parameter of its own (`t850`, `t500`, ...) with the usual
maps, scores and browser entries. Available are `t`, `z`, `ff`, `dd` and `td`
(VINO has no other TEMP parameters).

* `vino_getgeo` cannot select a level and `vino_genfc` has no upper-air
  parameters, so both steps are done by `vino_upper.py` (run with the modules
  `python3 ecmwf-toolbox vtb/vino`):
  * observations come from VINO (`level_type=pl`) into `<store>/obs/temp/`;
  * forecasts are always collocated locally, for every experiment (also
    `collocate: vino` ones such as operations): the pressure-level fields are
    retrieved from MARS on `upper_air.grid` (or the experiment's `grid`) and
    picked at the stations of the observation file valid at each step. They go
    to `<store>/fc/upper/<class>/<stream>/<expver>/<type>[...]/grid<grid>_<interpolation>_temp/`.
* `ff` and `dd` are computed from u/v and `td` from t and q. `z` is plotted
  as geopotential height (gpm).
* This happens in the `collocate` stage, which first downloads the TEMP
  observations it needs. Steps whose valid time has no TEMP observations
  (e.g. 06/18 UTC at most stations) are skipped.
* There is no height correction and no `max_dz` filter (no model height at
  the stations), so the Model − station height view is not shown.
* The observations are used as stored in VINO, with a gross-error check only
  for t/td (170–340 K) and ff (≤ 150 m/s). A few bad soundings can inflate the
  RMSE, e.g. for z500.
