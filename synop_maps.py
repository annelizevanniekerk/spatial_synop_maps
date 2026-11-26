#!/usr/bin/env python3
"""
synop_maps.py - collocate experiment forecasts with SYNOP observations from
VINO, and plot station-based verification maps that can be browsed in a
single HTML page.

Stages (each can be rerun safely; work already done is skipped):

    upload    vino_genfc   : MARS forecast fields -> station values in VINO
    collocate vino_genfc   : MARS forecast fields -> station values in local
                             geopoint files (no VINO upload; collocate: local)
              vino_upper.py: same for pressure-level parameters at TEMP stations
                             (always local, for every experiment)
    download  vino_getgeo  : VINO -> local geopoint files (obs + forecasts);
                             TEMP observations with vino_upper.py;
                             with observations.source: mars the SYNOP
                             observations come from MARS feedback (type ofb)
                             and the forecasts are picked at its stations by
                             vino_upper.py in the collocate stage
    plot      matplotlib   : station maps + scores + index.html browser

Geopoint files are kept in a store directory shared by all configs (`store`,
default ./store next to this script), so data computed or downloaded for one
run is reused by every other run that needs it.

Usage:
    python3 synop_maps.py submit    config.yaml   # submit all stages to slurm
    python3 synop_maps.py upload    config.yaml [--experiment LABEL]
    python3 synop_maps.py collocate config.yaml [--experiment LABEL]
    python3 synop_maps.py download  config.yaml
    python3 synop_maps.py plot     config.yaml   # also publishes if publish.site is set
    python3 synop_maps.py publish  config.yaml [--force]
    python3 synop_maps.py status   config.yaml

Run `module load python3` first (needs Python >= 3.8 with yaml, pandas,
matplotlib and cartopy).
"""

import sys

if sys.version_info < (3, 8):
    sys.exit("synop_maps.py needs Python >= 3.8: run 'module load python3' first")

import argparse
import datetime as dt
import fcntl
import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

SCRIPT = Path(__file__).resolve()
MISSING = 1e30  # geopoint missing value is 3e38
# parameters corrected for the model/station height difference (as quaver error_map)
LAPSE_PARAMS = {"2t", "tmax", "tmin", "mx2t", "mn2t"}
# gross-error limits for observations (quaver error_map "toss"); temperatures in K
OBS_RANGE = {"2t": (173.15, 333.15), "tmax": (173.15, 333.15), "tmin": (173.15, 333.15),
             "mx2t": (173.15, 333.15), "mn2t": (173.15, 333.15),
             "tp": (0, 900), "10ff": (0, 50), "tcc": (0, 100),
             "t": (170, 340), "td": (170, 340), "ff": (0, 150)}
ANGLE_PARAMS = {"10dd", "dd"}  # wind directions: angular errors and circular means
# pressure-level parameters available from the VINO TEMP data set
UPPER_PARAMS = {"t", "z", "ff", "dd", "td"}
G = 9.80665  # z is plotted as geopotential height (gpm)
# ODB varno of the surface parameters available from MARS feedback (observations.source: mars)
OFB_VARNO = {"2t": 39, "2d": 40, "10u": 41, "10v": 42, "10dd": 111, "10ff": 112}
# ODB varno of the pressure-level parameters (TEMP) in MARS feedback; z is geopotential
OFB_UPPER_VARNO = {"z": 1, "t": 2, "td": 59, "dd": 111, "ff": 112}
UPPER_ENV = "module load python3 ecmwf-toolbox vtb/vino >/dev/null 2>&1"
# MARS read buffer: one tco3999 GRIB field is ~128 MB, above the 64 MB default.
# Set here so every MARS call made by this script and its sub-tools inherits it.
os.environ.setdefault("MARS_READANY_BUFFER_SIZE", str(1024 * 1024 * 1024))

DEFAULTS = {
    "name": "synop_maps",
    "workdir": None,
    "store": None,  # shared geopoint store; default: ./store next to this script
    "collocate": "vino",  # vino: upload to + download from VINO; local: genfc -> files
    "times": [0],
    "sources": ["synop"],
    "observations": {
        "source": "vino",  # vino | mars: surface observations from VINO (sources) or MARS feedback (ODB)
        # MARS feedback request (source: mars); time and date are added per window
        "mars": {"class": "od", "expver": 1, "stream": "lwda", "type": "ofb", "obsgroup": "conv"},
        "window": [-3, 9],  # assimilation window of each `window_times` analysis, in hours
        "window_times": [0, 12],
        "tolerance": 30,  # minutes; reports within this of the valid hour, closest one used
        "obstype": [1],  # ODB obstype (1 = SYNOP, incl. SHIP and METAR)
        "codetype": None,  # e.g. [11, 14] for land SYNOP only; null = all
        "active_only": False,  # only observations active in the assimilation
        # pressure-level parameters (source: mars): radiosonde reports at the standard levels.
        # Their time is the launch time (about 1 h before the nominal hour), so the
        # report closest to a valid time within `tolerance` minutes is used.
        "upper_air": {"obstype": [5], "codetype": None, "tolerance": 90},
    },
    "interpolation": "nearest_point",
    "chunk_size": 15,
    "parallel": 4,  # concurrent vino_genfc / vino_getgeo calls
    "tools": {
        "genfc": "~moz/bin/vino_genfc",
        "getgeo": "~moz/bin/vino_getgeo",
        "mars": "mars",
        "fdb": "fdb",  # lists experiments with `database: fdb` (MARS cannot)
    },
    "genfc": {
        "ignore_missing": True,
        "get_orog_once": True,
        "extra_args": [],
    },
    "upper_air": {
        "sources": ["temp"],  # VINO data set of the pressure-level observations
        "grid": None,  # MARS grid of the forecast fields, e.g. O1280 (an experiment's grid wins)
        "interpolation": "nearest_point",  # nearest_point | bilinear
    },
    "plot": {
        "regions": {"global": None},
        "min_fraction": 0.5,  # station used if valid for >= this fraction of forecasts (error_map minpopul)
        "daily_maps": False,
        "processes": 4,
        "dpi": 150,
        "lapse_rate": 0.0065,  # K/m, height correction of 2t forecasts to station height
        "max_dz": 250,  # m, drop stations where |model - station elevation| exceeds this
        # weights of the stations in the regional scores (quaver spatial_mean_weights):
        # none (all equal) | station_density (1 / Gaussian-kernel station density, so
        # each station of a dense network counts less), kernel width scangle in degrees
        "spatial_mean_weights": "station_density",
        "scangle": 0.75,
        "valid_hours": [],  # e.g. [0, 12]: extra periods with only these valid hours (UTC)
        # station identifier classes (see station_class): wmo (numeric WMO block/station
        # number), icao (4 letters, METAR), noid (no identifier: '_' + hash in the ODB,
        # e.g. national networks), other (ships, mobile, 7-character call signs)
        "station_ids": None,  # classes to keep everywhere, e.g. [wmo, icao]; null = all
        "station_groups": {},  # extra periods per group, e.g. {wmo: [wmo], noid: [noid]}
    },
    "publish": {
        "site": None,  # e.g. https://sites.ecmwf.int/<user>/<site>/ ; null = do not publish
        "path": None,  # folder on the site, default: the run name
        "token_file": "~/.spatial_maps_sites",
        "parallel": 8,
    },
    "slurm": {
        "qos": "nf",
        "account": None,
        "upload": {"time": "24:00:00", "mem": "16G", "cpus": 4},
        "collocate": {"time": "24:00:00", "mem": "16G", "cpus": 4},
        "download": {"time": "12:00:00", "mem": "16G", "cpus": 4},
        "plot": {"time": "06:00:00", "mem": "32G", "cpus": 4},
    },
}

EXP_DEFAULTS = {
    "class": "rd",
    "stream": "oper",
    "type": "fc",
    "numbers": None,  # ensemble members for type pf, e.g. 1/to/50
    "control": False,  # type pf: add the control forecast (type cf) as a member
    "upload": True,
    "collocate": None,  # default: the global `collocate`
    "anoffset": None,
    "origin": None,
    "expver_label": None,
    "grid": None,
    "resol": None,
    "database": None,  # alternative MARS database, e.g. fdb for data not archived yet;
    # mars+fdb: each date from MARS, and from the FDB where MARS does not have it
    "genfc_args": [],
}


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def log(*args):
    stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{stamp}]", *args, flush=True)


def deep_update(base, new):
    out = dict(base)
    for k, v in (new or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_update(out[k], v)
        else:
            out[k] = v
    return out


def tokens(spec):
    """Split a MARS-style spec (str, int or list) into string tokens."""
    if spec is None:
        return []
    if isinstance(spec, (list, tuple)):
        out = []
        for s in spec:
            out += tokens(s)
        return out
    return [t for t in str(spec).replace("/", " ").split() if t]


def expand(spec, parse, step, fmt):
    """Expand 'a/to/b/by/n' style specs; parse/step/fmt define the value type."""
    tok = [t.lower() for t in tokens(spec)]
    out, i = [], 0
    while i < len(tok):
        if i + 2 < len(tok) and tok[i + 1] == "to":
            a, b = parse(tok[i]), parse(tok[i + 2])
            n = 1
            i += 3
            if i + 1 < len(tok) and tok[i] == "by":
                n = int(tok[i + 1])
                i += 2
            v = a
            while v <= b:
                out.append(v)
                v = step(v, n)
        else:
            out.append(parse(tok[i]))
            i += 1
    return sorted(set(out), key=fmt)


def expand_dates(spec):
    return expand(spec,
                  lambda s: dt.datetime.strptime(s, "%Y%m%d").date(),
                  lambda d, n: d + dt.timedelta(days=n),
                  lambda d: d)


def expand_ints(spec):
    return expand(spec, int, lambda v, n: v + n, lambda v: v)


def ymd(d):
    return d.strftime("%Y%m%d")


def chunks(seq, n):
    seq = list(seq)
    for i in range(0, len(seq), max(1, n)):
        yield seq[i:i + n]


@contextmanager
def lock(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        try:
            fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            sys.exit(f"another run holds {path}; is a job for this stage already running?")
        yield


@contextmanager
def store_lock(path):
    """Wait for exclusive use of one part of the shared store, so that runs
    of different configs never fetch or compute the same files at once."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        yield


def move_into(src, dst):
    """Move all .geo files below src to the same relative place below dst.
    The store only ever sees complete files (os.replace is atomic)."""
    n = 0
    for f in Path(src).rglob("*.geo"):
        target = Path(dst) / f.relative_to(src)
        target.parent.mkdir(parents=True, exist_ok=True)
        os.replace(f, target)
        n += 1
    return n


def run(cmd):
    log("RUN", " ".join(str(c) for c in cmd))
    env = dict(os.environ)
    # metview (used inside the vino tools) can be slow to start on busy nodes
    env.setdefault("METVIEW_PYTHON_START_TIMEOUT", "120")
    # the vino tools retrieve into $SCRATCHDIR with fixed file names
    # (mars_request_N.grib, tmpout), so concurrent calls need their own
    # directory or they overwrite each other's data
    base = os.environ.get("SCRATCHDIR") or os.environ.get("TMPDIR") or None
    if base:
        os.makedirs(base, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="vino_", dir=base) as tmp:
        env["SCRATCHDIR"] = env["TMPDIR"] = tmp
        env.pop("VINO_TMPDIR", None)
        return subprocess.run([str(c) for c in cmd], env=env).returncode


# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------

OPEN_END = "-"
MAX_GAP = 7  # days; an open period ends at the first longer gap in the MARS dates


MARS_FDB = "mars+fdb"  # `database` value: MARS first, the FDB for the dates MARS lacks


def db_attempts(spec):
    """The specs to retrieve an experiment with, in order: with
    `database: mars+fdb` first from MARS, then from the FDB."""
    if spec.get("database") == MARS_FDB:
        return [dict(spec, database=None), dict(spec, database="fdb")]
    return [spec]


def db_where(e):
    db = e.get("database")
    return "MARS" if not db else "MARS+FDB" if db == MARS_FDB else "FDB" if db == "fdb" else f"database {db}"


def period_spec(e, per):
    """Experiment settings for one period (values given as {period: value} resolved)."""
    return {k: ({str(a): b for a, b in v.items()}.get(per) if isinstance(v, dict) else v)
            for k, v in e.items()}


def list_dates(cfg, e, t, step, start, cap):
    """Base dates in [start, cap] at time t for which experiment e has 2t at step.
    None if the listing failed."""
    if e.get("database") == MARS_FDB:
        found = [list_dates(cfg, sp, t, step, start, cap) for sp in db_attempts(e)]
        if all(f is None for f in found):
            return None
        return sorted({d for f in found if f for d in f})
    req = {"class": e["class"], "stream": e["stream"], "expver": e["expver"],
           "type": e["type"], "levtype": "sfc", "param": "167", "time": f"{t:02d}00",
           "step": step}
    for k in ("anoffset", "origin"):
        if e.get(k) not in (None, ""):
            req[k] = e[k]
    if e.get("database") == "fdb":
        # MARS cannot LIST the FDB; fdb list takes no date ranges, so filter here
        text = ",".join(f"{k}={v}" for k, v in req.items())
        log(f"fdb list: {text}")
        r = subprocess.run([cfg.tools["fdb"], "list", text], capture_output=True, text=True)
        pattern = r"\bdate=(\d{4})(\d{2})(\d{2})\b"
    elif e.get("database"):
        sys.exit(f"experiment {e['expver']}: open end dates ('-') are not supported with "
                 f"`database: {e['database']}`, give an end date")
    else:
        req.update(date=f"{ymd(start)}/to/{ymd(cap)}", output="tree")
        text = "list," + ",".join(f"{k}={v}" for k, v in req.items()) + "\n"
        log(f"MARS list: {text.strip()}")
        r = subprocess.run([cfg.tools["mars"]], input=text, capture_output=True, text=True)
        pattern = r"\b(\d{4})-(\d{2})-(\d{2})\b"
    if r.returncode:
        log("listing failed:", (r.stderr or r.stdout).strip()[-500:])
        return None
    return sorted({d for d in (dt.date(int(y), int(m), int(dd))
                               for y, m, dd in re.findall(pattern, r.stdout))
                   if start <= d <= cap})


def mars_last_date(cfg, e, start):
    """Last base date >= start for which MARS (or the FDB) has the longest
    configured step of experiment e at every configured base time (2t is used
    as the probe), and whose observations for that step already exist (valid
    time before yesterday). None if there is nothing."""
    maxstep = max(cfg.steps)
    cap = dt.date.today() - dt.timedelta(days=maxstep // 24 + 2)
    if cap < start:
        return None
    last = []
    for t in cfg.times:
        dates = list_dates(cfg, e, t, maxstep, start, cap)
        # the data must start near the period start and is followed only up to
        # the first gap of more than MAX_GAP days, so that e.g. a run of the
        # same expver in another year does not extend the period
        if not dates or (dates[0] - start).days > MAX_GAP:
            where = db_where(e)
            log(f"{e.get('label', e['expver'])} ({e['expver']}) {t:02d}z: no forecasts in {where} from "
                f"{ymd(start)}" + (f" (first date found {ymd(dates[0])})" if dates else ""))
            return None
        end = dates[0]
        for d in dates[1:]:
            if (d - end).days > MAX_GAP:
                break
            end = d
        last.append(end)
    return min(last)


class Param:
    def __init__(self, d, level=None):
        if isinstance(d, str):
            d = {"param": d}
        self.param = str(d["param"])
        self.period = int(d.get("period", 0) or 0)
        self.level = level  # pressure level in hPa (TEMP), None for surface parameters
        self.local = self.upper  # collocated by vino_upper.py at the stations of the obs files
        # observation source of this parameter (vino | mars); None: observations.source
        self.obs = d.get("obs")
        if self.obs not in (None, "vino", "mars"):
            sys.exit(f"{self.param}: obs must be vino or mars")
        self.mars = False  # observations from MARS feedback (set by Config)
        # SEEPS (Rodwell et al. 2010) with the station climatology Quaver uses
        self.seeps = bool(d.get("seeps", False))
        if self.seeps and (self.param != "tp" or self.period != 24 or level is not None):
            sys.exit(f"{self.param}: seeps is available for tp with period: 24 only")
        if level is not None and self.param not in UPPER_PARAMS:
            sys.exit(f"{self.param}: pressure levels are available for {sorted(UPPER_PARAMS)} only")
        default = f"{self.param} {level} hPa" if level is not None else self.key
        self.label = str(d.get("label", default)).format(level=level)

    @classmethod
    def expand(cls, d):
        """One Param per level of {param: t, levels: [850, 500]}."""
        if isinstance(d, dict) and d.get("levels") is not None:
            return [cls(d, int(lev)) for lev in expand_ints(d["levels"])]
        return [cls(d)]

    @property
    def upper(self):
        return self.level is not None

    @property
    def key(self):
        """Name used in file names, e.g. 2t, tp24, t850."""
        if self.upper:
            return f"{self.param}{self.level}"
        return f"{self.param}{self.period}" if self.period else self.param


class Config:
    def __init__(self, path, refresh_open=False):
        self.path = Path(path).resolve()
        with open(self.path) as fh:
            raw = yaml.safe_load(fh)
        c = deep_update(DEFAULTS, raw)
        self.raw = c
        self.name = c["name"]
        self.workdir = Path(c["workdir"] or self.path.parent / "runs" / self.name).expanduser().resolve()
        self.store = Path(os.path.expandvars(str(c["store"] or SCRIPT.parent / "store"))).expanduser().resolve()
        self.times = expand_ints(c["times"])
        self.steps = expand_ints(c["steps"])
        self.params = [q for p in c["parameters"] for q in Param.expand(p)]
        if len({p.key for p in self.params}) != len(self.params):
            sys.exit("a parameter is configured twice")
        self.upper = c["upper_air"]
        self.upper["sources"] = tokens(self.upper["sources"])
        if self.upper["interpolation"] not in ("nearest_point", "bilinear"):
            sys.exit("upper_air.interpolation must be nearest_point or bilinear")
        self.tools = {k: os.path.expanduser(v) for k, v in c["tools"].items()}
        # named verification periods (e.g. DJF, JJA); plain `dates` is one period "all".
        # An end date "-" (e.g. 20241201/to/-/by/1) means: up to the last date in MARS.
        raw_periods = ({str(k): v for k, v in c["periods"].items()} if c.get("periods")
                       else {"all": c["dates"]})
        ends = self.resolve_open_ends(raw_periods, c["experiments"], refresh_open)
        self.periods = {}
        for k, v in raw_periods.items():
            ds = expand_dates([ymd(ends[k]) if tok == OPEN_END else tok for tok in tokens(v)]
                              if k in ends else v)
            if not ds:
                sys.exit(f"period {k}: no dates ({v}, last date in MARS {ymd(ends[k])})")
            self.periods[k] = ds
        self.dates = sorted({d for ds in self.periods.values() for d in ds})
        self.sources = tokens(c["sources"])
        self.obs = c["observations"]
        if self.obs["source"] not in ("vino", "mars"):
            sys.exit("observations.source must be vino or mars")
        self.interpolation = c["interpolation"]
        for p in self.params:
            p.mars = (p.obs or self.obs["source"]) == "mars"
        if self.obs_mars:
            if self.interpolation != "nearest_point":
                sys.exit("observations.source: mars needs interpolation: nearest_point")
            self.obs["upper_air"]["obstype"] = tokens(self.obs["upper_air"]["obstype"])
            # reports are told apart by obstype (ff/dd share their varno with 10ff/10dd)
            if (any(p.upper for p in self.mars_params) and any(not p.upper for p in self.mars_params)
                    and set(tokens(self.obs["obstype"])) & set(self.obs["upper_air"]["obstype"])):
                sys.exit("observations.obstype and observations.upper_air.obstype must not overlap")
            for p in self.mars_params:
                if p.upper and p.param not in OFB_UPPER_VARNO:
                    sys.exit(f"{p.key}: not available with observations.source: mars "
                             f"(pressure levels: {sorted(OFB_UPPER_VARNO)})")
                if not p.upper and (p.param not in OFB_VARNO or p.period):
                    sys.exit(f"{p.key}: not available with observations.source: mars "
                             f"(available: {sorted(OFB_VARNO)})")
                p.local = True
        self.chunk_size = int(c["chunk_size"])
        self.parallel = max(1, int(c["parallel"]))
        self.genfc = c["genfc"]
        self.plot = c["plot"]
        self.slurm = c["slurm"]
        self.publish = c["publish"]
        self.experiments = []
        for e in c["experiments"]:
            e = deep_update(EXP_DEFAULTS, e)
            if e["collocate"] is None:
                e["collocate"] = c["collocate"]
            if "label" not in e:
                if isinstance(e["expver"], dict):
                    sys.exit("experiments with a per-period expver need a label")
                e["label"] = e["expver"]
            e["label"] = str(e["label"])
            e["variants"] = self.exp_variants(e)
            for spec, _ in e["variants"]:
                if spec["collocate"] not in ("vino", "local"):
                    sys.exit(f"experiment {e['label']}: collocate must be 'vino' or 'local'")
            for spec, _ in e["variants"]:
                if (spec["type"] == "pf" and not spec["numbers"]
                        and (uploads(spec) or collocated(spec))):
                    sys.exit(f"experiment {e['label']}: type pf needs `numbers` (e.g. 1/to/50)")
                if spec["control"] and spec["type"] != "pf":
                    sys.exit(f"experiment {e['label']}: `control: true` needs type pf")
            if any(p.upper for p in self.params):
                for spec, _ in e["variants"]:
                    if not (spec["grid"] or self.upper["grid"]):
                        sys.exit(f"experiment {e['label']}: pressure-level parameters need a MARS "
                                 "grid, set upper_air.grid (e.g. O1280) or the experiment's grid")
            e["by_date"] = {d: spec for spec, ds in e["variants"] for d in ds}
            self.experiments.append(e)
        labels = [e["label"] for e in self.experiments]
        if len(set(labels)) != len(labels):
            sys.exit("experiment labels must be unique")
        if not self.experiments:
            sys.exit("no experiments configured")
        # ensemble experiments add fair CRPS maps and scores
        self.ensemble = any(spec["type"] == "pf" for e in self.experiments for spec, _ in e["variants"])
        self.metrics = {m: v for m, v in METRICS.items()
                        if (self.ensemble or not m.startswith("fcrps"))
                        and (any(p.seeps for p in self.params) or not m.startswith("seeps"))}

    def resolve_open_ends(self, raw_periods, experiments, refresh):
        """{period: end date} for periods whose end is "-": the last base date
        for which every experiment has the longest step in MARS. The result is
        kept in <workdir>/open_periods.json and MARS is asked again only when
        `refresh` is set (submit) or the period/experiments changed, so all
        stages and jobs of one submission use the same dates."""
        open_periods = {k: v for k, v in raw_periods.items() if OPEN_END in tokens(v)}
        if not open_periods:
            return {}
        cache_file = self.workdir / "open_periods.json"
        try:
            cache = json.loads(cache_file.read_text())
        except (OSError, ValueError):
            cache = {}
        ends, todo = {}, {}
        for per, spec in open_periods.items():
            start = expand_dates(tokens(spec)[0])[0]
            exps = [period_spec(deep_update(EXP_DEFAULTS, e), per) for e in experiments]
            key = json.dumps([str(spec), exps, self.times, max(self.steps)], sort_keys=True, default=str)
            hit = cache.get(per, {})
            if not refresh and hit.get("key") == key:
                ends[per] = dt.datetime.strptime(hit["end"], "%Y%m%d").date()
            else:
                todo[per] = (start, exps, key)
        # all MARS list requests at once (each takes tens of seconds)
        from concurrent.futures import ThreadPoolExecutor
        queries = [(per, e) for per, (_, exps, _) in todo.items() for e in exps]
        with ThreadPoolExecutor(max(1, min(8, len(queries)))) as pool:
            found = list(pool.map(lambda q: mars_last_date(self, q[1], todo[q[0]][0]), queries))
        for per, (start, exps, key) in todo.items():
            last = [d for (p, _), d in zip(queries, found) if p == per]
            for e, d in zip(exps, last):
                where = db_where(e)
                log(f"period {per}: {e.get('label', e['expver'])} ({e['expver']}) "
                    + (f"is in {where} up to {ymd(d)}" if d else
                       f"NOT usable: no step {max(self.steps)} from {ymd(start)} at every "
                       f"base time {self.times} in {where} (see the lines above)"))
            if None in last:
                sys.exit(f"period {per}: an experiment has no forecasts from {ymd(start)} "
                         f"at every base time {self.times}")
            ends[per] = min(last)
            log(f"period {per}: using {ymd(start)} to {ymd(ends[per])}")
            cache[per] = {"key": key, "end": ymd(ends[per]),
                          "checked": dt.datetime.now().strftime("%Y-%m-%d %H:%M")}
        if todo:
            self.workdir.mkdir(parents=True, exist_ok=True)
            tmp = cache_file.with_suffix(".tmp")
            tmp.write_text(json.dumps(cache, indent=1))
            os.replace(tmp, cache_file)
        return ends

    def exp_variants(self, e):
        """[(spec, dates)] for an experiment. Any setting (expver, stream, ...)
        can be given per period as {period: value}; spec has it resolved."""
        per_period = [k for k, v in e.items() if isinstance(v, dict) and k != "variants"]
        for k in per_period:
            unknown = set(map(str, e[k])) - set(self.periods)
            missing = set(self.periods) - set(map(str, e[k]))
            if unknown or missing:
                sys.exit(f"experiment {e['label']}: '{k}' is given per period but "
                         f"periods {sorted(missing)} are missing / {sorted(unknown)} are unknown")
        if not per_period:
            spec = dict(e, expver=str(e["expver"]))
            return [(spec, self.dates)]
        groups, owner = {}, {}
        for per, dates in self.periods.items():
            spec = {k: ({str(a): b for a, b in v.items()}[per] if k in per_period else v)
                    for k, v in e.items() if k != "variants"}
            spec["expver"] = str(spec["expver"])
            key = json.dumps(spec, sort_keys=True, default=str)
            for d in dates:
                if owner.setdefault(d, key) != key:
                    sys.exit(f"experiment {e['label']}: {ymd(d)} is in two periods with different settings")
            groups.setdefault(key, (spec, set()))[1].update(dates)
        return [(spec, sorted(ds)) for spec, ds in groups.values()]

    # ---- paths
    @property
    def obs_mars(self):
        """True when any parameter's observations come from MARS feedback."""
        return any(p.mars for p in self.params)

    @property
    def mars_params(self):
        """Parameters whose observations come from MARS feedback."""
        return [p for p in self.params if p.mars]

    def obs_dir(self, p):
        if p.mars:
            # named after everything that selects the observations
            o = self.obs
            sel = o["upper_air"] if p.upper else o
            parts = ["ofbtemp" if p.upper else "ofb"] + [f"{k}{v}" for k, v in sorted(o["mars"].items())]
            parts += ["win" + "x".join(map(str, o["window"])), "at" + "x".join(map(str, tokens(o["window_times"]))),
                      f"tol{sel['tolerance']}", "obstype" + "x".join(tokens(sel["obstype"]))]
            if sel["codetype"]:
                parts.append("codetype" + "x".join(tokens(sel["codetype"])))
            if o["active_only"]:
                parts.append("active")
            return self.store / "obs" / "_".join(str(x).replace("/", "x") for x in parts)
        return self.store / "obs" / "_".join(self.upper["sources"] if p.upper else self.sources)

    def fc_dir(self, spec, p=None):
        """Store directory of one experiment variant. It is named after
        everything that defines the data (not the run-specific label), so
        all configs that use the same forecasts share it. Pressure-level
        parameters (p.upper) are always collocated locally by vino_upper.py."""
        if p is not None and p.upper:
            parts = ["upper", spec["class"], spec["stream"], spec["expver"], spec["type"]]
            for k in ("anoffset", "origin", "numbers"):
                if spec[k] not in (None, ""):
                    parts.append(f"{k}{spec[k]}")
            parts += [f"grid{spec['grid'] or self.upper['grid']}",
                      f"{self.upper['interpolation']}_{self.obs_dir(p).name}"]
            return self.store.joinpath("fc", *(str(x).replace("/", "x").replace(" ", "") for x in parts))
        if p is not None and p.local:
            # picked at the stations of the MARS feedback observations
            parts = ["ofb", spec["class"], spec["stream"], spec["expver"], spec["type"]]
            for k in ("anoffset", "origin", "grid", "numbers"):
                if spec[k] not in (None, ""):
                    parts.append(f"{k}{spec[k]}")
            orog = orog_file(self, spec)
            if orog:
                parts.append("orog" + hashlib.md5(orog.encode()).hexdigest()[:8])
            parts.append(f"{self.interpolation}_{self.obs_dir(p).name}")
            return self.store.joinpath("fc", *(str(x).replace("/", "x").replace(" ", "") for x in parts))
        if spec["collocate"] == "vino":
            parts = ["vino", spec["class"], spec["stream"], vino_expver(spec), spec["type"]]
        else:
            parts = ["local", spec["class"], spec["stream"], spec["expver"], spec["type"]]
        for k in ("anoffset", "origin", "grid", "resol", "numbers"):
            if spec[k] not in (None, ""):
                parts.append(f"{k}{spec[k]}")
        if spec["collocate"] == "local":
            extra = list(self.genfc.get("extra_args") or []) + list(spec["genfc_args"] or [])
            if extra:
                parts.append("args" + hashlib.md5(" ".join(map(str, extra)).encode()).hexdigest()[:8])
        parts.append(f"{self.interpolation}_{'_'.join(self.sources)}")
        return self.store.joinpath("fc", *(str(x).replace("/", "x").replace(" ", "") for x in parts))

    def fc_path(self, e, p, d, t, s):
        """Forecast files of experiment e for base date d, time t, step s: one,
        or pf + cf for an ensemble with its control (None if any is absent)."""
        spec = e["by_date"].get(d)
        if not spec:
            return None
        files = [fc_file(self.fc_dir(ds, p), p, d, t, s) for ds in data_specs(spec)]
        return None if None in files else files

    def store_lock(self, root, name):
        """Lock for the files called `name` below store directory `root`."""
        rel = Path(root).relative_to(self.store)
        return store_lock(self.store / "locks" / rel / f"{name}.lock")

    def store_tmp(self):
        tmproot = self.store / "tmp"
        tmproot.mkdir(parents=True, exist_ok=True)
        return tempfile.TemporaryDirectory(dir=tmproot)

    @property
    def plot_dir(self):
        return self.workdir / "plots"

    def plot_periods(self):
        """{name: (dates, valid hour or None, station classes or None)}: the
        periods, plus one per period and plot.valid_hours entry (e.g. DJF_00UTC)
        with only the forecasts valid at that hour (UTC), for a view of the
        diurnal cycle, and one per plot.station_groups entry for each of these
        (e.g. DJF_wmo, DJF_00UTC_wmo) with only stations of those classes."""
        hours = [None] + expand_ints(self.plot["valid_hours"] or [])
        groups = [(None, None)] + [(g, list(c)) for g, c in (self.plot["station_groups"] or {}).items()]
        for g, c in groups[1:]:
            bad = set(c) - set(STATION_CLASSES)
            if bad:
                raise SystemExit(f"plot.station_groups.{g}: unknown station class(es) {sorted(bad)}, "
                                 f"use {list(STATION_CLASSES)}")
        out = {}
        for per, ds in self.periods.items():
            for h in hours:
                for g, c in groups:
                    parts = ([] if per == "all" and (h is not None or g) else [per]) \
                        + ([] if h is None else [f"{h:02d}UTC"]) + ([g] if g else [])
                    out["_".join(parts)] = (ds, h, c)
        return out

    def steps_for(self, p):
        # accumulated/period parameters need step >= period
        return [s for s in self.steps if s >= p.period]


def data_specs(spec):
    """The data sets that make up one experiment variant: the forecast itself
    and, for an ensemble with `control: true`, its control forecast (cf)."""
    out = [spec]
    if spec["control"]:
        out.append(dict(spec, type="cf", numbers=None, control=False))
    return out


def day_dir(root, d):
    return root / f"{d:%Y}" / f"{d:%m}" / f"{d:%d}"


def obs_file(cfg, p, vt):
    return day_dir(cfg.obs_dir(p), vt.date()) / f"{p.key}_obs_{vt:%Y%m%d%H}.geo"


def fc_file(root, p, d, t, s):
    """Forecast file for base date d, time t, step s (None if absent)."""
    hits = glob.glob(str(day_dir(root, d) / f"{p.key}_mod_{ymd(d)}{t:02d}_{s}_*.geo"))
    return Path(sorted(hits)[0]) if hits else None


def valid_time(d, t, s):
    return dt.datetime(d.year, d.month, d.day, t) + dt.timedelta(hours=s)


def orog_file(cfg, spec):
    """Orography file given to vino_genfc with --orog_from_file, if any."""
    args = list(map(str, list(cfg.genfc.get("extra_args") or []) + list(spec["genfc_args"] or [])))
    for i, a in enumerate(args[:-1]):
        if a == "--orog_from_file":
            return os.path.expanduser(args[i + 1])
    return None


def vino_expver(exp):
    return exp["expver_label"] or exp["expver"]


def getgeo_model_args(cfg, exp, p):
    a = ["--table", "model", "--parameter", p.param,
         "--class", exp["class"], "--stream", exp["stream"],
         "--expver", vino_expver(exp), "--type", exp["type"],
         "--interpolation", cfg.interpolation]
    if p.period:
        a += ["--period", p.period]
    if exp["anoffset"] is not None:
        a += ["--anoffset", exp["anoffset"]]
    if exp["origin"]:
        a += ["--origin", exp["origin"]]
    return a


# --------------------------------------------------------------------------
# upload (vino_genfc)
# --------------------------------------------------------------------------

def in_vino(cfg, exp, p, t, dates):
    """Return the subset of dates already present in VINO for (exp, p, t).

    A date counts as present when both the first and the last step can be
    read back from VINO.
    """
    steps = cfg.steps_for(p)
    probe = sorted({steps[0], steps[-1]})
    present = set()
    tmproot = cfg.workdir / "tmp"
    tmproot.mkdir(parents=True, exist_ok=True)
    for chunk in chunks(dates, 60):
        tmp = Path(tempfile.mkdtemp(dir=tmproot, prefix="probe_"))
        try:
            run([cfg.tools["getgeo"], *getgeo_model_args(cfg, exp, p),
                 "--dates", *map(ymd, chunk), "--times", t,
                 "--steps", *probe, "--outdir", tmp])
            for d in chunk:
                if all(fc_file(tmp, p, d, t, s) for s in probe):
                    present.add(d)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    return present


def uploads(spec):
    return spec["upload"] and spec["collocate"] == "vino"


def stage_upload(cfg, args):
    exps = [e for e in cfg.experiments if any(uploads(s) for s, _ in e["variants"])]
    if args.experiment:
        exps = [e for e in exps if e["label"] in args.experiment]
    if not exps:
        log("nothing to upload (no experiments with upload: true selected)")
        return 0
    failures = []
    for exp in exps:
        with lock(cfg.workdir / "locks" / f"upload_{exp['label']}.lock"):
            # one combination per (expver variant, parameter, time)
            combos = []
            for spec, dates in ((ds, dates) for vs, dates in exp["variants"] for ds in data_specs(vs)):
                if not uploads(spec):
                    continue
                for p in cfg.params:
                    if p.local:  # not in VINO; collocated locally
                        continue
                    if not cfg.steps_for(p):
                        log(f"{p.key}: no steps >= period {p.period}, skipped")
                        continue
                    combos += [(spec, dates, p, t) for t in cfg.times]
            # 1) what is already in VINO, 2) upload the rest; both in parallel
            if args.force:
                done = [set() for _ in combos]
            else:
                done = parallel_map(cfg, lambda c: in_vino(cfg, c[0], c[2], c[3], c[1]), combos)
            tasks = []
            for (spec, dates, p, t), have in zip(combos, done):
                todo = [d for d in dates if d not in have]
                log(f"[{exp['label']} {spec['expver']} {p.key} {t:02d}z] {len(have)} dates "
                    f"already in VINO, {len(todo)} to upload")
                tasks += [(spec, p, t, chunk) for chunk in chunks(todo, cfg.chunk_size)]
            for missing in parallel_map(cfg, lambda a: upload_chunk(cfg, *a), tasks):
                failures += missing
    if failures:
        log(f"upload finished with {len(failures)} missing date/time combinations "
            "(missing in MARS?); rerun to retry")
    else:
        log("upload finished, everything is in VINO")
    return 0


def upload_chunk(cfg, exp, p, t, chunk):
    """Upload one chunk of dates and return the (label, param, time, date)
    still missing from VINO. With `database: mars+fdb` the dates MARS does
    not have are uploaded from the FDB."""
    todo = list(chunk)
    attempts = db_attempts(exp)
    for i, spec in enumerate(attempts):
        if i:
            log(f"[{exp['label']} {p.key} {t:02d}z] {len(todo)} dates not in MARS, "
                "trying the FDB: " + " ".join(map(ymd, todo)))
        todo = upload_dates(cfg, spec, p, t, todo, warn=i == len(attempts) - 1)
        if not todo:
            break
    return [(exp["label"], p.key, t, ymd(d)) for d in todo]


def upload_dates(cfg, exp, p, t, chunk, warn=True):
    """Upload dates with vino_genfc and return those still missing from VINO."""
    steps = cfg.steps_for(p)
    rc = run(genfc_cmd(cfg, exp, p, t, chunk, steps))
    failed = []
    if rc:
        # vino_genfc stops at the first error (a date missing in MARS, or a
        # DODB race between concurrent calls on the "updated" table), leaving
        # the rest of the chunk unwritten or holding data from an earlier run
        # that in_vino cannot tell apart, so redo the dates one at a time
        if len(chunk) > 1:
            log(f"[{exp['label']} {p.key} {t:02d}z] rc={rc}, retrying "
                f"{ymd(chunk[0])}-{ymd(chunk[-1])} one date at a time")
            failed = [d for d in chunk if run(genfc_cmd(cfg, exp, p, t, [d], steps))]
        else:
            failed = list(chunk)
    # the store copies of these dates may predate what is in VINO now
    root = cfg.fc_dir(exp)
    with cfg.store_lock(root, f"{p.key}_{t:02d}"):
        for d in chunk:
            for s in steps:
                f = fc_file(root, p, d, t, s)
                if f:
                    f.unlink()
    present = in_vino(cfg, exp, p, t, chunk)
    missing = [d for d in chunk if d in failed or d not in present]
    if missing and warn:
        log(f"[{exp['label']} {p.key} {t:02d}z] WARNING rc={rc}, not in VINO after upload: "
            + " ".join(map(ymd, missing)))
    return missing


def parallel_map(cfg, fn, items):
    """fn over items with cfg.parallel threads (the work is in subprocesses)."""
    from concurrent.futures import ThreadPoolExecutor
    if cfg.parallel == 1 or len(items) < 2:
        return [fn(i) for i in items]
    with ThreadPoolExecutor(min(cfg.parallel, len(items))) as pool:
        return list(pool.map(fn, items))


def genfc_cmd(cfg, exp, p, t, dates, steps):
    g = cfg.genfc
    cmd = [cfg.tools["genfc"], "--driven_by_forecast_date",
           "--parameter", p.param,
           "--dates", *map(ymd, dates), "--times", t,
           "--steps", *steps,
           "--class", exp["class"], "--stream", exp["stream"],
           "--expver", exp["expver"], "--type", exp["type"],
           "--interpolation", cfg.interpolation,
           "--observations_sources", *cfg.sources]
    if p.period:
        cmd += ["--period", p.period]
    for k in ("anoffset", "origin", "expver_label", "grid", "resol"):
        if exp[k] not in (None, ""):
            cmd += [f"--{k}", exp[k]]
    if exp["numbers"] not in (None, ""):
        cmd += ["--numbers", "/".join(tokens(exp["numbers"]))]
    if exp["database"]:
        cmd += ["--mars_database", exp["database"]]
    if g.get("ignore_missing"):
        cmd.append("--ignore_missing")
    if g.get("get_orog_once"):
        cmd.append("--get_orog_once")
    if cfg.interpolation in ("nearest_land_point", "nearest_sea_point"):
        cmd.append("--get_lsm_once")
    return cmd + list(g.get("extra_args") or []) + list(exp["genfc_args"] or [])


# --------------------------------------------------------------------------
# collocate (vino_genfc -> geopoint files, no VINO upload)
# --------------------------------------------------------------------------

def collocated(spec, p=None):
    """Collocated locally: collocate: local, and every pressure-level parameter."""
    return spec["collocate"] == "local" or (p is not None and p.local)


def collocates(cfg, e):
    return any(collocated(s) for s, _ in e["variants"]) or any(p.local for p in cfg.params)


def stage_collocate(cfg, args):
    exps = [e for e in cfg.experiments if collocates(cfg, e)]
    if args.experiment:
        exps = [e for e in exps if e["label"] in args.experiment]
    if not exps:
        log("nothing to collocate (no experiments with collocate: local selected)")
        return 0
    if cfg.obs_mars:
        # the forecasts are picked at the stations of the observations
        parallel_map(cfg, lambda chunk: download_obs_mars(cfg, chunk),
                     list(chunks(mars_obs_missing(cfg), cfg.chunk_size * 2)))
    tasks = [(spec, dates, p, t) for e in exps for vs, dates in e["variants"]
             for spec in data_specs(vs)
             for p in cfg.params if cfg.steps_for(p) and collocated(spec, p) for t in cfg.times]
    failures = []
    for missing in parallel_map(cfg, lambda a: collocate_set(cfg, args.force, *a), tasks):
        failures += missing
    if failures:
        log(f"collocate finished with {len(failures)} missing date/time combinations "
            "(missing in MARS?); rerun to retry")
    else:
        log("collocate finished, everything is in the store")
    return 0


def collocate_set(cfg, force, spec, dates, p, t):
    """Collocate the dates of one (experiment variant, parameter, time) that
    are not in the store yet; return the ones still missing."""
    root = cfg.fc_dir(spec, p)
    steps = cfg.steps_for(p)
    to_store = upper_to_store if p.local else genfc_to_store
    if p.upper and not p.mars:
        # the forecasts are picked at the stations of the observations, so
        # these must be in the store first (MARS feedback: done in stage_collocate)
        needed = defaultdict(set)
        for d in dates:
            for s in steps:
                vt = valid_time(d, t, s)
                needed[vt.hour].add(vt.date())
        for hour, ds in sorted(needed.items()):
            download_obs(cfg, p, hour, ds)

    def complete(d):
        # at observation stations: only steps with observations can be collocated
        return all(fc_file(root, p, d, t, s) for s in steps
                   if not p.local or obs_file(cfg, p, valid_time(d, t, s)).exists())

    with cfg.store_lock(root, f"{p.key}_{t:02d}"):
        todo = list(dates) if force else [d for d in dates if not complete(d)]
        log(f"[{spec['label']} {spec['expver']} {p.key} {t:02d}z] {len(dates) - len(todo)} dates "
            f"already in the store, {len(todo)} to collocate")
        missing = todo
        # with `database: mars+fdb` the dates MARS does not have come from the FDB
        for i, sp in enumerate(db_attempts(spec)):
            if not missing:
                break
            if i:
                log(f"[{spec['label']} {p.key} {t:02d}z] {len(missing)} dates not in MARS, "
                    "trying the FDB: " + " ".join(map(ymd, missing)))
            failed = []
            for chunk in chunks(missing, cfg.chunk_size):
                if to_store(cfg, sp, p, t, chunk, steps) == 0:
                    continue
                # vino_genfc stops at the first error; redo the dates one at a time
                if len(chunk) > 1:
                    log(f"[{spec['label']} {p.key} {t:02d}z] retrying "
                        f"{ymd(chunk[0])}-{ymd(chunk[-1])} one date at a time")
                    failed += [d for d in chunk if to_store(cfg, sp, p, t, [d], steps)]
                else:
                    failed += chunk
            missing = [d for d in missing if d in failed or not complete(d)]
    if missing:
        log(f"[{spec['label']} {p.key} {t:02d}z] WARNING not collocated: "
            + " ".join(map(ymd, missing)))
    return [(spec["label"], p.key, t, ymd(d)) for d in missing]


def genfc_to_store(cfg, spec, p, t, dates, steps):
    """Run vino_genfc with file output only and move its results into the
    store under the same names vino_getgeo uses. Returns the exit code; the
    output of a failed call is discarded."""
    with cfg.store_tmp() as tmp:
        out, geo = Path(tmp) / "genfc", Path(tmp) / "geo"
        rc = run(genfc_cmd(cfg, spec, p, t, dates, steps)
                 + ["--store_in_database", 0, "--store_geopoint_files", 1, "--rootdir", out])
        if rc:
            return rc
        for gz in out.rglob("*.gz"):
            m = re.search(r"_mod_Global_(\d{10})_(\d+)_", gz.name)
            if not m:
                continue
            base, step = m.groups()
            d = dt.datetime.strptime(base, "%Y%m%d%H")
            name = f"{p.key}_mod_{base}_{step}_{spec['stream']}_{spec['expver']}_{spec['type']}.geo"
            target = day_dir(geo, d.date()) / name
            target.parent.mkdir(parents=True, exist_ok=True)
            genfc_geo_to_store(gz, target, p, spec, d, int(step))
        move_into(geo, cfg.fc_dir(spec))
    return 0


def upper_cmd(*args):
    """Command line of vino_upper.py, run with the modules it needs."""
    return ["bash", "-c", f'{UPPER_ENV}; exec python3 "$0" "$@"',
            SCRIPT.with_name("vino_upper.py"), *args]


def upper_to_store(cfg, spec, p, t, dates, steps):
    """Collocate a pressure-level parameter, or a surface parameter at the
    MARS feedback stations, with vino_upper.py into the store."""
    a = ["fc", "--parameter", p.param, "--key", p.key,
         "--dates", *map(ymd, dates), "--time", t, "--steps", *steps,
         "--class", spec["class"], "--stream", spec["stream"], "--expver", spec["expver"],
         "--type", spec["type"], "--obsdir", cfg.obs_dir(p), "--mars", cfg.tools["mars"]]
    if p.upper:
        a += ["--level", p.level, "--grid", spec["grid"] or cfg.upper["grid"],
              "--interpolation", cfg.upper["interpolation"]]
    else:  # surface: native grid unless the experiment sets one, model orography at the stations
        if spec["grid"]:
            a += ["--grid", spec["grid"]]
        if orog_file(cfg, spec):
            a += ["--orog_file", orog_file(cfg, spec)]
    for k in ("anoffset", "origin", "database"):
        if spec[k] not in (None, ""):
            a += [f"--{k}", spec[k]]
    if spec["numbers"] not in (None, ""):
        a += ["--numbers", "/".join(map(str, expand_ints(spec["numbers"])))]
    with cfg.store_tmp() as tmp:
        rc = run(upper_cmd(*a, "--outdir", tmp))
        if rc == 0:
            move_into(tmp, cfg.fc_dir(spec, p))
    return rc


def genfc_geo_to_store(src, target, p, spec, base, step):
    """Copy a vino_genfc geopoint file and add the metadata vino_getgeo
    files carry. The values and units are the same as in VINO."""
    import gzip
    with gzip.open(src, "rt") as fh:
        text = fh.read()
    meta = {"class": spec["class"], "stream": spec["stream"], "expver": spec["expver"],
            "type": spec["type"], "anoffset": spec["anoffset"], "origin": spec["origin"],
            "param": p.param, "period": p.period, "date": f"{base:%Y%m%d}",
            "time": f"{base:%H}00", "step": step,
            "source": "vino_genfc (local collocation)"}
    lines = [f"{k}={v}" for k, v in meta.items() if v not in (None, "")]
    if "#METADATA" not in text:
        text = text.replace("#DATA\n", "#METADATA\n" + "\n".join(lines) + "\n#DATA\n", 1)
    target.write_text(text)


# --------------------------------------------------------------------------
# download (vino_getgeo)
# --------------------------------------------------------------------------

def stage_download(cfg, args):
    with lock(cfg.workdir / "locks" / "download.lock"):
        tasks = []
        if cfg.obs_mars:
            tasks += [(download_obs_mars, (cfg, chunk))
                      for chunk in chunks(mars_obs_missing(cfg), cfg.chunk_size * 2)]
        for p in cfg.params:
            if not p.local:
                tasks += [(download_obs, (cfg, p, hour, dates))
                          for hour, dates in sorted(obs_needed(cfg, p).items())]
            for exp in cfg.experiments:
                tasks += [(download_fc, (cfg, spec, vdates, p, t))
                          for vs, vdates in exp["variants"] for spec in data_specs(vs)
                          if not collocated(spec, p)
                          for t in cfg.times]
        log(f"checking {len(tasks)} observation/forecast sets, {cfg.parallel} at a time")
        parallel_map(cfg, lambda task: task[0](*task[1]), tasks)
    log("download finished")
    print_status(cfg)
    return 0


def fetch(cfg, cmd, root):
    """Run vino_getgeo into a temporary directory of the store, then move the
    complete files into place."""
    with cfg.store_tmp() as tmp:
        run(cmd + ["--outdir", tmp])
        move_into(tmp, root)


def obs_needed(cfg, p):
    """{hour: dates} of the observation valid times this config uses."""
    needed = defaultdict(set)
    for d in cfg.dates:
        for t in cfg.times:
            for s in cfg.steps_for(p):
                vt = valid_time(d, t, s)
                needed[vt.hour].add(vt.date())
    return needed


def download_obs(cfg, p, hour, dates):
    root = cfg.obs_dir(p)
    with cfg.store_lock(root, f"{p.key}_{hour:02d}"):
        todo = [d for d in sorted(dates)
                if not obs_file(cfg, p, dt.datetime(d.year, d.month, d.day, hour)).exists()]
        if not todo:
            return
        log(f"[obs {p.key} {hour:02d}z] {len(todo)} files to download")
        for chunk in chunks(todo, cfg.chunk_size * 4):
            if p.upper:
                fetch(cfg, upper_cmd("obs", "--parameter", p.param, "--level", p.level, "--key", p.key,
                                     "--dates", *map(ymd, chunk), "--hour", hour,
                                     "--sources", *cfg.upper["sources"]), root)
                continue
            cmd = [cfg.tools["getgeo"], "--parameter", p.param,
                   "--dates", *map(ymd, chunk), "--times", hour, "--sources", *cfg.sources]
            if p.period:
                cmd += ["--period", p.period]
            fetch(cfg, cmd, root)


def mars_obs_missing(cfg):
    """Valid times for which an observation file from MARS feedback is missing."""
    missing = set()
    for p in cfg.mars_params:
        missing |= {vt for hour, dates in obs_needed(cfg, p).items() for d in dates
                    for vt in [dt.datetime(d.year, d.month, d.day, hour)]
                    if not obs_file(cfg, p, vt).exists()}
    return sorted(missing)


def obs_tolerance(cfg):
    """Largest tolerance (minutes) of the observations this config uses."""
    tol = [float(cfg.obs["tolerance"])]
    if any(p.upper for p in cfg.mars_params):
        tol.append(float(cfg.obs["upper_air"]["tolerance"]))
    return max(tol)


def obs_windows(cfg, vt):
    """Analysis (date, time) of the assimilation windows that hold reports
    within the observation tolerance of valid time vt. The window of analysis
    time T covers (T + window[0], T + window[1]] hours."""
    o = cfg.obs
    tol = dt.timedelta(minutes=obs_tolerance(cfg))
    w0, w1 = (dt.timedelta(hours=float(h)) for h in o["window"])
    out = set()
    for days in range(-2, 3):
        for h in expand_ints(o["window_times"]):
            an = dt.datetime.combine(vt.date() + dt.timedelta(days=days), dt.time(h))
            if an + w0 < vt + tol and an + w1 >= vt - tol:
                out.add(an)
    return out


def ofb_where(cfg):
    """SQL condition selecting the feedback reports of all parameters: SYNOP
    for the surface ones, radiosondes at the configured levels for the others."""
    o, conds = cfg.obs, []
    sfc = [p for p in cfg.mars_params if not p.upper]
    upper = [p for p in cfg.mars_params if p.upper]
    if sfc:
        c = (f"obstype in ({','.join(tokens(o['obstype']))}) and "
             f"varno in ({','.join(map(str, sorted({OFB_VARNO[p.param] for p in sfc})))})")
        if o["codetype"]:
            c += f" and codetype in ({','.join(tokens(o['codetype']))})"
        conds.append(c)
    if upper:
        u = o["upper_air"]
        c = (f"obstype in ({','.join(tokens(u['obstype']))}) and vertco_type=1 and "
             f"varno in ({','.join(map(str, sorted({OFB_UPPER_VARNO[p.param] for p in upper})))}) and "
             f"vertco_reference_1 in ({','.join(str(lev * 100) for lev in sorted({p.level for p in upper}))})")
        if u["codetype"]:
            c += f" and codetype in ({','.join(tokens(u['codetype']))})"
        conds.append(c)
    return " or ".join(f"({c})" for c in conds)


def nearest_time(times, t):
    """Element of the sorted datetime64 array `times` closest to each of t."""
    i = np.searchsorted(times, t).clip(1, max(1, len(times) - 1))
    if len(times) == 1:
        return np.repeat(times, len(t))
    before = np.abs(t - times[i - 1]) <= np.abs(t - times[i])
    return times[np.where(before, i - 1, i)]


def download_obs_mars(cfg, vts):
    """Observations of all parameters at the valid times vts from MARS
    feedback (observations.source: mars) into the store: SYNOP for the
    surface parameters, radiosondes (TEMP) for the pressure levels. One MARS
    request per analysis time fetches all windows the valid times need."""
    import pyodc
    o = cfg.obs
    params = cfg.mars_params
    with cfg.store_lock(cfg.store / "obs", f"mars_{vts[0]:%Y%m%d%H}"):
        vts = [vt for vt in vts if any(not obs_file(cfg, p, vt).exists() for p in params)]
        if not vts:
            return
        windows = set().union(*(obs_windows(cfg, vt) for vt in vts))
        log(f"[obs mars] {len(vts)} valid times {vts[0]:%Y%m%d%H}-{vts[-1]:%Y%m%d%H}, "
            f"{len(windows)} assimilation windows to retrieve")
        sql = ("select statid,lat,lon,stalt,date,time,obstype,varno,vertco_reference_1,"
               f"obsvalue,datum_status@body where {ofb_where(cfg)}")
        # the feedback (GBs for a month) goes to scratch, the geopoint files to the store
        scratch = os.environ.get("SCRATCHDIR") or os.environ.get("TMPDIR") or None
        with tempfile.TemporaryDirectory(dir=scratch) as tmp, cfg.store_tmp() as out:
            tmp, out = Path(tmp), Path(out)
            req, targets = [], []
            for h in sorted({an.hour for an in windows}):
                dates = sorted({an.date() for an in windows if an.hour == h})
                targets.append(tmp / f"ofb_{h:02d}.odb")
                fields = dict(o["mars"], date="/".join(map(ymd, dates)), time=f"{h:02d}",
                              filter=f'"{sql}"', target=f'"{targets[-1]}"')
                req.append("retrieve," + ",".join(f"{k}={v}" for k, v in fields.items()))
            (tmp / "request").write_text("\n".join(req) + "\n")
            rc = run([cfg.tools["mars"], tmp / "request"])
            if rc:
                log(f"[obs mars] WARNING MARS rc={rc}, using whatever was retrieved")
            # keep only the reports near the wanted valid times (the windows hold every report)
            want = np.array(sorted(vts), dtype="datetime64[ns]")
            upper_types = {int(x) for x in tokens(o["upper_air"]["obstype"])}
            frames = []
            for f in targets:
                if not f.exists() or not f.stat().st_size:
                    continue
                with open(f, "rb") as fh:
                    for frame in pyodc.Reader(fh).frames:
                        df = frame.dataframe()
                        df.columns = [c.split("@")[0] for c in df.columns]
                        df = df[np.isfinite(df["obsvalue"]) & (df["obsvalue"].abs() < MISSING)]
                        if o["active_only"]:
                            df = df[(df["datum_status"].fillna(0).astype(int) & 1) == 1]
                        when = pd.to_datetime(df["date"].astype(str) + df["time"].astype(int).map("{:06d}".format),
                                              format="%Y%m%d%H%M%S")
                        df["vt"] = nearest_time(want, when.to_numpy())
                        df["off"] = np.abs((when.to_numpy() - df["vt"].to_numpy()) / np.timedelta64(1, "m"))
                        df["upper"] = df["obstype"].astype(int).isin(upper_types) & any(p.upper for p in params)
                        tol = np.where(df["upper"], float(o["upper_air"]["tolerance"]), float(o["tolerance"]))
                        frames.append(df[df["off"] <= tol])
            if not frames:
                log("[obs mars] WARNING no observations retrieved")
                return
            df = pd.concat(frames, ignore_index=True)
            # station ids as VINO has them: no padding, no leading zeros in WMO numbers
            sid = df["statid"].astype(str).str.strip()
            df["stnid"] = sid.where(~sid.str.fullmatch(r"\d+"), sid.str.lstrip("0").replace("", "0"))
            df["level"] = np.where(df["upper"], (df["vertco_reference_1"] / 100).round(), 0).astype(int)
            # one report per station: the closest in time, active ones first
            # (e.g. a radiosonde comes both as TAC and as BUFR report)
            df["passive"] = (df["datum_status"].fillna(0).astype(int) & 1) == 0
            df = df.sort_values(["off", "passive"]).drop_duplicates(["vt", "upper", "level", "varno", "stnid"])
            for p in params:
                if p.upper:
                    sel = df[df["upper"] & (df["level"] == p.level) & (df["varno"] == OFB_UPPER_VARNO[p.param])]
                else:
                    sel = df[~df["upper"] & (df["varno"] == OFB_VARNO[p.param])]
                for vt, g in sel.groupby("vt"):
                    vt = vt.to_pydatetime()
                    if obs_file(cfg, p, vt).exists():
                        continue
                    write_obs_geo(out / obs_file(cfg, p, vt).relative_to(cfg.store), p, vt, g)
            n = move_into(out, cfg.store)
        log(f"[obs mars] {n} observation files written")


def write_obs_geo(path, p, vt, df):
    """Observation geopoint file laid out like those of vino_getgeo."""
    units = {"2t": "K", "2d": "K", "10dd": "deg", "t": "K", "td": "K", "z": "m2/s2",
             "dd": "deg"}.get(p.param, "m/s")
    level, ltype = (p.level, "pl") if p.upper else (0, "sfc")
    lines = ["#GEO", "#FORMAT NCOLS", "#COLUMNS",
             "stnid\tlatitude\tlongitude\tlevel\tdate\ttime\televation\tvalue_0",
             "# Missing values represented by 3e+38", "#METADATA",
             f"date={vt:%Y%m%d}", f"level={level}", f"leveltype={ltype}", f"param={p.param}", "period=0",
             "step=0", f"time={vt:%H}00", f"units={units}", "source=MARS feedback (ofb)", "#DATA"]
    elev = df["stalt"].where(np.isfinite(df["stalt"]) & (df["stalt"].abs() < MISSING))
    for sid, lat, lon, z, v in zip(df["stnid"], df["lat"], df["lon"], elev, df["obsvalue"]):
        z = "3e+38" if pd.isna(z) else f"{z:.0f}"
        lines.append(f"{sid}\t{lat:.4f}\t{lon:.4f}\t{level}\t{vt:%Y%m%d}\t{vt:%H}00\t{z}\t{v:.6g}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def download_fc(cfg, spec, vdates, p, t):
    root = cfg.fc_dir(spec)
    with cfg.store_lock(root, f"{p.key}_{t:02d}"):
        missing = [(d, s) for d in vdates for s in cfg.steps_for(p)
                   if fc_file(root, p, d, t, s) is None]
        if not missing:
            return
        dates = sorted({d for d, _ in missing})
        steps = sorted({s for _, s in missing})
        log(f"[{spec['label']} {spec['expver']} {p.key} {t:02d}z] {len(missing)} files to download")
        for chunk in chunks(dates, cfg.chunk_size):
            fetch(cfg, [cfg.tools["getgeo"], *getgeo_model_args(cfg, spec, p),
                        "--dates", *map(ymd, chunk), "--times", t, "--steps", *steps], root)


def print_status(cfg):
    log(f"data in the store {cfg.store}")
    for p in cfg.params:
        steps = cfg.steps_for(p)
        n = len(cfg.dates) * len(cfg.times) * len(steps)
        vts = {valid_time(d, t, s) for d in cfg.dates for t in cfg.times for s in steps}
        nobs = sum(obs_file(cfg, p, vt).exists() for vt in vts)
        print(f"  {p.key:8s} obs        {nobs:6d} / {len(vts)} files")
        for exp in cfg.experiments:
            nfc = sum(cfg.fc_path(exp, p, d, t, s) is not None
                      for d in cfg.dates for t in cfg.times for s in steps)
            print(f"  {p.key:8s} {exp['label']:10s} {nfc:6d} / {n} files")
    for exp in cfg.experiments:
        for vs, _ in exp["variants"]:
            for spec in data_specs(vs):
                print(f"  {exp['label']} {spec['expver']} {spec['type']}: {cfg.fc_dir(spec)}")
                if any(p.upper for p in cfg.params):
                    print(f"  {exp['label']} {spec['expver']} {spec['type']} (pressure levels): "
                          f"{cfg.fc_dir(spec, next(p for p in cfg.params if p.upper))}")
                if any(not p.upper for p in cfg.mars_params):
                    print(f"  {exp['label']} {spec['expver']} {spec['type']} (MARS feedback stations): "
                          f"{cfg.fc_dir(spec, next(p for p in cfg.mars_params if not p.upper))}")


def stage_status(cfg, args):
    print_status(cfg)
    return 0


# --------------------------------------------------------------------------
# plot
# --------------------------------------------------------------------------

@lru_cache(maxsize=512)
def read_geo(path):
    """Read a Metview NCOLS geopoint file -> (DataFrame, metadata dict)."""
    meta, cols, skip = {}, None, 0
    with open(path) as fh:
        section = None
        for i, line in enumerate(fh):
            line = line.rstrip("\n")
            if line.startswith("#COLUMNS"):
                section = "columns"
            elif line.startswith("#METADATA"):
                section = "meta"
            elif line.startswith("#DATA"):
                skip = i + 1
                break
            elif line.startswith("#"):
                continue
            elif section == "columns":
                cols = line.split()
                section = None
            elif section == "meta" and "=" in line:
                k, v = line.split("=", 1)
                meta[k.strip()] = v.strip()
    df = pd.read_csv(path, sep=r"\s+", skiprows=skip, names=cols, header=None,
                     dtype={"stnid": str}, engine="c")
    # one value column per ensemble member (value_0 ... value_<n-1>)
    vals = sorted((c for c in df.columns if re.fullmatch(r"value_\d+", str(c))),
                  key=lambda c: int(c[6:])) or ["value_0"]
    for c in vals + ["elevation"]:
        if c not in df:
            df[c] = np.nan
        df[c] = pd.to_numeric(df[c], errors="coerce")
        df[c] = df[c].where(df[c].abs() < MISSING)
    df["elevation"] = df["elevation"].where(df["elevation"] <= 9000)  # obs use 99999
    df = df.dropna(subset=vals)
    return df[["stnid", "latitude", "longitude", "elevation"] + vals], meta


STATION_CLASSES = ("wmo", "icao", "noid", "other")


def station_class(ids):
    """Identifier class of station ids (as in the store: WMO numbers without
    leading zeros): wmo = digits only, icao = 4 letters (METAR), noid = no
    identifier (the ODB gives these reports '_' + a hash, e.g. national
    networks reported without a WMO number), other = anything else (ships,
    mobile stations)."""
    ids = pd.Series(ids, dtype=str)
    return pd.Series(np.select([ids.str.fullmatch(r"\d+"), ids.str.fullmatch(r"[A-Z]{4}"),
                                ids.str.startswith("_") | (ids.str.strip() == "")],
                               ["wmo", "icao", "noid"], "other"), index=ids.index)


def key_class(keys):
    """station_class of station keys (stnid_lat_lon)."""
    return station_class(pd.Series(keys, dtype=str).str.rsplit("_", n=2).str[0])


def station_key(df):
    return (df["stnid"].astype(str) + "_" + df["latitude"].round(2).astype(str)
            + "_" + df["longitude"].round(2).astype(str))


def matched(cfg, p, step, dates=None):
    """Forecast/observation pairs for one parameter and step, all experiments.

    Returns a DataFrame with columns key, lat, lon, base, obs and one forecast
    column per experiment, restricted to the sample common to all
    experiments, plus the units string. Only pairs where the observation and
    every forecast are valid are kept, with observations outside OBS_RANGE
    rejected. For 2t (and tmax/tmin) the forecasts are corrected from model
    orography to station height with plot.lapse_rate. Stations whose model and
    station elevations are both known and differ by more than plot.max_dz are
    dropped, and so are stations whose identifier class (station_class) is not
    in plot.station_ids (when set). This follows the processing of quaver's error_map. The columns
    dz_<label> hold the model minus station elevation (NaN when unknown).
    For ensembles (several members per file) the forecast column is the
    ensemble mean, and crps_<label> holds the fair CRPS of each pair (the
    absolute error for single forecasts).
    """
    frames, units = [], ""
    lapse = float(cfg.plot["lapse_rate"] or 0) if p.param in LAPSE_PARAMS else 0.0
    max_dz = cfg.plot["max_dz"]
    keep_ids = cfg.plot.get("station_ids")
    if keep_ids is not None and set(keep_ids) - set(STATION_CLASSES):
        raise SystemExit(f"plot.station_ids: unknown class(es) {sorted(set(keep_ids) - set(STATION_CLASSES))}, "
                         f"use {list(STATION_CLASSES)}")
    for d in cfg.dates if dates is None else dates:
        for t in cfg.times:
            vt = valid_time(d, t, step)
            of = obs_file(cfg, p, vt)
            if not of.exists():
                continue
            fcs = [cfg.fc_path(e, p, d, t, step) for e in cfg.experiments]
            if any(f is None for f in fcs):
                continue
            obs, ometa = read_geo(str(of))
            m = obs.assign(key=station_key(obs)).rename(columns={"value_0": "obs", "elevation": "z_obs"})
            m = m.drop_duplicates("key")[["key", "latitude", "longitude", "z_obs", "obs"]]
            if keep_ids is not None:
                m = m[key_class(m["key"]).isin(keep_ids).to_numpy()]
            labels = [e["label"] for e in cfg.experiments]
            members = {}
            for lab, files in zip(labels, fcs):
                # members of all files of the experiment (pf, plus cf with control: true)
                members[lab] = []
                for i, f in enumerate(files):
                    fc, meta = read_geo(str(f))
                    units = units or meta.get("units", "") or ometa.get("units", "")
                    fc = fc.assign(key=station_key(fc)).drop_duplicates("key")
                    vals = [c for c in fc.columns if c.startswith("value_")]
                    names = [f"{lab}#{len(members[lab]) + j}" for j in range(len(vals))]
                    members[lab] += names
                    cols = {"elevation": f"z_{lab}", **dict(zip(vals, names))}
                    keep = ["key"] + (["elevation"] if i == 0 else []) + vals
                    m = m.merge(fc[keep].rename(columns=cols), on="key", how="inner")
            need = ["obs"] + [c for lab in labels for c in members[lab]]
            if p.param in OBS_RANGE:
                lo, hi = OBS_RANGE[p.param]
                m = m[m["obs"].between(lo, hi)]
            if lapse:
                need += ["z_obs"] + [f"z_{lab}" for lab in labels]
            if max_dz is not None:
                for lab in labels:
                    m = m[~((m[f"z_{lab}"] - m["z_obs"]).abs() > float(max_dz))]
            if lapse:
                for lab in labels:
                    # model surface above the station -> station is warmer
                    m[members[lab]] = m[members[lab]].add(lapse * (m[f"z_{lab}"] - m["z_obs"]), axis=0)
            m = m.dropna(subset=need)
            if p.param == "z" and p.upper:  # geopotential -> geopotential height
                cols = ["obs"] + [c for lab in labels for c in members[lab]]
                m[cols] = m[cols] / G
                units = "gpm"
            obs_v = m["obs"].to_numpy(float)
            for lab in labels:
                ens = m[members[lab]].to_numpy(float)
                m[lab] = ens_mean(p, ens)
                m[f"crps_{lab}"] = fair_crps(p, ens, obs_v)
                m[f"dz_{lab}"] = m[f"z_{lab}"] - m["z_obs"]
            m = m[["key", "latitude", "longitude", "obs"] + labels
                  + [f"crps_{lab}" for lab in labels] + [f"dz_{lab}" for lab in labels]]
            m["base"] = f"{ymd(d)}{t:02d}"
            frames.append(m)
    if not frames:
        return None, units
    return pd.concat(frames, ignore_index=True), units


def error(p, fc, obs):
    e = fc - obs
    if p.param in ANGLE_PARAMS:  # wind direction: shortest angular difference
        e = (e + 180.0) % 360.0 - 180.0
    return e


def ens_mean(p, ens):
    """Ensemble mean of the members (columns of ens); circular for 10dd."""
    if ens.shape[1] == 1:
        return ens[:, 0]
    if p.param in ANGLE_PARAMS:
        r = np.deg2rad(ens)
        return np.rad2deg(np.arctan2(np.sin(r).mean(1), np.cos(r).mean(1))) % 360.0
    return ens.mean(1)


def fair_crps(p, ens, obs):
    """Fair CRPS (Ferro 2014) of each row of ens (members in columns):
    mean |x_i - y| - sum_{i<j} |x_i - x_j| / (m (m - 1)).
    For a single member it is the absolute error."""
    m = ens.shape[1]
    skill = np.abs(error(p, ens, obs[:, None])).mean(1)
    if m == 1:
        return skill
    if p.param in ANGLE_PARAMS:
        spread = np.abs(error(p, ens[:, :, None], ens[:, None, :])).sum((1, 2)) / 2
    else:
        # sum_{i<j} |x_i - x_j| = sum_k (2k - m - 1) x_(k) over the sorted members
        spread = (np.sort(ens, axis=1) * (2 * np.arange(1, m + 1) - m - 1)).sum(1)
    return skill - spread / (m * (m - 1))


# SEEPS climatology from the STVL database (DODB), as read by Quaver's
# stationclimatology (climat rodw, 1980-2009, monthly per station). Per
# station: value_0 p1 (dry probability), value_1 dry/light threshold,
# value_3 light/heavy threshold, value_4..12 the 3x3 SEEPS error matrix
# (row = forecast category, column = observed category).
SEEPS_CLIMATE = "1980-2009"
SEEPS_P1 = (0.1, 0.85)  # stations used where lower <= p1 < upper (Quaver seeps_filter)
SEEPS_FETCH = r'''
import sys, pandas as pd
from vtb.media import dodb
out = []
for mo in range(1, 13):
    req = dict(table="climatology", parameter="tp", level_type="sfc", level=0, period=86400,
               reference_datetimes=pd.to_datetime([f"2001-{mo:02d}-15"]), forecast_lengths=[0],
               model={"class": "od", "stream": "oper", "type": "ec", "climat": "rodw",
                      "averaged": "monthlystation"},
               climatology={"climate_period": sys.argv[2], "category": "equitable_categories"})
    for c in dodb.DODB().read(**req):
        df = dodb.DODB.container_to_geodf(c)
        out.append(pd.DataFrame(df).assign(month=mo))
pd.concat(out, ignore_index=True).to_csv(sys.argv[1], index=False)
'''


def seeps_climatology(cfg):
    """Station SEEPS climatology of all months (fetched once into the store):
    one row per stnid and month with p1, t1, t2 and m0..m8."""
    path = cfg.store / "clim" / f"seeps_tp24_{SEEPS_CLIMATE}.csv"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(f".{os.getpid()}.tmp")
        log(f"fetching the SEEPS station climatology {SEEPS_CLIMATE} from STVL")
        subprocess.run(["bash", "-c", f'{UPPER_ENV}; exec python3 -c "$0" "$@"',
                        SEEPS_FETCH, str(tmp), SEEPS_CLIMATE], check=True)
        os.replace(tmp, path)
    c = pd.read_csv(path, dtype={"stnid": str})
    # several records for a station: keep the one from most observations (as Quaver)
    c = c.sort_values(["obscount", "enddate"]).drop_duplicates(["stnid", "month"], keep="last")
    c = c.rename(columns={"value_0": "p1", "value_1": "t1", "value_3": "t2",
                          **{f"value_{i + 4}": f"m{i}" for i in range(9)}})
    return c[["stnid", "month", "p1", "t1", "t2"] + [f"m{i}" for i in range(9)]]


def add_seeps(df, clim, step, labels):
    """Columns seeps_<label>: SEEPS of each pair (NaN without a usable
    climatology), as Quaver's seepsField: values rounded to 0.1 mm, categories
    dry (<= t1), light, heavy (>= t2), the climatology of the valid month."""
    vt = pd.to_datetime(df["base"], format="%Y%m%d%H") + pd.Timedelta(hours=step)
    m = df.assign(stnid=df["key"].str.split("_").str[0], month=vt.dt.month)
    m = m[["stnid", "month"]].merge(clim, on=["stnid", "month"], how="left")
    lo, hi = SEEPS_P1
    ok = ((m["p1"] >= lo) & (m["p1"] < hi)).to_numpy()
    mat = np.nan_to_num(m[[f"m{i}" for i in range(9)]].to_numpy(float))
    t1, t2 = m["t1"].to_numpy(float), m["t2"].to_numpy(float)

    def cat(v):
        v = np.round(np.asarray(v, float), 1)
        return (v > t1).astype(int) + (v >= t2).astype(int)

    ob = cat(df["obs"])
    out = df.copy()
    for lab in labels:
        idx = cat(df[lab]) * 3 + ob
        out[f"seeps_{lab}"] = np.where(ok, mat[np.arange(len(m)), idx], np.nan)
    return out


def in_region(df, extent):
    if not extent:
        return df
    w, e, s, n = extent
    return df[(df.longitude >= w) & (df.longitude <= e)
              & (df.latitude >= s) & (df.latitude <= n)]


def station_stats(cfg, p, df, nbases):
    """Per-station bias/rmse/mae per experiment, and rmse/mae differences.

    As in quaver's error_map, a station is kept only when it has valid pairs
    for at least plot.min_fraction of the nbases forecasts.
    """
    labels = [e["label"] for e in cfg.experiments]
    work = df[["key", "latitude", "longitude"] + [f"dz_{lab}" for lab in labels]].copy()
    for lab in labels:
        err = error(p, df[lab], df["obs"])
        work[f"bias_{lab}"] = err
        work[f"sq_{lab}"] = err ** 2
        work[f"abs_{lab}"] = err.abs()
        work[f"crps_{lab}"] = df[f"crps_{lab}"]
    seeps = f"seeps_{labels[0]}" in df
    if seeps:
        for lab in labels:
            work[f"seeps_{lab}"] = df[f"seeps_{lab}"]
    g = work.groupby("key")
    st = g[["latitude", "longitude"]].first()
    st["n"] = g.size()
    for lab in labels:
        st[f"dz_{lab}"] = g[f"dz_{lab}"].mean()
    for lab in labels:
        st[f"bias_{lab}"] = g[f"bias_{lab}"].mean()
        st[f"rmse_{lab}"] = np.sqrt(g[f"sq_{lab}"].mean())
        st[f"mae_{lab}"] = g[f"abs_{lab}"].mean()
        st[f"fcrps_{lab}"] = g[f"crps_{lab}"].mean()
        if seeps:  # stations without a usable climatology stay NaN
            st[f"seeps_{lab}"] = g[f"seeps_{lab}"].mean()
    ref = labels[0]
    if seeps:
        for lab in labels[1:]:
            st[f"seepsdiff_{lab}"] = st[f"seeps_{lab}"] - st[f"seeps_{ref}"]
    for lab in labels[1:]:
        st[f"mef_{lab}"] = st[f"bias_{lab}"] - st[f"bias_{ref}"]
        st[f"rmsediff_{lab}"] = st[f"rmse_{lab}"] - st[f"rmse_{ref}"]
        st[f"maediff_{lab}"] = st[f"mae_{lab}"] - st[f"mae_{ref}"]
        rel = 100.0 * (st[f"rmse_{lab}"] / st[f"rmse_{ref}"] - 1.0)
        st[f"rmserel_{lab}"] = rel.replace([np.inf, -np.inf], np.nan)
        st[f"fcrpsdiff_{lab}"] = st[f"fcrps_{lab}"] - st[f"fcrps_{ref}"]
        rel = 100.0 * (st[f"fcrps_{lab}"] / st[f"fcrps_{ref}"] - 1.0)
        st[f"fcrpsrel_{lab}"] = rel.replace([np.inf, -np.inf], np.nan)
    min_count = max(1, float(cfg.plot["min_fraction"]) * nbases)
    return st[st["n"] >= min_count].reset_index()


_DENSITY_CACHE = {}


def station_density_weights(lat, lon, scangle=0.75):
    """Station-density weights as quaver's spatial_mean_weights=station_density
    (vtb.geo.station_density_weights): w_i = 1 / sum_j exp(-(d_ij / scangle)^2)
    over all stations j within 4 * scangle (great-circle distance, degrees),
    normalised to sum 1. Isolated stations get weight ~1, a station among k
    close neighbours ~1/k."""
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    key = (float(scangle), hashlib.md5(np.round(np.c_[lat, lon], 5).tobytes()).hexdigest())
    if key in _DENSITY_CACHE:
        return _DENSITY_CACHE[key]
    sc = np.deg2rad(float(scangle))
    cutoff = 4.0 * sc
    la, lo = np.deg2rad(lat), np.deg2rad(lon)
    xyz = np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
    dens = np.zeros(len(lat))
    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(xyz)
        nbrs = tree.query_ball_point(xyz, 2.0 * np.sin(0.5 * cutoff) * (1 + 1e-9))
    except ImportError:
        nbrs = None
    block = 512
    for i0 in range(0, len(lat), block):
        rows = xyz[i0:i0 + block]
        if nbrs is None:  # blocked brute force: fine for O(10^4) stations
            chord = np.sqrt(np.maximum(((rows[:, None, :] - xyz[None, :, :]) ** 2).sum(-1), 0.0))
            d = 2.0 * np.arcsin(np.minimum(chord / 2.0, 1.0))
            k = np.where(d > cutoff, 0.0, np.exp(-(d / sc) ** 2))
            dens[i0:i0 + block] = k.sum(1)
        else:
            for j, idx in enumerate(nbrs[i0:i0 + block]):
                chord = np.linalg.norm(xyz[idx] - rows[j], axis=1)
                d = 2.0 * np.arcsin(np.minimum(chord / 2.0, 1.0))
                dens[i0 + j] = np.exp(-(d[d <= cutoff] / sc) ** 2).sum()
    w = 1.0 / dens
    w /= w.sum()
    _DENSITY_CACHE[key] = w
    return w


def station_weights(cfg, st):
    """Weights of the stations in the regional scores (plot.spatial_mean_weights)."""
    mode = str(cfg.plot.get("spatial_mean_weights") or "none").lower()
    if mode in ("none", "generic"):
        return np.full(len(st), 1.0 / max(len(st), 1))
    if mode == "station_density":
        return station_density_weights(st["latitude"], st["longitude"],
                                       cfg.plot.get("scangle", 0.75))
    sys.exit(f"plot.spatial_mean_weights: unknown value {mode!r} (none | station_density)")


def wmean(values, w):
    """Weighted mean over the stations with a valid value."""
    v = np.asarray(values, dtype=float)
    ok = np.isfinite(v)
    return float((v[ok] * w[ok]).sum() / w[ok].sum()) if ok.any() else np.nan


def overall_scores(cfg, st):
    """Regional scores as means over the station scores, like the domain means
    of quaver's error_map: RMSE = sqrt(mean of station RMSE^2). The stations are
    weighted equally or by plot.spatial_mean_weights (e.g. station_density)."""
    w = station_weights(cfg, st)
    out = {"n": int(st["n"].sum()), "stations": int(len(st))}
    for e in cfg.experiments:
        lab = e["label"]
        out[lab] = {"bias": wmean(st[f"bias_{lab}"], w),
                    "rmse": float(np.sqrt(wmean(st[f"rmse_{lab}"] ** 2, w))),
                    "mae": wmean(st[f"mae_{lab}"], w),
                    "fcrps": wmean(st[f"fcrps_{lab}"], w)}
        if f"seeps_{lab}" in st:
            out[lab]["seeps"] = wmean(st[f"seeps_{lab}"], w)
    return out


def compute_step(args):
    cfg_path, per, pi, step = args
    cfg = Config(cfg_path)
    p = cfg.params[pi]
    dates, hour, classes = cfg.plot_periods()[per]
    df, units = matched(cfg, p, step, dates)
    if df is not None and hour is not None:
        df = df[(df["base"].str[8:].astype(int) + step) % 24 == hour]
    if df is not None and classes is not None:
        df = df[key_class(df["key"]).isin(classes).to_numpy()]
    if df is None or df.empty:
        return per, pi, step, None
    if p.seeps:
        df = add_seeps(df, seeps_climatology(cfg), step, [e["label"] for e in cfg.experiments])
    res = {"units": units, "bases": sorted(df["base"].unique()), "regions": {}}
    for rname, extent in cfg.plot["regions"].items():
        sub = in_region(df, extent)
        if sub.empty:
            continue
        st = station_stats(cfg, p, sub, len(res["bases"]))
        if st.empty:
            continue
        res["regions"][rname] = {"stations": st,
                                 "scores": overall_scores(cfg, st),
                                 "daily": sub if cfg.plot["daily_maps"] else None}
    return per, pi, step, res


def nice_limit(values, q=98):
    v = np.abs(np.asarray(values, dtype=float))
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 1.0
    x = np.percentile(v, q)
    if x <= 0:
        return 1.0
    mag = 10 ** np.floor(np.log10(x))
    for f in (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10):
        if f * mag >= x:
            return f * mag
    return 10 * mag


METRICS = {
    # name: (colormap, diverging, description)
    "bias": ("RdBu_r", True, "mean error (fc - obs)"),
    "rmse": ("viridis", False, "root mean square error"),
    "mae": ("viridis", False, "mean absolute error"),
    "mef": ("RdBu_r", True, "mean bias difference vs reference (experiment - reference)"),
    "rmsediff": ("RdBu_r", True, "RMSE difference vs reference (red = worse)"),
    "maediff": ("RdBu_r", True, "MAE difference vs reference (red = worse)"),
    "rmserel": ("RdBu_r", True, "relative RMSE change vs reference (red = worse)"),
    # only with ensemble experiments; for single forecasts the CRPS is the MAE
    "fcrps": ("viridis", False, "fair CRPS"),
    "fcrpsdiff": ("RdBu_r", True, "fair CRPS difference vs reference (red = worse)"),
    "fcrpsrel": ("RdBu_r", True, "relative fair CRPS change vs reference (red = worse)"),
    # only with `seeps: true` (tp24); stations with 0.1 <= p1 < 0.85
    "seeps": ("viridis", False, "SEEPS (0 = perfect)"),
    "seepsdiff": ("RdBu_r", True, "SEEPS difference vs reference (red = worse)"),
}
PERCENT = {"rmserel", "fcrpsrel"}
UNITLESS = {"seeps", "seepsdiff"}
NAMES = {"mef": "MEF", "rmserel": "RMSE change", "fcrps": "FCRPS",
         "fcrpsdiff": "FCRPS diff", "fcrpsrel": "FCRPS change"}
OROG = ("BrBG_r", "model minus station elevation")  # step-independent map orog_<label>.png  # metrics plotted in % rather than the parameter's units
STYLE = 3  # bump to redraw all images after a change in the plotting code

LAND, OCEAN = "#f3f1ec", "#e3ecf3"
SERIES = ["#1f4e79", "#d1495b", "#2e8b57", "#e09f3e", "#6a4c93", "#3d9ec9"]


def map_note(m, val):
    """Short summary written in the corner of a station map."""
    v = np.asarray(val, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return ""
    unit = " %" if m in PERCENT else ""
    fmt = "+.2f" if METRICS[m][1] else ".2f"
    note = f"mean {v.mean():{fmt}}{unit}   median {np.median(v):{fmt}}{unit}"
    if m.endswith("diff") or m in PERCENT:
        note += f"\nbetter at {100 * (v < 0).mean():.0f}% / worse at {100 * (v > 0).mean():.0f}% of stations"
    return note


def draw_map(path, lon, lat, val, title, cmap, vmin, vmax, extent, label, dpi, subtitle="", note=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    from matplotlib.ticker import MaxNLocator
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    lon, lat, val = map(np.asarray, (lon, lat, val))
    ok = np.isfinite(val)
    lon, lat, val = lon[ok], lat[ok], val[ok]
    w, e, s, n = extent or (-180, 180, -90, 90)
    area = (e - w) * (n - s)
    data_crs = ccrs.PlateCarree()
    if extent:
        proj = ccrs.PlateCarree(central_longitude=0.5 * (w + e))
        aspect = (n - s) / max(e - w, 1)
    else:
        proj = ccrs.Robinson()
        aspect = 0.5
    fig = plt.figure(figsize=(11, max(4.0, min(10.0, 11 * aspect + 1.6))), layout="constrained")
    ax = fig.add_subplot(projection=proj)
    if extent:
        ax.set_extent(extent, crs=data_crs)
    else:
        ax.set_global()
    res = "50m" if area < 5000 else "110m"
    try:
        ax.add_feature(cfeature.OCEAN.with_scale(res), facecolor=OCEAN, edgecolor="none", zorder=0)
        ax.add_feature(cfeature.LAND.with_scale(res), facecolor=LAND, edgecolor="none", zorder=0)
        ax.add_feature(cfeature.BORDERS.with_scale(res), linewidth=0.3, edgecolor="0.6", zorder=1)
        ax.coastlines(resolution=res, linewidth=0.45, color="0.35", zorder=1)
    except Exception as exc:  # no Natural Earth data available offline
        log("coastlines unavailable:", exc)
    gl = ax.gridlines(crs=data_crs, draw_labels=True, linewidth=0.4, color="0.6",
                      alpha=0.5, linestyle=":", zorder=1)
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"size": 7, "color": "0.35"}
    ax.spines["geo"].set_edgecolor("0.4")
    ax.spines["geo"].set_linewidth(0.6)

    levels = MaxNLocator(nbins=12, symmetric=vmin < 0).tick_values(vmin, vmax)
    cm = plt.get_cmap(cmap)
    norm = BoundaryNorm(levels, cm.N, extend="both")
    order = np.argsort(np.abs(val) if vmin < 0 else val)  # largest values on top
    size = float(np.clip(5 * np.sqrt(64800.0 / max(area, 1)), 6, 40))
    sc = ax.scatter(lon[order], lat[order], c=val[order], s=size, cmap=cm, norm=norm,
                    edgecolors="0.25", linewidths=0.15, transform=data_crs, zorder=2)
    cb = fig.colorbar(sc, ax=ax, orientation="horizontal", extend="both", shrink=0.6,
                      aspect=45, pad=0.02, ticks=levels[::2] if len(levels) > 9 else levels)
    cb.set_label(label, fontsize=9)
    cb.ax.tick_params(labelsize=8, length=2)
    cb.outline.set_linewidth(0.5)

    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
    if subtitle:
        ax.set_title(subtitle, loc="right", fontsize=8, color="0.35")
    if note:
        ax.text(0.01, 0.015, note, transform=ax.transAxes, fontsize=7.5, va="bottom", ha="left",
                zorder=3, bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.7", lw=0.5, alpha=0.9))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.08, facecolor="white")
    plt.close(fig)


def draw_scores(path, p, steps, scores, labels, units, region, dpi, ensemble=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ref = labels[0]
    seeps = any("seeps" in sc[ref] for sc in scores.values())
    panels = [("bias", f"Bias [{units}]"), ("rmse", f"RMSE [{units}]"), ("mae", f"MAE [{units}]")]
    if len(labels) > 1:
        panels.append(("mef", f"Mean bias difference vs {ref} [{units}]"))
    if ensemble:
        panels.append(("fcrps", f"FCRPS [{units}]"))
    if len(labels) > 1:
        panels.append(("rmserel", f"RMSE change vs {ref} [%]"))
        if ensemble:
            panels.append(("fcrpsrel", f"FCRPS change vs {ref} [%]"))
    if seeps:
        panels.append(("seeps", "SEEPS [-]"))
        if len(labels) > 1:
            panels.append(("seepsdiff", f"SEEPS difference vs {ref} [-]"))
    ncols = 3 if len(panels) > 4 else len(panels)
    nrows = -(-len(panels) // ncols)
    fig, axs = plt.subplots(nrows, ncols, figsize=(4.4 * ncols, 4.0 * nrows), sharex=True,
                            squeeze=False, layout="constrained")
    axs = axs.ravel()
    for ax in axs[len(panels):]:
        ax.set_visible(False)
    colors = {lab: SERIES[i % len(SERIES)] for i, lab in enumerate(labels)}

    def get(s, lab, m):
        return scores[s][lab].get(m, np.nan) if s in scores else np.nan

    for ax, (m, ylabel) in zip(axs, panels):
        for lab in labels:
            if m.endswith("rel"):
                if lab == ref:
                    continue
                y = [100.0 * (get(s, lab, m[:-3]) / get(s, ref, m[:-3]) - 1.0) for s in steps]
            elif m == "mef":
                if lab == ref:
                    continue
                y = [get(s, lab, "bias") - get(s, ref, "bias") for s in steps]
            elif m.endswith("diff"):
                if lab == ref:
                    continue
                y = [get(s, lab, m[:-4]) - get(s, ref, m[:-4]) for s in steps]
            else:
                y = [get(s, lab, m) for s in steps]
            ax.plot(steps, y, marker="o", ms=3.5, lw=1.6, color=colors[lab], label=lab)
        if m in ("bias", "mef") or m.endswith("rel") or m.endswith("diff"):
            ax.axhline(0, color="0.3", lw=0.8)
        if m.endswith("rel"):
            lo, hi = ax.get_ylim()
            ax.axhspan(0, max(hi, 0), color="#d1495b", alpha=0.06, lw=0)
            ax.axhspan(min(lo, 0), 0, color="#2e8b57", alpha=0.06, lw=0)
            ax.set_ylim(lo, hi)
        ax.set_title(ylabel.split(" [")[0], loc="left", fontsize=11, fontweight="bold")
        ax.set_xlabel("forecast step [h]")
        ax.set_ylabel(ylabel)
        ax.set_xticks(steps if len(steps) <= 12 else steps[::max(1, len(steps) // 10)])
        ax.grid(alpha=0.3, lw=0.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    handles, names = axs[0].get_legend_handles_labels()
    fig.legend(handles, names, loc="outside lower center", ncol=len(labels), frameon=False)
    fig.suptitle(f"{p.label} – {region}", fontsize=12, fontweight="bold", x=0.01, ha="left")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, facecolor="white")
    plt.close(fig)


def plot_task(job):
    kind, path, kw = job
    if kind == "map":
        draw_map(Path(path), **kw)
    else:
        draw_scores(Path(path), **kw)
    return path


def stage_plot(cfg, args):
    from multiprocessing import Pool

    with lock(cfg.workdir / "locks" / "plot.lock"):
        labels = [e["label"] for e in cfg.experiments]
        ref = labels[0]
        nproc = int(os.environ.get("SLURM_CPUS_PER_TASK", cfg.plot["processes"]))
        dpi = int(cfg.plot["dpi"])
        pperiods = cfg.plot_periods()
        periods = list(pperiods)
        tasks = [(str(cfg.path), per, pi, s) for per in periods
                 for pi, p in enumerate(cfg.params) for s in cfg.steps_for(p)]
        if any(p.seeps for p in cfg.params):
            seeps_climatology(cfg)  # fetched once here, not by each worker
        log(f"computing statistics for {len(tasks)} period/parameter/step combinations")
        results = {per: defaultdict(dict) for per in periods}  # [period][param index][step]
        with Pool(nproc) as pool:
            for per, pi, step, res in pool.imap_unordered(compute_step, tasks):
                if res is not None:
                    results[per][pi][step] = res

        sig_file = cfg.plot_dir / "signatures.json"
        old_sigs = json.loads(sig_file.read_text()) if sig_file.exists() and not args.force else {}
        sigs, jobs = {}, []
        manifest = {"title": cfg.name, "experiments": labels, "reference": ref,
                    "regions": list(cfg.plot["regions"]), "metrics": {}, "periods": {},
                    "generated": dt.datetime.now().strftime("%Y-%m-%d %H:%M")}
        for per in periods:
            d, hour, classes = pperiods[per]
            dates = (f"{ymd(d[0])}-{ymd(d[-1])}" + ("" if hour is None else f", valid {hour:02d} UTC only")
                     + ("" if classes is None else f", stations with {'/'.join(classes)} identifiers only"))
            manifest["periods"][per] = {"dates": dates, "params": {}}
        for m, (_, div, desc) in cfg.metrics.items():
            if cfg.ensemble and m in ("bias", "rmse", "mae"):
                desc += " of the ensemble mean for ensembles"
            elif m == "fcrps":
                desc += " (absolute error for single forecasts)"
            manifest["metrics"][m] = desc
        manifest["metrics"]["orog"] = (OROG[1] + " (brown = model surface above the station)."
                                       " Stations verified at any step are shown.")

        def want(path, sig):
            sigs[str(path)] = sig
            return args.force or old_sigs.get(str(path)) != sig or not path.exists()

        for pi, p in enumerate(cfg.params):
            if not any(results[per][pi] for per in periods):
                log(f"{p.key}: no matched forecast/observation pairs, skipped")
                continue
            for rname, extent in cfg.plot["regions"].items():
                found = {per: sorted(s for s, r in results[per][pi].items() if rname in r["regions"])
                         for per in periods}
                all_stats = [results[per][pi][s]["regions"][rname]["stations"]
                             for per in periods for s in found[per]]
                if not all_stats:
                    continue
                # common colour limits over all steps and periods so that maps are comparable
                lims = {}
                for m in cfg.metrics:
                    cols = [c for c in all_stats[0].columns if c.startswith(m + "_")]
                    vals = np.concatenate([st[c].values for st in all_stats for c in cols]) if cols else []
                    lims[m] = nice_limit(vals, 98 if METRICS[m][1] else 95)
                    if m in PERCENT:
                        lims[m] = min(lims[m], 100.0)
                dzcols = [f"dz_{lab}" for lab in labels]
                lims["orog"] = nice_limit(np.concatenate(
                    [st[c].values for st in all_stats for c in dzcols]), 98)
                for per in periods:
                    rsteps = found[per]
                    if not rsteps:
                        continue
                    res = results[per][pi]
                    units = res[rsteps[0]]["units"]
                    rdir = cfg.plot_dir / per / p.key / rname
                    stats = {s: res[s]["regions"][rname]["stations"] for s in rsteps}
                    scores = {s: res[s]["regions"][rname]["scores"] for s in rsteps}
                    ptag = "" if per == "all" else f"{per} · "
                    for s in rsteps:
                        st = stats[s]
                        bases = res[s]["bases"]
                        span = f"{bases[0][:8]}-{bases[-1][:8]}"
                        for m, (cmap, div, desc) in cfg.metrics.items():
                            for lab in labels:
                                col = f"{m}_{lab}"
                                if col not in st:
                                    continue
                                lim = lims[m]
                                vmin, vmax = (-lim, lim) if div else (0.0, lim)
                                path = rdir / f"{m}_{lab}_T{s:03d}.png"
                                vs = (f" − {ref}" if m.endswith("diff") or m == "mef"
                                      else f" vs {ref}" if m in PERCENT else "")
                                name = NAMES.get(m, m.upper().replace("DIFF", " diff"))
                                title = f"{p.label} {name}   {lab}{vs}   T+{s}h"
                                nst = int(np.isfinite(st[col]).sum())
                                subtitle = f"{ptag}{rname} · {len(bases)} forecasts {span} · {nst} stations"
                                note = map_note(m, st[col].values)
                                sig = hashlib.md5(json.dumps(
                                    [STYLE, dpi, bases, len(st), float(st[col].sum()), vmin, vmax,
                                     title, subtitle, note]).encode()).hexdigest()
                                if want(path, sig):
                                    jobs.append(("map", str(path), dict(
                                        lon=st.longitude.values, lat=st.latitude.values,
                                        val=st[col].values, title=title, subtitle=subtitle, note=note,
                                        cmap=cmap, vmin=vmin, vmax=vmax, extent=extent,
                                        label=f"{desc} [{'%' if m in PERCENT else '-' if m in UNITLESS else units}]",
                                        dpi=dpi)))
                        if cfg.plot["daily_maps"]:
                            jobs += daily_jobs(p, s, rname, extent, res[s]["regions"][rname]["daily"],
                                               units, dpi, want, labels, rdir / "daily")
                    # model - station height, over the stations verified at any step
                    dz = pd.concat([stats[s][["key", "longitude", "latitude"] + dzcols] for s in rsteps])
                    dz = dz.drop_duplicates("key")
                    nbases = len({b for s in rsteps for b in res[s]["bases"]})
                    orog = []
                    for lab in labels:
                        val = dz[f"dz_{lab}"].values
                        ok = np.isfinite(val)
                        if not ok.any():
                            continue
                        orog.append(lab)
                        path = rdir / f"orog_{lab}.png"
                        lim = lims["orog"]
                        title = f"Model − station height   {lab}"
                        subtitle = f"{ptag}{rname} · {p.label} stations · {ok.sum()} stations"
                        v = val[ok]
                        note = (f"mean {v.mean():+.0f} m   median {np.median(v):+.0f} m   "
                                f"mean |dz| {np.abs(v).mean():.0f} m\n"
                                f"model above station at {100 * (v > 0).mean():.0f}% / "
                                f"below at {100 * (v < 0).mean():.0f}% of stations")
                        sig = hashlib.md5(json.dumps(
                            [STYLE, dpi, nbases, int(ok.sum()), float(v.sum()), lim,
                             title, subtitle, note]).encode()).hexdigest()
                        if want(path, sig):
                            jobs.append(("map", str(path), dict(
                                lon=dz.longitude.values, lat=dz.latitude.values, val=val,
                                title=title, subtitle=subtitle, note=note, cmap=OROG[0],
                                vmin=-lim, vmax=lim, extent=extent,
                                label=f"{OROG[1]} [m]", dpi=dpi)))
                    path = rdir / "scores.png"
                    region_title = rname if per == "all" else f"{rname} – {per}"
                    if str(cfg.plot.get("spatial_mean_weights") or "none").lower() == "station_density":
                        region_title += " (station-density weighted)"
                    sig = hashlib.md5(json.dumps([STYLE, dpi, region_title, scores],
                                                 sort_keys=True, default=str).encode()).hexdigest()
                    if want(path, sig):
                        jobs.append(("scores", str(path), dict(p=p, steps=rsteps, scores=scores,
                                                               labels=labels, units=units,
                                                               region=region_title, dpi=dpi,
                                                               ensemble=cfg.ensemble)))
                    pinfo = manifest["periods"][per]["params"].setdefault(
                        p.key, {"label": p.label, "units": units, "regions": {}})
                    pinfo["regions"][rname] = {
                        "steps": rsteps,
                        "orog": orog,
                        "scores": {str(s): scores[s] for s in rsteps},
                        "daily": sorted({b for s in rsteps for b in res[s]["bases"]})
                        if cfg.plot["daily_maps"] else [],
                    }

        log(f"{len(jobs)} images to draw ({len(sigs) - len(jobs)} unchanged, skipped)")
        with Pool(nproc) as pool:
            for i, _ in enumerate(pool.imap_unordered(plot_task, jobs, chunksize=4), 1):
                if i % 50 == 0 or i == len(jobs):
                    log(f"  {i}/{len(jobs)} images")
        sig_file.parent.mkdir(parents=True, exist_ok=True)
        sig_file.write_text(json.dumps(sigs))
        write_html(cfg, manifest)
    if cfg.publish["site"]:
        return stage_publish(cfg, args)
    return 0


def daily_jobs(p, s, rname, extent, df, units, dpi, want, labels, rdir):
    """Maps of the observed value and each experiment's error per forecast."""
    jobs = []
    lo, hi = np.nanpercentile(df["obs"], [2, 98]) if len(df) else (0, 1)
    elim = nice_limit(np.concatenate([error(p, df[lab], df["obs"]).values for lab in labels]), 98)
    for base, g in df.groupby("base"):
        vt = dt.datetime.strptime(base, "%Y%m%d%H") + dt.timedelta(hours=s)
        items = [("obs", g["obs"], "viridis", lo, hi, f"observed value [{units}]")]
        items += [(f"err_{lab}", error(p, g[lab], g["obs"]), "RdBu_r",
                   -elim, elim, f"error fc - obs [{units}]") for lab in labels]
        for name, val, cmap, vmin, vmax, clabel in items:
            path = rdir / base / f"{name}_T{s:03d}.png"
            who = "observations" if name == "obs" else name[4:] + " error"
            title = f"{p.label} {who}   T+{s}h"
            subtitle = (f"{rname} · base {base[:8]} {base[8:]}z · valid {vt:%Y-%m-%d %H}z"
                        f" · {len(g)} stations")
            note = "" if name == "obs" else map_note("bias", val.values)
            sig = hashlib.md5(json.dumps([STYLE, dpi, len(g), float(val.sum()), vmin, vmax,
                                          title, subtitle]).encode()).hexdigest()
            if want(path, sig):
                jobs.append(("map", str(path), dict(lon=g.longitude.values, lat=g.latitude.values,
                                                    val=val.values, title=title, subtitle=subtitle,
                                                    note=note, cmap=cmap,
                                                    vmin=float(vmin), vmax=float(vmax), extent=extent,
                                                    label=clabel, dpi=dpi)))
    return jobs


def write_html(cfg, manifest):
    html = (Path(__file__).with_name("browser_template.html")).read_text()
    html = html.replace("__MANIFEST__", json.dumps(manifest)).replace("__TITLE__", cfg.name)
    out = cfg.workdir / "index.html"
    out.write_text(html)
    log(f"browser written: {out}")


# --------------------------------------------------------------------------
# publish on sites.ecmwf.int
# --------------------------------------------------------------------------

def put_file(site, token, local, remote):
    """PUT one file through the sites.ecmwf.int files API; True on success.
    The token is passed on stdin so that it does not show up in `ps`."""
    url = f"{site}s/api/v1/files/{remote}"
    r = subprocess.run(["/usr/bin/curl", "--silent", "--show-error", "--fail", "--retry", "3",
                        "-K", "-", "-X", "PUT", "-F", f"file=@{local}", url],
                       input=f'header = "Authorization: Bearer {token}"\n',
                       text=True, capture_output=True)
    if r.returncode:
        log(f"upload of {remote} failed: {r.stderr.strip()}")
    return r.returncode == 0


def stage_publish(cfg, args):
    """Upload index.html and the plots to the site; unchanged files are skipped."""
    site = cfg.publish["site"]
    if not site:
        log("publish.site is not set in the configuration, nothing to publish")
        return 1
    site = site.rstrip("/") + "/"
    token_file = Path(os.path.expanduser(cfg.publish["token_file"]))
    if not token_file.exists():
        log(f"token file {token_file} not found (copy the token from the site's Settings tab)")
        return 1
    token = token_file.read_text().strip()
    prefix = (cfg.publish["path"] or cfg.name).strip("/") + "/"
    record = cfg.workdir / "published.json"
    old = json.loads(record.read_text()) if record.exists() and not args.force else {}
    if old.get("_site") != site + prefix:
        old = {}

    sig_file = cfg.plot_dir / "signatures.json"
    current = json.loads(sig_file.read_text()) if sig_file.exists() else {}
    files = sorted(Path(f) for f in current if Path(f).exists())
    todo = []
    for f in files + [cfg.workdir / "index.html"]:
        rel = f.relative_to(cfg.workdir).as_posix()
        digest = hashlib.md5(f.read_bytes()).hexdigest()
        if old.get(rel) != digest:
            todo.append((f, rel, digest))
    log(f"publishing to {site}{prefix}: {len(todo)} files to upload "
        f"({len(files) + 1 - len(todo)} unchanged)")
    new = {k: v for k, v in old.items()}
    new["_site"] = site + prefix

    def one(item):
        f, rel, digest = item
        return rel, digest, put_file(site, token, f, prefix + rel)

    # images first, index.html last so that the page never points at missing plots
    images = [x for x in todo if x[1] != "index.html"]
    index = [x for x in todo if x[1] == "index.html"]
    pcfg = argparse.Namespace(parallel=max(1, int(cfg.publish["parallel"])))
    failed = 0
    for batch in (images, index):
        for rel, digest, ok in parallel_map(pcfg, one, batch):
            if ok:
                new[rel] = digest
            else:
                failed += 1
        record.write_text(json.dumps(new))
    if failed:
        log(f"{failed} uploads failed; rerun 'publish' to retry")
        return 1
    log(f"published: {site}{prefix}index.html")
    return 0


# --------------------------------------------------------------------------
# slurm submission
# --------------------------------------------------------------------------

JOB_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --qos={qos}
#SBATCH --time={time}
#SBATCH --mem={mem}
#SBATCH --cpus-per-task={cpus}
#SBATCH --output={logdir}/%x-%j.out
{extra}
set -e
module load python3
cd {workdir}
python3 {script} {stage} {config} {args}
"""


def write_job(cfg, name, stage, args=""):
    sl = cfg.slurm
    res = sl[stage]
    jobdir = cfg.workdir / "jobs"
    logdir = cfg.workdir / "logs"
    jobdir.mkdir(parents=True, exist_ok=True)
    logdir.mkdir(parents=True, exist_ok=True)
    extra = f"#SBATCH --account={sl['account']}" if sl.get("account") else ""
    text = JOB_TEMPLATE.format(name=name, qos=sl["qos"], time=res["time"], mem=res["mem"],
                               cpus=res["cpus"], logdir=logdir, extra=extra,
                               workdir=cfg.workdir, script=SCRIPT, stage=stage,
                               config=cfg.path, args=args)
    path = jobdir / f"{name}.sh"
    path.write_text(text)
    return path


def sbatch(path, deps, dry):
    cmd = ["sbatch", "--parsable"]
    if deps:
        cmd.append("--dependency=afterany:" + ":".join(deps))
    cmd.append(str(path))
    if dry:
        print("DRY-RUN:", " ".join(cmd))
        return f"<{path.stem}>"
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip()
    jobid = out.split(";")[0]
    log(f"submitted {path.name} as job {jobid}")
    return jobid


def stage_submit(cfg, args):
    stages = args.stages
    force = " --force" if args.force else ""
    upload_ids, collocate_ids, download_ids = [], [], []
    if "upload" in stages:
        for exp in cfg.experiments:
            if any(uploads(s) for s, _ in exp["variants"]):
                job = write_job(cfg, f"{cfg.name}_upload_{exp['label']}", "upload",
                                f"--experiment {exp['label']}{force}")
                upload_ids.append(sbatch(job, [], args.dry_run))
    if "collocate" in stages:
        for exp in cfg.experiments:
            if collocates(cfg, exp):
                job = write_job(cfg, f"{cfg.name}_collocate_{exp['label']}", "collocate",
                                f"--experiment {exp['label']}{force}")
                collocate_ids.append(sbatch(job, [], args.dry_run))
    if "download" in stages:
        job = write_job(cfg, f"{cfg.name}_download", "download")
        download_ids.append(sbatch(job, upload_ids, args.dry_run))
    if "plot" in stages:
        job = write_job(cfg, f"{cfg.name}_plot", "plot", force.strip())
        sbatch(job, (download_ids or upload_ids) + collocate_ids, args.dry_run)
    log(f"job logs: {cfg.workdir / 'logs'};  monitor with: squeue -u $USER")
    return 0


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="stage", required=True)
    for name in ("submit", "upload", "collocate", "download", "plot", "publish", "status"):
        sp = sub.add_parser(name)
        sp.add_argument("config")
        if name in ("upload", "collocate", "plot", "publish", "submit"):
            sp.add_argument("--force", action="store_true",
                            help="upload: re-upload even if in VINO; collocate: recompute even "
                                 "if in the store; plot: redraw all images; "
                                 "publish: re-upload all files")
        if name in ("upload", "collocate"):
            sp.add_argument("--experiment", nargs="+", help="only these experiment labels")
        if name == "submit":
            sp.add_argument("--stages", nargs="+", default=["upload", "collocate", "download", "plot"],
                            choices=["upload", "collocate", "download", "plot"])
            sp.add_argument("--dry-run", action="store_true", help="write job scripts, do not submit")
    args = ap.parse_args()
    # submit checks MARS again for periods ending in "-"; its jobs reuse that answer
    cfg = Config(args.config, refresh_open=args.stage == "submit")
    cfg.workdir.mkdir(parents=True, exist_ok=True)
    stage = {"submit": stage_submit, "upload": stage_upload, "collocate": stage_collocate,
             "download": stage_download,
             "plot": stage_plot, "publish": stage_publish, "status": stage_status}[args.stage]
    sys.exit(stage(cfg, args))


if __name__ == "__main__":
    main()
