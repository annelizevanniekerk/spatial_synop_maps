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
                             TEMP observations with vino_upper.py
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
UPPER_ENV = "module load python3 ecmwf-toolbox vtb/vino >/dev/null 2>&1"

DEFAULTS = {
    "name": "synop_maps",
    "workdir": None,
    "store": None,  # shared geopoint store; default: ./store next to this script
    "collocate": "vino",  # vino: upload to + download from VINO; local: genfc -> files
    "times": [0],
    "sources": ["synop"],
    "interpolation": "nearest_point",
    "chunk_size": 15,
    "parallel": 4,  # concurrent vino_genfc / vino_getgeo calls
    "tools": {
        "genfc": "~moz/bin/vino_genfc",
        "getgeo": "~moz/bin/vino_getgeo",
        "mars": "mars",
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
    "database": None,  # alternative MARS database, e.g. fdb for data not archived yet
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


def period_spec(e, per):
    """Experiment settings for one period (values given as {period: value} resolved)."""
    return {k: ({str(a): b for a, b in v.items()}.get(per) if isinstance(v, dict) else v)
            for k, v in e.items()}


def mars_last_date(cfg, e, start):
    """Last base date >= start for which MARS has the longest configured step
    of experiment e at every configured base time (2t is used as the probe),
    and whose observations for that step already exist (valid time before
    yesterday). None if there is nothing."""
    if e.get("database"):
        sys.exit(f"experiment {e['expver']}: open end dates ('-') are not supported with "
                 f"`database: {e['database']}` (MARS cannot list it), give an end date")
    maxstep = max(cfg.steps)
    cap = dt.date.today() - dt.timedelta(days=maxstep // 24 + 2)
    if cap < start:
        return None
    last = []
    for t in cfg.times:
        req = {"class": e["class"], "stream": e["stream"], "expver": e["expver"],
               "type": e["type"], "levtype": "sfc", "param": "167", "time": f"{t:02d}",
               "step": maxstep, "date": f"{ymd(start)}/to/{ymd(cap)}", "output": "tree"}
        for k in ("anoffset", "origin"):
            if e.get(k) not in (None, ""):
                req[k] = e[k]
        text = "list," + ",".join(f"{k}={v}" for k, v in req.items()) + "\n"
        log(f"MARS list: {text.strip()}")
        r = subprocess.run([cfg.tools["mars"]], input=text, capture_output=True, text=True)
        dates = sorted({dt.date(int(y), int(m), int(d))
                        for y, m, d in re.findall(r"\b(\d{4})-(\d{2})-(\d{2})\b", r.stdout)})
        if r.returncode:
            log("MARS list failed:", (r.stderr or r.stdout).strip()[-500:])
        # the data must start near the period start and is followed only up to
        # the first gap of more than MAX_GAP days, so that e.g. a run of the
        # same expver in another year does not extend the period
        if r.returncode or not dates or (dates[0] - start).days > MAX_GAP:
            log(f"{e.get('label', e['expver'])} ({e['expver']}) {t:02d}z: no forecasts in MARS from "
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
        self.interpolation = c["interpolation"]
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
        self.metrics = {m: v for m, v in METRICS.items() if self.ensemble or not m.startswith("fcrps")}

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
                log(f"period {per}: {e.get('label', e['expver'])} ({e['expver']}) is in MARS up to "
                    + (ymd(d) if d else "nothing after " + ymd(start)))
            if None in last:
                sys.exit(f"period {per}: an experiment has no forecasts in MARS from {ymd(start)}")
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
    def obs_dir(self, p):
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
                      f"{self.upper['interpolation']}_{'_'.join(self.upper['sources'])}"]
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
                    if p.upper:  # not in VINO; collocated locally
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
    """Upload one chunk of dates and return the dates still missing from VINO."""
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
    if missing:
        log(f"[{exp['label']} {p.key} {t:02d}z] WARNING rc={rc}, not in VINO after upload: "
            + " ".join(map(ymd, missing)))
    return [(exp["label"], p.key, t, ymd(d)) for d in missing]


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
    return spec["collocate"] == "local" or (p is not None and p.upper)


def collocates(cfg, e):
    return any(collocated(s) for s, _ in e["variants"]) or any(p.upper for p in cfg.params)


def stage_collocate(cfg, args):
    exps = [e for e in cfg.experiments if collocates(cfg, e)]
    if args.experiment:
        exps = [e for e in exps if e["label"] in args.experiment]
    if not exps:
        log("nothing to collocate (no experiments with collocate: local selected)")
        return 0
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
    to_store = upper_to_store if p.upper else genfc_to_store
    if p.upper:
        # the forecasts are picked at the stations of the observations, so
        # these must be in the store first
        needed = defaultdict(set)
        for d in dates:
            for s in steps:
                vt = valid_time(d, t, s)
                needed[vt.hour].add(vt.date())
        for hour, ds in sorted(needed.items()):
            download_obs(cfg, p, hour, ds)

    def complete(d):
        # pressure levels: only steps with observations can be collocated
        return all(fc_file(root, p, d, t, s) for s in steps
                   if not p.upper or obs_file(cfg, p, valid_time(d, t, s)).exists())

    with cfg.store_lock(root, f"{p.key}_{t:02d}"):
        todo = list(dates) if force else [d for d in dates if not complete(d)]
        log(f"[{spec['label']} {spec['expver']} {p.key} {t:02d}z] {len(dates) - len(todo)} dates "
            f"already in the store, {len(todo)} to collocate")
        failed = []
        for chunk in chunks(todo, cfg.chunk_size):
            if to_store(cfg, spec, p, t, chunk, steps) == 0:
                continue
            # vino_genfc stops at the first error; redo the dates one at a time
            if len(chunk) > 1:
                log(f"[{spec['label']} {p.key} {t:02d}z] retrying "
                    f"{ymd(chunk[0])}-{ymd(chunk[-1])} one date at a time")
                failed += [d for d in chunk if to_store(cfg, spec, p, t, [d], steps)]
            else:
                failed += chunk
        missing = [d for d in todo if d in failed or not complete(d)]
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
    """Collocate a pressure-level parameter with vino_upper.py into the store."""
    a = ["fc", "--parameter", p.param, "--level", p.level, "--key", p.key,
         "--dates", *map(ymd, dates), "--time", t, "--steps", *steps,
         "--class", spec["class"], "--stream", spec["stream"], "--expver", spec["expver"],
         "--type", spec["type"], "--grid", spec["grid"] or cfg.upper["grid"],
         "--interpolation", cfg.upper["interpolation"], "--obsdir", cfg.obs_dir(p),
         "--mars", cfg.tools["mars"]]
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
        for p in cfg.params:
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
    dropped. This follows the processing of quaver's error_map. The columns
    dz_<label> hold the model minus station elevation (NaN when unknown).
    For ensembles (several members per file) the forecast column is the
    ensemble mean, and crps_<label> holds the fair CRPS of each pair (the
    absolute error for single forecasts).
    """
    frames, units = [], ""
    lapse = float(cfg.plot["lapse_rate"] or 0) if p.param in LAPSE_PARAMS else 0.0
    max_dz = cfg.plot["max_dz"]
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
    ref = labels[0]
    for lab in labels[1:]:
        st[f"rmsediff_{lab}"] = st[f"rmse_{lab}"] - st[f"rmse_{ref}"]
        st[f"maediff_{lab}"] = st[f"mae_{lab}"] - st[f"mae_{ref}"]
        rel = 100.0 * (st[f"rmse_{lab}"] / st[f"rmse_{ref}"] - 1.0)
        st[f"rmserel_{lab}"] = rel.replace([np.inf, -np.inf], np.nan)
        st[f"fcrpsdiff_{lab}"] = st[f"fcrps_{lab}"] - st[f"fcrps_{ref}"]
        rel = 100.0 * (st[f"fcrps_{lab}"] / st[f"fcrps_{ref}"] - 1.0)
        st[f"fcrpsrel_{lab}"] = rel.replace([np.inf, -np.inf], np.nan)
    min_count = max(1, float(cfg.plot["min_fraction"]) * nbases)
    return st[st["n"] >= min_count].reset_index()


def overall_scores(cfg, st):
    """Regional scores as unweighted means over the station scores, like the
    domain means of quaver's error_map: RMSE = sqrt(mean of station RMSE^2)."""
    out = {"n": int(st["n"].sum()), "stations": int(len(st))}
    for e in cfg.experiments:
        lab = e["label"]
        out[lab] = {"bias": float(st[f"bias_{lab}"].mean()),
                    "rmse": float(np.sqrt((st[f"rmse_{lab}"] ** 2).mean())),
                    "mae": float(st[f"mae_{lab}"].mean()),
                    "fcrps": float(st[f"fcrps_{lab}"].mean())}
    return out


def compute_step(args):
    cfg_path, per, pi, step = args
    cfg = Config(cfg_path)
    p = cfg.params[pi]
    df, units = matched(cfg, p, step, cfg.periods[per])
    if df is None or df.empty:
        return per, pi, step, None
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
    "rmsediff": ("RdBu_r", True, "RMSE difference vs reference (red = worse)"),
    "maediff": ("RdBu_r", True, "MAE difference vs reference (red = worse)"),
    "rmserel": ("RdBu_r", True, "relative RMSE change vs reference (red = worse)"),
    # only with ensemble experiments; for single forecasts the CRPS is the MAE
    "fcrps": ("viridis", False, "fair CRPS"),
    "fcrpsdiff": ("RdBu_r", True, "fair CRPS difference vs reference (red = worse)"),
    "fcrpsrel": ("RdBu_r", True, "relative fair CRPS change vs reference (red = worse)"),
}
PERCENT = {"rmserel", "fcrpsrel"}
NAMES = {"rmserel": "RMSE change", "fcrps": "FCRPS", "fcrpsdiff": "FCRPS diff",
         "fcrpsrel": "FCRPS change"}
OROG = ("BrBG_r", "model minus station elevation")  # step-independent map orog_<label>.png  # metrics plotted in % rather than the parameter's units
STYLE = 2  # bump to redraw all images after a change in the plotting code

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
    panels = [("bias", f"Bias [{units}]"), ("rmse", f"RMSE [{units}]"), ("mae", f"MAE [{units}]")]
    if ensemble:
        panels.append(("fcrps", f"FCRPS [{units}]"))
    if len(labels) > 1:
        panels.append(("rmserel", f"RMSE change vs {ref} [%]"))
        if ensemble:
            panels.append(("fcrpsrel", f"FCRPS change vs {ref} [%]"))
    ncols = 3 if len(panels) > 4 else len(panels)
    nrows = -(-len(panels) // ncols)
    fig, axs = plt.subplots(nrows, ncols, figsize=(4.4 * ncols, 4.0 * nrows), sharex=True,
                            squeeze=False, layout="constrained")
    axs = axs.ravel()
    for ax in axs[len(panels):]:
        ax.set_visible(False)
    colors = {lab: SERIES[i % len(SERIES)] for i, lab in enumerate(labels)}

    def get(s, lab, m):
        return scores[s][lab][m] if s in scores else np.nan

    for ax, (m, ylabel) in zip(axs, panels):
        for lab in labels:
            if m.endswith("rel"):
                if lab == ref:
                    continue
                y = [100.0 * (get(s, lab, m[:-3]) / get(s, ref, m[:-3]) - 1.0) for s in steps]
            else:
                y = [get(s, lab, m) for s in steps]
            ax.plot(steps, y, marker="o", ms=3.5, lw=1.6, color=colors[lab], label=lab)
        if m == "bias" or m.endswith("rel"):
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
        periods = list(cfg.periods)
        tasks = [(str(cfg.path), per, pi, s) for per in periods
                 for pi, p in enumerate(cfg.params) for s in cfg.steps_for(p)]
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
            d = cfg.periods[per]
            manifest["periods"][per] = {"dates": f"{ymd(d[0])}-{ymd(d[-1])}", "params": {}}
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
                                vs = f" − {ref}" if m.endswith("diff") else f" vs {ref}" if m in PERCENT else ""
                                name = NAMES.get(m, m.upper().replace("DIFF", " diff"))
                                title = f"{p.label} {name}   {lab}{vs}   T+{s}h"
                                subtitle = f"{ptag}{rname} · {len(bases)} forecasts {span} · {len(st)} stations"
                                note = map_note(m, st[col].values)
                                sig = hashlib.md5(json.dumps(
                                    [STYLE, dpi, bases, len(st), float(st[col].sum()), vmin, vmax,
                                     title, subtitle, note]).encode()).hexdigest()
                                if want(path, sig):
                                    jobs.append(("map", str(path), dict(
                                        lon=st.longitude.values, lat=st.latitude.values,
                                        val=st[col].values, title=title, subtitle=subtitle, note=note,
                                        cmap=cmap, vmin=vmin, vmax=vmax, extent=extent,
                                        label=f"{desc} [{'%' if m in PERCENT else units}]", dpi=dpi)))
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
