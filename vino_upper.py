#!/usr/bin/env python3
"""
vino_upper.py - upper-air (TEMP) helper for synop_maps.py.

vino_getgeo cannot select a pressure level and vino_genfc has no upper-air
parameters, so this script does both jobs for pressure-level parameters:

    obs  TEMP observations of one parameter/level from VINO -> geopoint files
    fc   MARS pressure-level forecasts -> values at the stations of the
         observation file valid at each step -> geopoint files

`fc` without --level does the same for surface parameters (2t, 2d, 10ff, ...),
e.g. at the stations of observations from MARS feedback, and also writes the
model orography at the stations (elevation column, for the height correction).

Files are named and laid out like those of vino_getgeo
(<outdir>/YYYY/MM/DD/<key>_obs_YYYYMMDDHH.geo and
<key>_mod_YYYYMMDDHH_<step>_<stream>_<expver>_<type>.geo, with key e.g. t850).
synop_maps.py runs it with the modules python3, ecmwf-toolbox and vtb loaded.
"""

import argparse
import datetime as dt
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

MISSING = 3e38
G = 9.80665
# MARS parameters needed for each VINO parameter
MARS_PARAMS = {"t": ["t"], "z": ["z"], "ff": ["u", "v"], "dd": ["u", "v"], "td": ["t", "q"]}
SFC_PARAMS = {"2t": ["2t"], "2d": ["2d"], "10u": ["10u"], "10v": ["10v"],
              "10ff": ["10u", "10v"], "10dd": ["10u", "10v"]}
UNITS = {"t": "K", "td": "K", "z": "m2/s2", "ff": "m/s", "dd": "deg",
         "2t": "K", "2d": "K", "10u": "m/s", "10v": "m/s", "10ff": "m/s", "10dd": "deg"}


def log(*args):
    print(f"[{dt.datetime.now():%Y-%m-%d %H:%M:%S}]", *args, flush=True)


def day_dir(root, d):
    return Path(root) / f"{d:%Y}" / f"{d:%m}" / f"{d:%d}"


def write_geo(path, rows, ncols, meta):
    """rows: (stnid, lat, lon, elevation, value_0 ... value_<ncols-1>), NaN = missing."""
    cols = ["stnid", "latitude", "longitude", "level", "date", "time", "elevation"]
    cols += [f"value_{i}" for i in range(ncols)]
    date, time = meta["date"], int(meta["time"]) // 100
    lines = ["#GEO", "#FORMAT NCOLS", "#COLUMNS", "\t".join(cols),
             "# Missing values represented by 3e+38", "#METADATA"]
    lines += [f"{k}={v}" for k, v in meta.items()]
    lines.append("#DATA")
    for stnid, lat, lon, elev, *vals in rows:
        elev, *vals = ["3e+38" if not np.isfinite(v) else f"{v:.6g}" for v in (elev, *vals)]
        lines.append("\t".join([stnid, f"{lat:.5f}", f"{lon:.5f}", "0", date, str(time), elev, *vals]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def read_stations(path):
    """(stnids, lats, lons) of a geopoint file written by `obs`."""
    ids, lats, lons = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#DATA"):
                break
        for line in fh:
            f = line.split()
            if len(f) >= 3:
                ids.append(f[0])
                lats.append(float(f[1]))
                lons.append(float(f[2]))
    return ids, np.array(lats), np.array(lons)


def unit_xyz(lat, lon):
    lat, lon = np.deg2rad(lat), np.deg2rad(lon)
    return np.c_[np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]


# --------------------------------------------------------------------------
# obs
# --------------------------------------------------------------------------

def cmd_obs(a):
    from vtb.media.stvl import retrieve_to_geopoints
    import metview as mv

    for d in a.dates:
        vt = dt.datetime.strptime(d, "%Y%m%d").replace(hour=a.hour)
        data = list(retrieve_to_geopoints(
            table="observation", reference_datetimes=[vt], forecast_lengths=[0],
            parameter=a.parameter, period=0, sources=a.sources,
            level_type="pl", level=a.level))
        if not data:
            log(f"obs {a.key} {vt:%Y%m%d%H}: no data found")
            continue
        out = day_dir(a.outdir, vt) / f"{a.key}_obs_{vt:%Y%m%d%H}.geo"
        out.parent.mkdir(parents=True, exist_ok=True)
        mv.write(str(out), data[0])
        log(f"obs {a.key} {vt:%Y%m%d%H}: {len(data[0])} stations")


# --------------------------------------------------------------------------
# fc
# --------------------------------------------------------------------------

def derive(param, fields, level):
    """VINO parameter from the MARS fields (dict shortName -> values)."""
    if param in ("t", "z", "2t", "2d", "10u", "10v"):
        return fields[param]
    if param in ("ff", "10ff"):
        u, v = (fields["10u"], fields["10v"]) if param == "10ff" else (fields["u"], fields["v"])
        return np.hypot(u, v)
    if param in ("dd", "10dd"):  # direction the wind blows from
        u, v = (fields["10u"], fields["10v"]) if param == "10dd" else (fields["u"], fields["v"])
        return np.rad2deg(np.arctan2(-u, -v)) % 360.0
    if param == "td":  # inverse Magnus formula over water (IFS constants)
        q, p = fields["q"], level * 100.0
        e = np.maximum(q * p / (0.622 + 0.378 * q), 1e-3)
        x = np.log(e / 611.21)
        return (17.502 * 273.16 - 32.19 * x) / (17.502 - x)
    raise ValueError(param)


def mars(a, req):
    text = "retrieve," + ",".join(f"{k}={v}" for k, v in req.items()) + "\n"
    return subprocess.run([a.mars], input=text, text=True).returncode


def model_orography(a, req, tmp):
    """(KD-tree of the grid points, model surface height in m) from
    --orog_file, or else z at step 0 of the forecast itself. (None, None)
    if there is none: the stations then get no model height."""
    import metview as mv
    from scipy.spatial import cKDTree
    path = a.orog_file
    if not path:
        path = os.path.join(tmp, "orog.grib")
        zreq = dict(req, param="z", step=0, target=f'"{path}"')
        if a.numbers:
            zreq["number"] = a.numbers.split("/")[0]
        if mars(a, zreq) or not os.path.exists(path) or not os.path.getsize(path):
            log("WARNING no orography (z at step 0) in MARS, the stations get no model height")
            return None, None
    f = mv.read(path)[0]
    tree = cKDTree(unit_xyz(f.latitudes(), f.longitudes()))
    return tree, np.asarray(f.values(), float) / G


def cmd_fc(a):
    import metview as mv

    members = [int(n) for n in a.numbers.split("/")] if a.numbers else [0]
    sfc = a.level is None
    params = (SFC_PARAMS if sfc else MARS_PARAMS)[a.parameter]
    orog = None  # (KD-tree, heights in m) of the model orography, surface parameters only
    failed = 0
    for d in a.dates:
        base = dt.datetime.strptime(d, "%Y%m%d").replace(hour=a.time)
        # stations: those of the observation file valid at each step
        steps = {}
        for s in a.steps:
            vt = base + dt.timedelta(hours=s)
            of = day_dir(a.obsdir, vt) / f"{a.key}_obs_{vt:%Y%m%d%H}.geo"
            if of.exists():
                steps[s] = of
        if not steps:
            log(f"fc {a.key} {base:%Y%m%d%H}: no observation files, skipped")
            continue
        with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
            target = os.path.join(tmp, "fc.grib")
            req = {"class": a.mclass, "stream": a.stream, "expver": a.expver, "type": a.mtype,
                   "levtype": "sfc" if sfc else "pl", "param": "/".join(params), "date": d,
                   "time": f"{a.time:02d}", "step": "/".join(map(str, sorted(steps))),
                   "target": f'"{target}"'}
            if not sfc:
                req["levelist"] = a.level
            if a.grid:
                req["grid"] = a.grid
            if a.numbers:
                req["number"] = a.numbers
            for k in ("anoffset", "origin", "database"):
                if getattr(a, k):
                    req[k] = getattr(a, k)
            rc = mars(a, req)
            if sfc and orog is None and rc == 0:
                orog = model_orography(a, req, tmp)
            try:
                if rc or not os.path.getsize(target):
                    raise RuntimeError(f"mars rc={rc}")
                fs = mv.read(target)
                keys = fs.grib_get(["shortName", "step", "number"])
            except Exception as exc:
                log(f"fc {a.key} {base:%Y%m%d%H}: MARS retrieval failed: {exc}")
                failed += 1
                continue
            index = {(k[0], int(k[1]), int(k[2]) if a.numbers else 0): i for i, k in enumerate(keys)}
            tree = None
            for s, of in sorted(steps.items()):
                need = [(p, s, n) for n in members for p in params]
                if any(k not in index for k in need):
                    log(f"fc {a.key} {base:%Y%m%d%H} +{s}: fields missing in MARS")
                    continue
                ids, lats, lons = read_stations(of)
                if not ids:
                    continue
                if a.interpolation == "bilinear":
                    gp = mv.read(str(of))
                    pick = lambda f: np.asarray(mv.interpolate(f, gp).values(), float)
                else:
                    if tree is None:
                        from scipy.spatial import cKDTree
                        f0 = fs[index[need[0]]]
                        tree = cKDTree(unit_xyz(f0.latitudes(), f0.longitudes()))
                    _, idx = tree.query(unit_xyz(lats, lons))
                    pick = lambda f: np.asarray(f.values(), float)[idx]
                vals = []
                for n in members:
                    fields = {p: pick(fs[index[(p, s, n)]]) for p in params}
                    vals.append(derive(a.parameter, fields, a.level))
                vals = np.array(vals).T
                vals[~np.isfinite(vals) | (np.abs(vals) > 1e30)] = np.nan
                elev = np.full(len(ids), np.nan)
                if orog and orog[0] is not None:
                    elev = orog[1][orog[0].query(unit_xyz(lats, lons))[1]]
                meta = {"class": a.mclass, "stream": a.stream, "expver": a.expver, "type": a.mtype,
                        "param": a.parameter, "level": 0 if sfc else a.level,
                        "leveltype": "sfc" if sfc else "pl",
                        "date": d, "time": f"{a.time:02d}00", "step": s, "grid": a.grid or "native",
                        "units": UNITS[a.parameter], "source": f"vino_upper.py ({a.interpolation})"}
                name = f"{a.key}_mod_{base:%Y%m%d%H}_{s}_{a.stream}_{a.expver}_{a.mtype}.geo"
                write_geo(day_dir(a.outdir, base) / name,
                          [(i, la, lo, z, *v) for i, la, lo, z, v in zip(ids, lats, lons, elev, vals)],
                          len(members), meta)
            log(f"fc {a.key} {base:%Y%m%d%H}: {len(steps)} steps")
    if failed:
        log(f"{failed} dates failed")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("obs", "fc"):
        sp = sub.add_parser(name)
        sp.add_argument("--parameter", required=True, choices=sorted(MARS_PARAMS | SFC_PARAMS))
        sp.add_argument("--level", type=int, help="pressure level in hPa (none: surface)")
        sp.add_argument("--key", required=True, help="file name key, e.g. t850")
        sp.add_argument("--dates", required=True, nargs="+", help="YYYYMMDD ...")
        sp.add_argument("--outdir", required=True)
    o = sub.choices["obs"]
    o.add_argument("--hour", required=True, type=int)
    o.add_argument("--sources", nargs="+", default=["temp"])
    f = sub.choices["fc"]
    f.add_argument("--time", required=True, type=int)
    f.add_argument("--steps", required=True, nargs="+", type=int)
    f.add_argument("--class", dest="mclass", required=True)
    f.add_argument("--stream", required=True)
    f.add_argument("--expver", required=True)
    f.add_argument("--type", dest="mtype", required=True)
    f.add_argument("--anoffset")
    f.add_argument("--origin")
    f.add_argument("--database", help="alternative MARS database, e.g. fdb")
    f.add_argument("--numbers", help="ensemble members, e.g. 1/2/3")
    f.add_argument("--grid", help="MARS grid, e.g. O1280 (default: native)")
    f.add_argument("--orog_file", help="GRIB file of z for the model height at the stations "
                   "(surface; default: z at step 0 of the forecast)")
    f.add_argument("--interpolation", default="nearest_point", choices=["nearest_point", "bilinear"])
    f.add_argument("--obsdir", required=True, help="store directory of the observations")
    f.add_argument("--mars", default="mars", help="MARS client")
    a = ap.parse_args()
    {"obs": cmd_obs, "fc": cmd_fc}[a.cmd](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
