"""Frequency of 2d > 2t at SYNOP stations: jc0v member 1 (first pf) vs observations.
Uses the collocated store files of run AIFS_v3_MJJA26 (MJJA 2026, 00/12 UTC, steps 24..216)."""
import sys, glob, datetime as dt
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeature

sys.path.insert(0, "/perm/pavn/spatial_synop_maps")
from synop_maps import read_geo, station_key, expand_dates

ROOT = Path("/perm/pavn/spatial_synop_maps/store")
FC = ROOT / "fc/local/rd/enfo/jc0v/pf/numbers1xtox8/args823d719a/nearest_point_synop"
OBS = ROOT / "obs/synop"
OUT = Path(sys.argv[1])
dates = expand_dates("20260513/to/20260831/by/1")
steps = range(24, 241, 48)
MIN_N = 20

def load(path, col):
    df, _ = read_geo(str(path))
    df = df.assign(key=station_key(df)).drop_duplicates("key")
    return df.set_index("key")[["latitude", "longitude", col]]

fc_rows, obs_cache, obs_rows = [], {}, []
for d in dates:
    for t in (0, 12):
        for s in steps:
            day = FC / f"{d:%Y/%m/%d}"
            f2t = glob.glob(str(day / f"2t_mod_{d:%Y%m%d}{t:02d}_{s}_*.geo"))
            f2d = glob.glob(str(day / f"2d_mod_{d:%Y%m%d}{t:02d}_{s}_*.geo"))
            if not (f2t and f2d):
                continue
            a = load(f2t[0], "value_0").rename(columns={"value_0": "t"})
            b = load(f2d[0], "value_0")[["value_0"]].rename(columns={"value_0": "td"})
            m = a.join(b, how="inner").dropna()
            m["step"] = s
            fc_rows.append(m.reset_index())
            vt = dt.datetime(d.year, d.month, d.day, t) + dt.timedelta(hours=s)
            if vt not in obs_cache:
                obs_cache[vt] = True
                o2t = OBS / f"{vt:%Y/%m/%d}" / f"2t_obs_{vt:%Y%m%d%H}.geo"
                o2d = OBS / f"{vt:%Y/%m/%d}" / f"2d_obs_{vt:%Y%m%d%H}.geo"
                if o2t.exists() and o2d.exists():
                    a = load(o2t, "value_0").rename(columns={"value_0": "t"})
                    b = load(o2d, "value_0")[["value_0"]].rename(columns={"value_0": "td"})
                    o = a.join(b, how="inner").dropna()
                    o = o[o["t"].between(173.15, 333.15) & o["td"].between(173.15, 333.15)]
                    obs_rows.append(o.reset_index())

fc = pd.concat(fc_rows, ignore_index=True)
obs = pd.concat(obs_rows, ignore_index=True)
fc["sup"] = fc["td"] > fc["t"]
obs["sup"] = obs["td"] > obs["t"]
print(f"fc pairs {len(fc)}, {fc['sup'].mean()*100:.2f}% with 2d>2t; "
      f"obs pairs {len(obs)}, {obs['sup'].mean()*100:.2f}%")

def per_station(df):
    g = df.groupby("key").agg(lat=("latitude", "first"), lon=("longitude", "first"),
                              n=("sup", "size"), f=("sup", "mean"))
    g = g[g["n"] >= MIN_N]
    g["f"] *= 100
    return g

st_fc, st_obs = per_station(fc), per_station(obs)
bystep = fc.groupby("step")["sup"].mean() * 100

levels = [0, 0.5, 1, 2, 5, 10, 20, 50, 100]
cmap = plt.get_cmap("YlOrRd", len(levels) - 1)
norm = matplotlib.colors.BoundaryNorm(levels, cmap.N)

def draw(ax, st, title, extent=None):
    ax.add_feature(cfeature.LAND, facecolor="#eeeeee")
    ax.coastlines(lw=0.4, color="#666666")
    if extent:
        ax.set_extent(extent, ccrs.PlateCarree())
    else:
        ax.set_global()
    st = st.sort_values("f")
    zero = st[st["f"] == 0]
    pos = st[st["f"] > 0]
    ms = 1.5 if extent is None else 6
    ax.scatter(zero["lon"], zero["lat"], s=ms * 0.5, c="#bbbbbb", lw=0, transform=ccrs.PlateCarree())
    sc = ax.scatter(pos["lon"], pos["lat"], c=pos["f"], s=ms, cmap=cmap, norm=norm, lw=0,
                    transform=ccrs.PlateCarree())
    ax.set_title(title, fontsize=10)
    return sc

for region, extent in [("global", None), ("europe", [-25, 45, 33, 72])]:
    fig = plt.figure(figsize=(13, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[3.3, 1.6], hspace=0.45, wspace=0.08)
    ax1 = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree())
    ax2 = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())
    sc = draw(ax1, st_fc, f"jc0v member 1: {fc['sup'].mean()*100:.1f}% of all station forecasts\n"
              f"({len(st_fc)} stations, all steps 24–216 h)", extent)
    draw(ax2, st_obs, f"SYNOP observations: {obs['sup'].mean()*100:.1f}% of all reports\n"
         f"({len(st_obs)} stations, same valid times)", extent)
    cb = fig.colorbar(sc, ax=[ax1, ax2], orientation="horizontal", fraction=0.04, pad=0.04,
                      ticks=levels, extend="neither")
    cb.ax.set_xticklabels([f"{l:g}" for l in levels])
    cb.set_label("Frequency of 2d > 2t (%) per station   (grey: never)")
    ax3 = fig.add_subplot(gs[1, :])
    ax3.bar([str(s) for s in bystep.index], bystep.values, color="#d95f0e", width=0.6)
    ax3.axhline(obs["sup"].mean() * 100, color="#333333", ls="--", lw=1,
                label=f"observed ({obs['sup'].mean()*100:.2f}%)")
    for x, v in zip(range(len(bystep)), bystep.values):
        ax3.text(x, v, f"{v:.2f}%", ha="center", va="bottom", fontsize=9)
    ax3.set_xlabel("Forecast step (h)")
    ax3.set_ylabel("% of station forecasts\nwith 2d > 2t")
    ax3.spines[["top", "right"]].set_visible(False)
    ax3.set_ylim(0, bystep.max() * 1.2)
    ax3.legend(frameon=False, loc="lower right", bbox_to_anchor=(1, 1.0))
    fig.suptitle("Frequency of dew point above temperature, SYNOP stations, MJJA 2026 "
                 "(13 May–31 Aug, 00/12 UTC)", fontsize=12, y=0.97)
    path = OUT / f"freq_2d_gt_2t_jc0v_m1_{region}.png"
    fig.savefig(path, dpi=130, bbox_inches="tight")
    print(path)
