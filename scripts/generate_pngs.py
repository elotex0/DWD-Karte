import sys
import cfgrib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import pandas as pd
import os
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.patheffects as path_effects
from zoneinfo import ZoneInfo
import numpy as np
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

# ------------------------------
# Eingabe-/Ausgabe
# ------------------------------
data_dir = sys.argv[1]        # z.B. "output"
output_dir = sys.argv[2]      # z.B. "output/maps"
var_type = sys.argv[3]        # 't2m', 'ww', 'tp', 'cape_ml', 'dbz_cmax'
os.makedirs(output_dir, exist_ok=True)

# ------------------------------
# Geo-Daten
# ------------------------------
bundeslaender = gpd.read_file("scripts/bundeslaender.geojson")
cities = pd.DataFrame({
    'name': ['Berlin', 'Hamburg', 'München', 'Köln', 'Frankfurt', 'Dresden', 'Stuttgart', 'Düsseldorf',
             'Nürnberg', 'Erfurt', 'Leipzig', 'Bremen', 'Saarbrücken', 'Hannover'],
    'lat': [52.52, 53.55, 48.14, 50.94, 50.11, 51.05, 48.78, 51.23,
            49.45, 50.98, 51.34, 53.08, 49.24, 52.37],
    'lon': [13.40, 9.99, 11.57, 6.96, 8.68, 13.73, 9.18, 6.78,
            11.08, 11.03, 12.37, 8.80, 6.99, 9.73]
})

ignore_codes = {4}

# ------------------------------
# WW-Farben
# ------------------------------
ww_colors_base = {
    0: "#FFFFFF", 1: "#D3D3D3", 2: "#A9A9A9", 3: "#696969",
    45: "#FFFF00", 48: "#FFD700",
    56: "#FFA500", 57: "#C06A00",
    51: "#A3FFA3", 53: "#33FF33", 55: "#006600",
    61: "#33FF33", 63: "#009900", 65: "#006600",
    80: "#33FF33", 81: "#009900", 82: "#006600",
    66: "#FF6347", 67: "#8B0000",
    71: "#ADD8E6", 73: "#6495ED", 75: "#00008B",
    95: "#FF77FF", 96: "#C71585", 99: "#C71585"
}
ww_categories = {
    "Bewölkung": [0, 1 , 2, 3],
    "Nebel": [45],
    "Schneeregen": [56, 57],
    "Regen": [51, 61, 63, 65],
    "gefr. Regen": [66, 67],
    "Schnee": [71, 73, 75],
    "Gewitter": [95,96],
}

# ------------------------------
# Temperatur-Farben
# ------------------------------
t2m_bounds = list(range(-20, 41, 2))
t2m_colors = [
    "#000080", "#0010A8", "#0020D0", "#0030F8", "#0050FF",
    "#007FFF", "#00AFFF", "#00DFFF", "#00FFBF", "#00FF9F",
    "#00FF7F", "#30FF60", "#60FF40", "#90FF20", "#BFFF00",
    "#E0FF00", "#FFFF00", "#FFE000", "#FFC000", "#FFA000",
    "#FF8000", "#FF6000", "#FF4000", "#FF2000", "#FF0000",
    "#E00000", "#C00000", "#A00000", "#800000", "#600000"
]
t2m_cmap = ListedColormap(t2m_colors)
t2m_norm = mcolors.BoundaryNorm(t2m_bounds, t2m_cmap.N)

# ------------------------------
# Niederschlags-Farben (tp)
# ------------------------------
prec_bounds = [0.1, 0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
               12, 14, 16, 20, 24, 30, 40, 50, 60, 80, 100, 125]
prec_colors = ListedColormap([
    "#A8EFFF", "#77CBFF", "#4FAAFF", "#1E90FF", "#0077CC",
    "#32CD32", "#228B22", "#006400",
    "#FFFF99", "#FFD700", "#FFA500",
    "#FF8C00", "#FF4500", "#FF0000", "#B22222",
    "#800080", "#A0007F", "#C0009F", "#E000BF",
    "#CC6699", "#996666", "#665555", "#555555", "#323232"
])
prec_norm = mcolors.BoundaryNorm(prec_bounds, prec_colors.N)

# ------------------------------
# CAPE-Farben
# ------------------------------
cape_bounds = [0, 20, 40, 60, 80, 100, 200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000]
cape_colors = ListedColormap([
    "#676767", "#006400", "#008000", "#00CC00", "#66FF00", "#FFFF00", 
    "#FFCC00", "#FF9900", "#FF6600", "#FF3300", "#FF0000", "#FF0095", 
    "#FC439F", "#FF88D3", "#FF99FF"
])
cape_norm = mcolors.BoundaryNorm(cape_bounds, cape_colors.N)

# ------------------------------
# DBZ-CMAX Farben
# ------------------------------
dbz_bounds = np.concatenate(([2, 6], np.arange(6, 74, 2)))
dbz_bounds = np.unique(dbz_bounds)
dbz_colors = [
    "#676767", "#0050A0", "#009966", "#00AA00",
    "#AAAA00", "#AA5500", "#AA0000", "#550055"
]
dbz_cmap = LinearSegmentedColormap.from_list("dbz_cmap", dbz_colors, N=len(dbz_bounds)-1)
dbz_norm = BoundaryNorm(dbz_bounds, dbz_cmap.N)

# ------------------------------
# Kartenparameter
# ------------------------------
FIG_W_PX, FIG_H_PX = 880, 830
BOTTOM_AREA_PX = 179
TOP_AREA_PX = FIG_H_PX - BOTTOM_AREA_PX
TARGET_ASPECT = FIG_W_PX / TOP_AREA_PX

_minx, _miny, _maxx, _maxy = bundeslaender.total_bounds
_w = _maxx - _minx
_h = _maxy - _miny
ymin, ymax = _miny, _maxy
h = ymax - ymin
left_pad_factor = 0.56
right_pad_factor = 0.34
xmin = _minx - _w * left_pad_factor
xmax = _maxx + _w * right_pad_factor
needed_w = h * TARGET_ASPECT
current_w = xmax - xmin
if current_w < needed_w:
    extra = (needed_w - current_w) / 2
    xmin -= extra
    xmax += extra
extent = [xmin, xmax, ymin, ymax]

# ------------------------------
# WW-Legende Funktion
# ------------------------------
def add_ww_legend_bottom(fig, ww_categories, ww_colors_base):
    legend_height = 0.12
    legend_ax = fig.add_axes([0.05, 0.01, 0.9, legend_height])
    legend_ax.axis("off")
    for i, (label, codes) in enumerate(ww_categories.items()):
        n_colors = len(codes)
        block_width = 1.0 / len(ww_categories)
        gap = 0.05 * block_width
        x0 = i * block_width
        x1 = (i + 1) * block_width
        inner_width = x1 - x0 - gap
        color_width = inner_width / n_colors
        for j, c in enumerate(codes):
            color = ww_colors_base.get(c, "#FFFFFF")
            legend_ax.add_patch(mpatches.Rectangle((x0 + j * color_width, 0.3),
                                                  color_width, 0.6,
                                                  facecolor=color, edgecolor='black'))
        legend_ax.text((x0 + x1)/2, 0.05, label, ha='center', va='bottom', fontsize=10)

# ------------------------------
# Dateien durchgehen
# ------------------------------
for filename in sorted(os.listdir(data_dir)):
    if not filename.endswith(".grib2"):
        continue
    path = os.path.join(data_dir, filename)
    ds = cfgrib.open_dataset(path)

    # Daten je Typ
    if var_type == "t2m":
        if "t2m" not in ds:
            print(f"Keine t2m in {filename}")
            continue
        data = ds["t2m"].values - 273.15
    elif var_type == "ww":
        varname = next((vn for vn in ds.data_vars if vn.lower() in ["ww","weather"]), None)
        if varname is None:
            print(f"Keine WW in {filename}")
            continue
        data = ds[varname].values
    elif var_type == "tp":
        tp_var = None
        for vn in ["tp", "tot_prec"]:
            if vn in ds:
                tp_var = vn
                break
        if tp_var is None:
            print(f"Keine Niederschlagsvariable in {filename}")
            continue
        tp_all = ds[tp_var].values
        if tp_all.shape[0] > 1:
            data = tp_all[3] - tp_all[0]
        else:
            data = tp_all[0]
        data[data < 0.01] = np.nan
    elif var_type == "cape_ml":
        if "CAPE_ML" not in ds:
            print(f"Keine CAPE_ML-Variable in {filename}")
            continue
        data = ds["CAPE_ML"].values[0, :, :]
        data[data < 0] = np.nan
    elif var_type == "dbz_cmax":
        if "DBZ_CMAX" not in ds:
            print(f"Keine DBZ_CMAX in {filename}")
            continue
        dbz_all_steps = ds["DBZ_CMAX"].values
    
        # Zeitschritte auswählen und Differenz berechnen
        dbz_step0 = dbz_all_steps[0, :, :]
        dbz_step3 = dbz_all_steps[3, :, :]
        dbz_diff = dbz_step3 - dbz_step0  # Differenz zwischen Step 3 und Step 0
        dbz_diff[dbz_diff < -100] = np.nan  # Werte < -100 dBZ als fehlend markieren
    
        data = dbz_diff
    else:
        print(f"Unbekannter var_type {var_type}")
        continue

    if data.ndim == 3:
        data = data[0]

    lon = ds["longitude"].values
    lat = ds["latitude"].values
    run_time_utc = pd.to_datetime(ds["time"].values) if "time" in ds else None

    if "valid_time" in ds:
        valid_time_raw = ds["valid_time"].values
        valid_time_utc = pd.to_datetime(valid_time_raw[0]) if np.ndim(valid_time_raw) > 0 else pd.to_datetime(valid_time_raw)
    else:
        step = pd.to_timedelta(ds["step"].values[0])
        valid_time_utc = run_time_utc + step
    valid_time_local = valid_time_utc.tz_localize("UTC").astimezone(ZoneInfo("Europe/Berlin"))

    # --------------------------
    # Figure
    # --------------------------
    scale = 0.9
    fig = plt.figure(figsize=(FIG_W_PX/100*scale, FIG_H_PX/100*scale), dpi=100)
    shift_up = 0.02
    ax = fig.add_axes([0.0, BOTTOM_AREA_PX / FIG_H_PX + shift_up, 1.0, TOP_AREA_PX / FIG_H_PX],
                      projection=ccrs.PlateCarree())
    ax.set_extent(extent)
    ax.set_axis_off()
    ax.set_aspect('auto')

    # Plot
    if var_type == "t2m":
        im = ax.pcolormesh(lon, lat, data, cmap=t2m_cmap, norm=t2m_norm, shading="auto")
    elif var_type == "ww":
        valid_mask = np.isfinite(data)
        codes = np.unique(data[valid_mask]).astype(int)
        codes = [c for c in codes if c in ww_colors_base and c not in ignore_codes]
        codes.sort()
        cmap = ListedColormap([ww_colors_base[c] for c in codes])
        code2idx = {c: i for i, c in enumerate(codes)}
        idx_data = np.full_like(data, fill_value=np.nan, dtype=float)
        for c,i in code2idx.items():
            idx_data[data==c] = i
        im = ax.pcolormesh(lon, lat, idx_data, cmap=cmap, vmin=-0.5, vmax=len(codes)-0.5, shading="auto")
    elif var_type == "tp":
        im = ax.pcolormesh(lon, lat, data, cmap=prec_colors, norm=prec_norm, shading="auto")
    elif var_type == "cape_ml":
        im = ax.pcolormesh(lon, lat, data, cmap=cape_colors, norm=cape_norm, shading="auto")
    elif var_type == "dbz_cmax":
        im = ax.pcolormesh(lon, lat, data, cmap=dbz_cmap, norm=dbz_norm, shading="auto")

    # Bundesländer & Städte
    bundeslaender.boundary.plot(ax=ax, edgecolor="black", linewidth=1)
    for _, city in cities.iterrows():
        ax.plot(city["lon"], city["lat"], "o", markersize=6, markerfacecolor="black",
                markeredgecolor="white", markeredgewidth=1.5, zorder=5)
        txt = ax.text(city["lon"]+0.1, city["lat"]+0.1, city["name"], fontsize=9,
                      color="black", weight="bold", zorder=6)
        txt.set_path_effects([path_effects.withStroke(linewidth=1.5, foreground="white")])
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE)
    ax.add_patch(mpatches.Rectangle((0,0),1,1, transform=ax.transAxes, fill=False, color="black", linewidth=2))

    # Legende
    legend_h_px = 50
    legend_bottom_px = 45
    if var_type in ["t2m", "tp", "cape_ml", "dbz_cmax"]:
        bounds = t2m_bounds if var_type=="t2m" else prec_bounds if var_type=="tp" else cape_bounds if var_type=="cape_ml" else dbz_bounds
        cbar_ax = fig.add_axes([0.03, legend_bottom_px / FIG_H_PX, 0.94, legend_h_px / FIG_H_PX])
        cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal", ticks=bounds)
        cbar.ax.tick_params(colors="black", labelsize=7)
        cbar.outline.set_edgecolor("black")
        cbar.ax.set_facecolor("white")
    else:
        add_ww_legend_bottom(fig, ww_categories, ww_colors_base)

    # Footer
    footer_ax = fig.add_axes([0.0, (legend_bottom_px + legend_h_px)/FIG_H_PX, 1.0,
                              (BOTTOM_AREA_PX - legend_h_px - legend_bottom_px)/FIG_H_PX])
    footer_ax.axis("off")
    footer_texts = {
        "ww": "Signifikantes Wetter",
        "t2m": "Temperatur in 2m (in °C)",
        "tp": "Niederschlag, 1Std (in mm)",
        "cape_ml": "CAPE-Index (in J/kg)",
        "dbz_cmax": "Sim. max. Radarreflektivität (in dBZ)"
    }
    left_text = footer_texts.get(var_type, var_type)
    left_text += f"\nICON-D2 ({pd.to_datetime(run_time_utc).hour if run_time_utc else '??'}z), Deutscher Wetterdienst"
    footer_ax.text(0.01, 0.85, left_text, fontsize=12, fontweight="bold", va="top", ha="left")
    footer_ax.text(0.734, 0.92, "Prognose für:", fontsize=12, va="top", ha="left", fontweight="bold")
    footer_ax.text(0.99, 0.68, f"{valid_time_local:%d.%m.%Y %H:%M} Uhr", fontsize=12, va="top", ha="right", fontweight="bold")

    # Speichern
    outname = f"{var_type}_{valid_time_local:%Y%m%d_%H%M}.png"
    plt.savefig(os.path.join(output_dir, outname), dpi=100, bbox_inches=None, pad_inches=0)
    plt.close()
