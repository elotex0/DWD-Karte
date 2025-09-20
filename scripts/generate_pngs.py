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
from matplotlib.colors import ListedColormap

# Eingabe-/Ausgabe-Verzeichnisse
data_dir = sys.argv[1]
output_dir = sys.argv[2]
var_type = sys.argv[3]  # 't2m' oder 'ww'
os.makedirs(output_dir, exist_ok=True)

# Geo-Daten
bundeslaender = gpd.read_file("scripts/bundeslaender.geojson")
cities = pd.DataFrame({
    'name': ['Berlin', 'Hamburg', 'München', 'Köln', 'Frankfurt', 'Dresden', 'Stuttgart', 'Düsseldorf'],
    'lat': [52.52, 53.55, 48.14, 50.94, 50.11, 51.05, 48.78, 51.23],
    'lon': [13.40, 9.99, 11.57, 6.96, 8.68, 13.73, 9.18, 6.78]
})

ignore_codes = {4}

# WW-Farben
ww_colors_base = {
    0: "#FFFFFF", 1: "#D3D3D3", 2: "#A9A9A9", 3: "#696969",
    45: "#FFFF00", 48: "#FFD700",
    56: "#FFA500", 57: "#FF8C00",
    51: "#90EE90", 53: "#32CD32", 55: "#006400",
    61: "#90EE90", 63: "#32CD32", 65: "#006400",
    80: "#90EE90", 81: "#32CD32", 82: "#006400",
    66: "#FF6347", 67: "#8B0000",
    71: "#ADD8E6", 73: "#6495ED", 75: "#00008B",
    95: "#FF77FF", 96: "#C71585", 99: "#C71585"
}

ww_categories = {
    "Bewölkung": [0, 1 , 2, 3],
    "Nebel": [45, 48],
    "Schneeregen": [56, 57],
    "Regen": [51, 53, 55],
    "gefr. Regen": [66, 67],
    "Schnee": [71, 73, 75],
    "Gewitter": [95,96],
}

# Temperaturfarbskala
t2m_bounds = list(range(-30, 45, 5))
t2m_colors = [
    "#001070", "#0020c2", "#0040ff", "#0080ff", "#00c0ff", "#00ffff",
    "#80ff80", "#c0ff00", "#ffff00", "#ffcc00", "#ff8000", "#ff4000",
    "#ff0000", "#990000"
]
t2m_cmap = mcolors.ListedColormap(t2m_colors)
t2m_norm = mcolors.BoundaryNorm(t2m_bounds, t2m_cmap.N)

# Ziel-Pixelmaße und Ziel-Seitenverhältnis für die Kartenfläche oben
FIG_W_PX, FIG_H_PX = 880, 830
BOTTOM_AREA_PX = 179
TOP_AREA_PX = FIG_H_PX - BOTTOM_AREA_PX  # 651 px
TARGET_ASPECT = FIG_W_PX / TOP_AREA_PX   # ~1.3525

# Kartenextent: Deutschland + kleiner Rand, an TARGET_ASPECT angepasst
_minx, _miny, _maxx, _maxy = bundeslaender.total_bounds
_w = _maxx - _minx
_h = _maxy - _miny
base_pad_x = _w * 0.05  # seitlich etwas Rand
base_pad_y = _h * 0.03  # oben/unten etwas weniger
xmin = _minx - base_pad_x
xmax = _maxx + base_pad_x
ymin = _miny - base_pad_y
ymax = _maxy + base_pad_y
w = xmax - xmin
h = ymax - ymin
if w / h < TARGET_ASPECT:
    # zu schmal → Breite erweitern
    needed_w = h * TARGET_ASPECT
    extra = (needed_w - w) / 2
    xmin -= extra
    xmax += extra
else:
    # zu flach → Höhe erweitern
    needed_h = w / TARGET_ASPECT
    extra = (needed_h - h) / 2
    ymin -= extra
    ymax += extra
extent = [xmin, xmax, ymin, ymax]

# WW-Legende unten
def add_ww_legend_bottom(fig, ww_categories, ww_colors_base):
    legend_height = 0.05
    legend_ax = fig.add_axes([0.1, 0.01, 0.8, legend_height])
    legend_ax.axis("off")

    categories_present = [(label, codes) for label, codes in ww_categories.items()]
    n_categories = len(categories_present)
    if n_categories == 0:
        return

    total_width = 1.0
    block_width = total_width / n_categories
    gap = 0.02 * block_width

    for i, (label, codes_in_cat) in enumerate(categories_present):
        x0 = i * block_width
        x1 = (i + 1) * block_width
        block_inner_width = x1 - x0 - gap

        n_colors = len(codes_in_cat)
        color_width = block_inner_width / n_colors

        for j, c in enumerate(codes_in_cat):
            color = ww_colors_base.get(c, "#FFFFFF")
            legend_ax.add_patch(
                mpatches.Rectangle((x0 + j * color_width, 0.5),
                                   color_width, 0.5,
                                   facecolor=color, edgecolor='black')
            )

        legend_ax.text((x0 + x1)/2, 0.25, label, ha='center', va='center', fontsize=8)

# Schleife über Dateien
for filename in sorted(os.listdir(data_dir)):
    if not filename.endswith(".grib2"):
        continue
    path = os.path.join(data_dir, filename)
    ds = cfgrib.open_dataset(path)

    # Daten
    if var_type == "t2m":
        if "t2m" not in ds:
            print(f"Keine t2m-Variable in {filename}")
            continue
        data = ds["t2m"].values - 273.15
    else:
        varname = next((vn for vn in ds.data_vars if vn.lower() in ["ww", "weather"]), None)
        if varname is None:
            print(f"Keine WW-Variable in {filename}")
            continue
        data = ds[varname].values

    if data.ndim == 3:
        data = data[0, :, :]
    lon = ds["longitude"].values
    lat = ds["latitude"].values
    run_time_utc = pd.to_datetime(ds["time"].values) if "time" in ds else None
    valid_time_utc = pd.to_datetime(ds.valid_time.values) if "valid_time" in ds else None
    if isinstance(valid_time_utc, (np.ndarray, list)):
        valid_time_utc = valid_time_utc[0]
    valid_time_local = pd.to_datetime(valid_time_utc).tz_localize("UTC").astimezone(ZoneInfo("Europe/Berlin"))

    # Figur und Achsen-Layout mit exakten Pixeln: 880x830 gesamt, Kartenbereich oben 880x651, unten 179 px
    FIG_W_PX, FIG_H_PX = 880, 830
    TOP_H_PX = FIG_H_PX - BOTTOM_AREA_PX  # 651 px

    # Figure mit fixen Pixeln (dpi so wählen, dass figsize*dpi=Pixel)
    fig = plt.figure(figsize=(FIG_W_PX/100, FIG_H_PX/100), dpi=100)

    # Kartenachse füllt die gesamte Breite oben (keine weißen Ränder links/rechts)
    ax_left = 0.0
    ax_bottom = BOTTOM_AREA_PX / FIG_H_PX
    ax_width = 1.0
    ax_height = TOP_H_PX / FIG_H_PX

    ax = fig.add_axes([ax_left, ax_bottom, ax_width, ax_height], projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.set_axis_off()  # Achsen komplett ausblenden für randlose Karte

    # Plot
    if var_type == "t2m":
        im = ax.pcolormesh(lon, lat, data, cmap=t2m_cmap, norm=t2m_norm, shading="auto")
    else:
        valid_mask = np.isfinite(data)
        present_codes = np.unique(data[valid_mask]).astype(int)
        present_codes = [c for c in present_codes if c in ww_colors_base and c not in ignore_codes]
        present_codes.sort()
        colors = [ww_colors_base[c] for c in present_codes]
        cmap = ListedColormap(colors)
        code2idx = {c: i for i, c in enumerate(present_codes)}
        idx_data = np.full_like(data, fill_value=np.nan, dtype=float)
        for c, i in code2idx.items():
            idx_data[data == c] = i
        im = ax.pcolormesh(lon, lat, idx_data, cmap=cmap, vmin=-0.5, vmax=len(colors)-0.5, shading="auto")

    # Bundesländer & Städte
    bundeslaender.boundary.plot(ax=ax, edgecolor="black", linewidth=1)
    for _, city in cities.iterrows():
        ax.plot(city["lon"], city["lat"], "o", markersize=6, markerfacecolor="black",
                markeredgecolor="white", markeredgewidth=1.5, zorder=5)
        txt = ax.text(city["lon"] + 0.1, city["lat"] + 0.1, city["name"], fontsize=8, color="black", zorder=6)
        txt.set_path_effects([path_effects.withStroke(linewidth=1.5, foreground="white")])
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE)

    # Legende ganz unten (immer unten), darüber der Titel
    # Wir teilen die 179 px: Legende 69 px unten, Titel 110 px darüber
    legend_h_px = 69
    footer_h_px = BOTTOM_AREA_PX - legend_h_px  # 110 px

    # Legende (ganz unten)
    if var_type == "t2m":
        cbar_ax = fig.add_axes([
            0.03,
            0.0,
            0.94,
            legend_h_px / FIG_H_PX
        ])
        cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
        cbar.set_ticks(list(range(-30, 45, 5)))
        cbar.set_label("Temperatur 2m [°C]", color="black")
        cbar.ax.tick_params(colors="black", labelsize=8)
        cbar.outline.set_edgecolor("black")
        cbar.ax.set_facecolor("white")
    else:
        legend_ax = fig.add_axes([
            0.03,
            0.0,
            0.94,
            legend_h_px / FIG_H_PX
        ])
        legend_ax.axis("off")
        categories_present = [(label, codes) for label, codes in ww_categories.items()]
        n_categories = len(categories_present)
        if n_categories > 0:
            total_width = 1.0
            block_width = total_width / n_categories
            gap = 0.02 * block_width
            for i, (label, codes_in_cat) in enumerate(categories_present):
                x0 = i * block_width
                x1 = (i + 1) * block_width
                block_inner_width = x1 - x0 - gap
                n_colors = len(codes_in_cat)
                color_width = block_inner_width / n_colors
                for j, c in enumerate(codes_in_cat):
                    color = ww_colors_base.get(c, "#FFFFFF")
                    legend_ax.add_patch(
                        mpatches.Rectangle((x0 + j * color_width, 0.5),
                                           color_width, 0.5,
                                           facecolor=color, edgecolor='black')
                    )
                legend_ax.text((x0 + x1)/2, 0.25, label, ha='center', va='center', fontsize=8)

    # Footer (Titel) direkt über der Legende
    footer_ax = fig.add_axes([
        0.0,
        legend_h_px / FIG_H_PX,
        1.0,
        footer_h_px / FIG_H_PX
    ])
    footer_ax.axis("off")
    left_text = "Signifikantes Wetter" if var_type=="ww" else "Temperatur 2m"
    left_text += f"\nICON-D2 ({pd.to_datetime(run_time_utc).hour if run_time_utc else '??'}Z), Deutscher Wetterdienst"
    footer_ax.text(0.01, 0.85, left_text, fontsize=10, fontweight="bold", va="top", ha="left")
    footer_ax.text(0.99, 0.85, f"{valid_time_local:%d.%m.%Y %H:%M} Uhr", fontsize=10, va="top", ha="right")

    # Speichern mit exakter Pixelgröße 880x830
    outname = f"{var_type}_{pd.to_datetime(valid_time_utc):%Y%m%d_%H%M}.png"
    # Da Figure auf dpi=100 und figsize=(8.8, 8.3) gesetzt ist, ergibt das exakt 880x830
    plt.savefig(os.path.join(output_dir, outname), dpi=100, bbox_inches=None, pad_inches=0)
    plt.close()