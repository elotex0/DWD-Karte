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
from zoneinfo import ZoneInfo
import numpy as np
from matplotlib.colors import ListedColormap

# Eingabe-/Ausgabe-Verzeichnisse
data_dir = sys.argv[1]
output_dir = sys.argv[2]
var_type = sys.argv[3]  # 't2m' oder 'ww'
os.makedirs(output_dir, exist_ok=True)

# Geo-Daten laden
bundeslaender = gpd.read_file("scripts/bundeslaender.geojson")
cities = pd.DataFrame({
    'name': ['Berlin', 'Hamburg', 'München', 'Köln', 'Frankfurt', 'Dresden', 'Stuttgart', 'Düsseldorf'],
    'lat': [52.52, 53.55, 48.14, 50.94, 50.11, 51.05, 48.78, 51.23],
    'lon': [13.40, 9.99, 11.57, 6.96, 8.68, 13.73, 9.18, 6.78]
})

ignore_codes = {99}

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
    95: "#FF77FF", 96: "#C71585"
}

ww_categories = {
    "Bewölkung": [0, 1 , 2, 3],
    "Nebel": [45, 48],
    "Schneeregen": [56, 57],
    "Regen": [51, 61, 80, 53, 63, 81, 55, 65, 82],
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

# Kartenextent für Deutschland
germany_bounds = bundeslaender.total_bounds
extent = [germany_bounds[0]-1, germany_bounds[2]+1, germany_bounds[1]-1, germany_bounds[3]+1]

# Funktion für WW-Legende unterhalb der Karte
def add_ww_legend_bottom(fig, present_codes, ww_categories, ww_colors_base):
    codes_for_legend = []
    labels_for_legend = []
    colors_for_legend = []

    for label, codes in ww_categories.items():
        c = next((c for c in codes if c in present_codes), None)
        if c is not None:
            codes_for_legend.append(c)
            labels_for_legend.append(label)
            colors_for_legend.append(ww_colors_base[c])

    n = len(colors_for_legend)
    if n == 0:
        return

    # Achse unterhalb der Karte
    legend_height = 0.08
    legend_ax = fig.add_axes([0.1, 0.02, 0.8, legend_height])  # Unterhalb der Karte
    legend_ax.set_xlim(0, n)
    legend_ax.set_ylim(0, 2)
    legend_ax.axis("off")

    # Farbkästchen oben
    for i, color in enumerate(colors_for_legend):
        legend_ax.add_patch(mpatches.Rectangle((i, 1), 1, 1, facecolor=color, edgecolor='black'))

    # Beschriftungen darunter
    for i, label in enumerate(labels_for_legend):
        legend_ax.text(i + 0.5, 0.5, label, ha='center', va='center', fontsize=8)

# Schleife über Dateien
for filename in sorted(os.listdir(data_dir)):
    if not filename.endswith(".grib2"):
        continue
    path = os.path.join(data_dir, filename)
    ds = cfgrib.open_dataset(path)

    # Daten aus GRIB
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

    # Figur & Achse
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent(extent)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.15)  # Platz für WW-Legende

    if var_type == "t2m":
        im = ax.pcolormesh(lon, lat, data, cmap=t2m_cmap, norm=t2m_norm, shading="auto")
        cbar = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.04, pad=0.04)
        cbar.set_ticks(list(range(-30, 45, 5)))
        cbar.set_label("Temperatur 2m [°C]", color="black")
        cbar.ax.tick_params(colors="black", labelsize=8)
        cbar.outline.set_edgecolor("black")
        cbar.ax.set_facecolor("white")
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

        # WW-Legende unterhalb der Karte
        add_ww_legend_bottom(fig, present_codes, ww_categories, ww_colors_base)

    # Bundesländer & Städte
    bundeslaender.boundary.plot(ax=ax, edgecolor="black", linewidth=1)
    for _, city in cities.iterrows():
        ax.plot(city["lon"], city["lat"], "ko", markersize=4)
        ax.text(city["lon"] + 0.1, city["lat"] + 0.1, city["name"], fontsize=8)
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE)

    # Titel
    ax.set_title(
        f"ICON-D2 {var_type.upper()} Vorhersage\n"
        f"Lauf: {pd.to_datetime(run_time_utc).hour if run_time_utc is not None else '??'}Z, "
        f"gültig: {valid_time_local:%d.%m.%Y %H:%M} Uhr (MEZ/MESZ)",
        pad=20
    )

    # Speichern
    outname = f"{var_type}_{pd.to_datetime(valid_time_utc):%Y%m%d_%H%M}.png"
    plt.savefig(os.path.join(output_dir, outname), dpi=150, bbox_inches='tight')
    plt.close()
