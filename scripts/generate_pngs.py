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

bundeslaender = gpd.read_file("scripts/bundeslaender.geojson")

cities = pd.DataFrame({
    'name': ['Berlin','Hamburg','München','Köln','Frankfurt','Dresden','Stuttgart','Düsseldorf'],
    'lat': [52.52,53.55,48.14,50.94,50.11,51.05,48.78,51.23],
    'lon': [13.40,9.99,11.57,6.96,8.68,13.73,9.18,6.78]
})

# Codes, die ignoriert werden sollen
ignore_codes = {51, 53, 55}

# WW-Farben
ww_colors_base = {
    0:"#FFFFFF", 1:"#D3D3D3", 2:"#A9A9A9", 3:"#696969",
    45:"#FFFF00", 48:"#FFD700",
    56:"#FFA500", 57:"#FF8C00",
    51:"#90EE90", 53:"#32CD32", 55:"#006400",
    61:"#90EE90", 63:"#32CD32", 65:"#006400",
    80:"#90EE90", 81:"#32CD32",
    66:"#FF6347", 67:"#8B0000",
    71:"#ADD8E6", 73:"#6495ED", 75:"#00008B",
    95:"#FF77FF", 96:"#C71585"
}

# WW-Labels
ww_labels = {
    0:"klar",1:"leicht bewölkt",2:"teilweise bewölkt",3:"bedeckt",
    45:"Nebel",48:"Nebel mit Reif",
    56:"Schneeregen leicht",57:"Schneeregen stark",
    61:"Regen leicht",63:"Regen mäßig",81:"Regen stark",
    66:"gef. Regen leicht",67:"gef. Regen stark",
    71:"Schnee leicht",73:"Schnee mäßig",75:"Schnee stark",
    95:"Gewitter leicht/mäßig",96:"Gewitter stark"
}

# Temperaturfarbskala
t2m_bounds = list(range(-30, 45, 5))
t2m_colors = [
    "#001070","#0020c2","#0040ff","#0080ff","#00c0ff","#00ffff",
    "#80ff80","#c0ff00","#ffff00","#ffcc00","#ff8000","#ff4000",
    "#ff0000","#990000"
]
t2m_cmap = mcolors.ListedColormap(t2m_colors)
t2m_norm = mcolors.BoundaryNorm(t2m_bounds, t2m_cmap.N)

for filename in sorted(os.listdir(data_dir)):
    if not filename.endswith(".grib2"):
        continue

    path = os.path.join(data_dir, filename)
    ds = cfgrib.open_dataset(path)

    # Variable auswählen
    if var_type == "t2m":
        if "t2m" not in ds:
            print(f"Keine t2m-Variable in {filename} gefunden. Variablen: {list(ds.data_vars)}")
            continue
        data = ds["t2m"].values - 273.15  # Kelvin → °C
    else:  # WW
        varname = None
        if "WW" in ds:
            varname = "WW"
        elif "ww" in ds:
            varname = "ww"
        else:
            for vn in ds.data_vars:
                if vn.lower() == "ww" or "weather" in vn.lower():
                    varname = vn
                    break
        if varname is None:
            print(f"Keine WW-Variable in {filename} gefunden. Variablen: {list(ds.data_vars)}")
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

    fig, ax = plt.subplots(figsize=(12, 12), subplot_kw={"projection": ccrs.PlateCarree()})

    if var_type == "t2m":
        im = ax.pcolormesh(lon, lat, data, cmap=t2m_cmap, norm=t2m_norm, shading="auto")
    else:
        # WW-Codes filtern
        valid_mask = np.isfinite(data)
        present_codes = np.unique(data[valid_mask]).astype(int).tolist()
        present_codes = [c for c in present_codes if c not in ignore_codes and c in ww_colors_base]
        present_codes.sort()
        print(f"{filename} - gefundene WW-Codes (gefiltert): {present_codes}")

        # Farben + Mapping
        colors = [ww_colors_base[c] for c in present_codes]
        cmap = ListedColormap(colors)
        code2idx = {code: i for i, code in enumerate(present_codes)}

        # Index-Array (NaN für ignorierte Codes)
        idx_data = np.full_like(data, fill_value=np.nan, dtype=float)
        for code, idx in code2idx.items():
            idx_data[data == code] = idx

        im = ax.pcolormesh(lon, lat, idx_data, cmap=cmap,
                           vmin=-0.5, vmax=len(colors)-0.5, shading="auto")

    try:
        germany_bounds = bundeslaender.total_bounds
        ax.set_extent([germany_bounds[0]-1, germany_bounds[2]+1,
                       germany_bounds[1]-1, germany_bounds[3]+1])
    except Exception:
        pass

    bundeslaender.boundary.plot(ax=ax, edgecolor="black", linewidth=1)

    for _, city in cities.iterrows():
        ax.plot(city["lon"], city["lat"], "ko", markersize=4)
        ax.text(city["lon"] + 0.1, city["lat"] + 0.1, city["name"], fontsize=8)

    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE)

    if var_type == "t2m":
        cbar = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.05, aspect=60)
        cbar.set_ticks(list(range(-30, 45, 5)))
        cbar.set_label("Temperatur 2m [°C]", color="black")
        cbar.ax.tick_params(colors="black", labelsize=8)
        cbar.outline.set_edgecolor("black")
        cbar.ax.set_facecolor("white")
    else:
        handles = []
        for code in present_codes:
            label = f"{code}: {ww_labels.get(code, '')}".strip()
            if not label:
                continue
            color = ww_colors_base.get(code, "grey")
            handles.append(mpatches.Patch(color=color, label=label))
        ax.legend(handles=handles, loc="lower center",
                  bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=8)

    ax.set_title(
        f"ICON-D2 {var_type.upper()} Vorhersage\n"
        f"Lauf: {pd.to_datetime(run_time_utc).hour if run_time_utc is not None else '??'}Z, "
        f"gültig: {valid_time_local:%d.%m.%Y %H:%M} Uhr (MEZ/MESZ)"
    )

    outname = f"{var_type}_{pd.to_datetime(valid_time_utc):%Y%m%d_%H%M}.png"
    plt.savefig(os.path.join(output_dir, outname), dpi=150, bbox_inches="tight")
    plt.close()
