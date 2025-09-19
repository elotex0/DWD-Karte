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

# Deine Wunschfarben (nur die, die du festgelegt hast)
ww_colors_base = {
    0:"#FFFFFF", 1:"#D3D3D3", 2:"#A9A9A9", 3:"#696969",     # Bewölkung
    45:"#FFFF00", 48:"#FFD700",                               # Nebel
    56:"#FFA500", 57:"#FF8C00",                               # Schneeregen (du hattest diese)
    51:"#FFA07A", 53:"#FF8C69", 55:"#FF7043",                 # fallback für Nieselregen (eigene Wahl)
    61:"#90EE90", 63:"#32CD32", 65:"#006400",                # Regen
    66:"#FF6347", 67:"#8B0000",                               # gef. Regen
    71:"#ADD8E6", 73:"#6495ED", 75:"#00008B",                # Schneefall
    80:"#87CEFA", 81:"#4682B4",                              # Regenschauer (ergänzt)
    95:"#FF77FF", 96:"#C71585"                                # Gewitter
}

# Labeltexte (kurz) für Legende
ww_labels = {
    0:"klar",1:"leicht bewölkt",2:"teilweise bewölkt",3:"bedeckt",
    45:"Nebel",48:"Nebel mit Reif",
    51:"Niesel leicht",53:"Niesel mäßig",55:"Niesel stark",
    56:"Schneeregen leicht",57:"Schneeregen stark",
    61:"Regen leicht",63:"Regen mäßig",65:"Regen stark",
    66:"gef. Regen leicht",67:"gef. Regen stark",
    71:"Schnee leicht",73:"Schnee mäßig",75:"Schnee stark",
    80:"Schauer leicht",81:"Schauer mäßig",
    95:"Gewitter leicht/mäßig",96:"Gewitter stark"
}

for filename in sorted(os.listdir(data_dir)):
    if not filename.endswith(".grib2"):
        continue

    path = os.path.join(data_dir, filename)
    ds = cfgrib.open_dataset(path)

    # Variable finden (robust)
    varname = None
    if "WW" in ds:
        varname = "WW"
    elif "ww" in ds:
        varname = "ww"
    else:
        # fallback: finde den DataVar-Namen, der 'WW' enthält (case-insensitive)
        for vn in ds.data_vars:
            if vn.lower() == "ww" or "weather" in vn.lower():
                varname = vn
                break

    if varname is None:
        print(f"Keine WW-Variable in {filename} gefunden. Variablen: {list(ds.data_vars)}")
        continue

    data = ds[varname].values
    # 3D -> erste Zeitschicht
    if data.ndim == 3:
        data = data[0, :, :]

    lon = ds['longitude'].values
    lat = ds['latitude'].values

    # Bestimme alle echten Codes in der Datei (int)
    valid_mask = np.isfinite(data)
    present_codes = np.unique(data[valid_mask]).astype(int).tolist()
    present_codes.sort()
    print(f"{filename} - gefundene WW-Codes: {present_codes}")

    # Baue Farb-Liste: verwende deine festen Farben, andere automatisch aus einer colormap
    colors = []
    codes_sorted = present_codes
    # fallback colormap
    fallback_cmap = plt.get_cmap("tab20")
    fallback_i = 0
    for code in codes_sorted:
        if code in ww_colors_base:
            colors.append(ww_colors_base[code])
        else:
            # nimm nächste Farbe aus fallback (hex)
            colors.append(mcolors.to_hex(fallback_cmap(fallback_i)))
            fallback_i += 1

    cmap = ListedColormap(colors)

    # Mapping code -> index (0..N-1)
    code2idx = {code: i for i, code in enumerate(codes_sorted)}
    # Erzeuge Index-Array
    idx_data = np.full_like(data, fill_value=np.nan, dtype=float)
    for code, idx in code2idx.items():
        idx_data[data == code] = idx

    # Plot
    run_time_utc = pd.to_datetime(ds['time'].values) if 'time' in ds else None
    valid_time_utc = pd.to_datetime(ds.valid_time.values) if 'valid_time' in ds else None
    if isinstance(valid_time_utc, (np.ndarray, list)):
        valid_time_utc = valid_time_utc[0]
    valid_time_local = pd.to_datetime(valid_time_utc).tz_localize("UTC").astimezone(ZoneInfo("Europe/Berlin"))

    fig, ax = plt.subplots(figsize=(20,10), subplot_kw={'projection': ccrs.PlateCarree()})

    # Indizes so plotten, dass jede Indexnummer genau eine Farbe bekommt:
    im = ax.pcolormesh(lon, lat, idx_data, cmap=cmap, vmin=-0.5, vmax=len(colors)-0.5, shading='auto')

    germany_bounds = bundeslaender.total_bounds  # Achtung: evtl Tippfehler reparieren wenn variable anders heißt
    try:
        germany_bounds = bundeslaender.total_bounds
        ax.set_extent([germany_bounds[0]-1, germany_bounds[2]+1, germany_bounds[1]-1, germany_bounds[3]+1])
    except Exception:
        pass

    bundeslaender.boundary.plot(ax=ax, edgecolor='black', linewidth=1)

    for idx_city, city in cities.iterrows():
        ax.plot(city['lon'], city['lat'], 'ko', markersize=4)
        ax.text(city['lon'] + 0.1, city['lat'] + 0.1, city['name'], fontsize=8)

    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.COASTLINE)

    # Legende als farbige Patches
    handles = []
    labels = []
    for code in codes_sorted:
        color = colors[code2idx[code]]
        label = f"{code}: {ww_labels.get(code, '')}".strip()
        handles.append(mpatches.Patch(color=color, label=label))
    # Platziere die Legende unter der Karte (anpassen ncol falls viele Einträge)
    ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=8)

    ax.set_title(f"ICON-D2 WW Vorhersage\nLauf: {pd.to_datetime(run_time_utc).hour if run_time_utc is not None else '??'}Z, gültig: {valid_time_local:%d.%m.%Y %H:%M} Uhr (MEZ/MESZ)")

    outname = f"ww_{pd.to_datetime(valid_time_utc):%Y%m%d_%H%M}.png"
    plt.savefig(os.path.join(output_dir, outname), dpi=150, bbox_inches="tight")
    plt.close()
