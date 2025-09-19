import sys
import cfgrib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import pandas as pd
import os
import matplotlib.colors as mcolors
from zoneinfo import ZoneInfo
from matplotlib.patches import Patch

# Eingabe-/Ausgabe-Verzeichnisse
data_dir = sys.argv[1]
output_dir = sys.argv[2]
var_type = sys.argv[3]  # 't2m' oder 'ww'

# Ordner sicher erstellen
os.makedirs(output_dir, exist_ok=True)

# Bundesländergrenzen laden
bundeslaender = gpd.read_file("scripts/bundeslaender.geojson")

# Einige Städte in Deutschland
cities = pd.DataFrame({
    'name': ['Berlin','Hamburg','München','Köln','Frankfurt','Dresden','Stuttgart','Düsseldorf'],
    'lat': [52.52,53.55,48.14,50.94,50.11,51.05,48.78,51.23],
    'lon': [13.40,9.99,11.57,6.96,8.68,13.73,9.18,6.78]
})

if var_type == "t2m":
    # Temperaturbereiche von -30 bis +40 in 5°C Schritten
    bounds = list(range(-30, 45, 5))
    colors = [
        "#001070","#0020c2","#0040ff","#0080ff","#00c0ff","#00ffff",
        "#80ff80","#c0ff00","#ffff00","#ffcc00","#ff8000","#ff4000",
        "#ff0000","#990000"
    ]
    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
elif var_type == "WW":
    # Farben für WW-Karte (originale GRIB-Codes)
    ww_colors = {
        0:"#FFFFFF", 1:"#D3D3D3", 2:"#A9A9A9", 3:"#696969",
        45:"#FFFF00", 48:"#FFD700",
        51:"#90EE90", 61:"#32CD32", 63:"#228B22",
        65:"#FF6347", 66:"#FF0000",
        56:"#FFA500", 57:"#FF8C00",
        71:"#ADD8E6", 73:"#87CEEB", 75:"#4682B4",
        95:"#FF77FF", 96:"#C71585"
    }

    codes = sorted(ww_colors.keys())
    colors = [ww_colors[c] for c in codes]
    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(codes + [max(codes)+1], cmap.N)

    # Beschreibungen für die Legende
    legend_labels = {
        0: "klar/leicht bewölkt", 1: "bewölkt 1", 2: "bewölkt 2", 3: "bewölkt 3",
        45: "Fog", 48: "Fog Reifbildung",
        51: "Regen leicht", 61: "Regen mäßig", 63: "Regen stark",
        65: "Gef. Regen leicht", 66: "Gef. Regen stark",
        56: "Schneeregen leicht", 57: "Schneeregen stark",
        71: "Schneefall leicht", 73: "Schneefall mäßig", 75: "Schneefall stark",
        95: "Gewitter leicht/mäßig", 96: "Gewitter stark"
    }

# Loop über alle GRIB2-Dateien
for filename in sorted(os.listdir(data_dir)):
    if not filename.endswith(".grib2"):
        continue

    ds = cfgrib.open_dataset(os.path.join(data_dir, filename))

    # Daten extrahieren
    if var_type == "t2m":
        data = ds['t2m'] - 273.15  # K → °C
    else:
        data = ds['WW']

    # Falls 3D-Feld → erste Zeitschicht nehmen
    if len(data.shape) == 3:
        data = data[0, :, :]

    lon = ds['longitude'].values
    lat = ds['latitude'].values

    # Modell-Laufzeit
    run_time_utc = pd.to_datetime(ds['time'].values, unit='s', origin='unix')
    run_hour_z = f"{run_time_utc.hour:02d}Z"

    # Vorhersagezeit lokal in Berlin
    valid_time_utc = pd.to_datetime(ds.valid_time.values).tz_localize("UTC")
    valid_time_local = valid_time_utc.astimezone(ZoneInfo("Europe/Berlin"))

    # Figur erstellen
    fig, ax = plt.subplots(figsize=(20, 10), subplot_kw={'projection': ccrs.PlateCarree()})

    # Karte plotten
    im = ax.pcolormesh(lon, lat, data, cmap=cmap, norm=norm, shading='auto')

    # Deutschland + Nachbarländer automatisch als Ausschnitt
    germany_bounds = bundeslaender.total_bounds
    ax.set_extent([
        germany_bounds[0]-1,
        germany_bounds[2]+1,
        germany_bounds[1]-1,
        germany_bounds[3]+1
    ])

    # Bundesländergrenzen
    bundeslaender.boundary.plot(ax=ax, edgecolor='black', linewidth=1)

    # Städte
    for idx, city in cities.iterrows():
        ax.plot(city['lon'], city['lat'], 'ko', markersize=4)
        ax.text(city['lon'] + 0.1, city['lat'] + 0.1, city['name'], fontsize=8)

    # Küstenlinien & Ländergrenzen
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.COASTLINE)

    # Legende
    if var_type == "t2m":
        cbar = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=60, ticks=bounds)
        cbar.set_label("Temperatur 2m [°C]", color="black")
        cbar.ax.tick_params(colors="black", labelsize=8)
    else:
        legend_handles = [Patch(facecolor=ww_colors[c], edgecolor='black', label=legend_labels[c]) for c in codes]
        ax.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, -0.15),
                  ncol=4, fontsize=8, title="Wetterphänomene")

    # Titel
    ax.set_title(
        f"ICON-D2 {var_type} Vorhersage\nLauf: {run_hour_z}, "
        f"gültig: {valid_time_local:%d.%m.%Y %H:%M} Uhr (MEZ/MESZ)"
    )

    # Speichern
    outname = f"{var_type}_{valid_time_utc:%Y%m%d_%H%M}.png"
    plt.savefig(os.path.join(output_dir, outname), dpi=150, bbox_inches="tight")
    plt.close()
