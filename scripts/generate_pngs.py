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
import re

# Eingabe-/Ausgabe-Verzeichnisse
data_dir = sys.argv[1]
output_dir = sys.argv[2]

os.makedirs(output_dir, exist_ok=True)

bundeslaender = gpd.read_file("scripts/bundeslaender.geojson")

cities = pd.DataFrame({
    'name': ['Berlin','Hamburg','München','Köln','Frankfurt','Dresden','Stuttgart','Düsseldorf'],
    'lat': [52.52,53.55,48.14,50.94,50.11,51.05,48.78,51.23],
    'lon': [13.40,9.99,11.57,6.96,8.68,13.73,9.18,6.78]
})

bounds = list(range(-30, 45, 5))

colors = [
    "#001070", "#0020c2", "#0040ff", "#0080ff", "#00c0ff", "#00ffff",
    "#80ff80", "#c0ff00", "#ffff00", "#ffcc00", "#ff8000", "#ff4000",
    "#ff0000", "#990000",
]

cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(bounds, cmap.N)

for filename in sorted(os.listdir(data_dir)):
    if not filename.endswith(".grib2"):
        continue

    ds = cfgrib.open_dataset(os.path.join(data_dir, filename))
    t2m = ds['t2m'] - 273.15

    if len(t2m.shape) == 3:
        t2m = t2m[0, :, :]

    lon = ds['longitude']
    lat = ds['latitude']

    valid_time_utc = pd.to_datetime(ds.valid_time.values).tz_localize("UTC")
    valid_time_local = valid_time_utc.astimezone(ZoneInfo("Europe/Berlin"))

    # Modelllauf aus Dateiname extrahieren
    match = re.search(r'\((\d{2})', filename)
    lauf_str = match.group(1) + "z" if match else "??z"

    fig, ax = plt.subplots(figsize=(20, 10), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([5, 16, 47, 56])

    im = ax.pcolormesh(lon, lat, t2m, cmap=cmap, norm=norm, shading='auto')

    bundeslaender.boundary.plot(ax=ax, edgecolor='black', linewidth=1)

    for idx, city in cities.iterrows():
        ax.plot(city['lon'], city['lat'], 'ko', markersize=4)
        ax.text(city['lon'] + 0.1, city['lat'] + 0.1, city['name'], fontsize=8)

    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.COASTLINE)

    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=60, ticks=bounds)
    cbar.set_label("Temperatur 2m [°C]", color="black")
    cbar.ax.tick_params(colors="black", labelsize=8)
    cbar.outline.set_edgecolor("black")
    cbar.ax.set_facecolor("white")

    ax.set_title(f"ICON-D2 2m Temperatur - {valid_time_local:%d.%m.%Y %H:%M} Uhr ({lauf_str})")

    outname = f"output_{valid_time_utc:%Y%m%d_%H%M}.png"
    plt.savefig(os.path.join(output_dir, outname), dpi=150, bbox_inches="tight")
    plt.close()
