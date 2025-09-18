import sys
import cfgrib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import pandas as pd
import os

data_dir = sys.argv[1]
output_dir = sys.argv[2]

# Lade Bundesländergrenzen
bundeslaender = gpd.read_file("scripts/bundeslaender.geojson")

# Lade einige Städte in Deutschland
cities = pd.DataFrame({
    'name': ['Berlin','Hamburg','München','Köln','Frankfurt','Dresden','Stuttgart','Düsseldorf'],
    'lat': [52.52,53.55,48.14,50.94,50.11,51.05,48.78,51.23],
    'lon': [13.40,9.99,11.57,6.96,8.68,13.73,9.18,6.78]
})

for filename in sorted(os.listdir(data_dir)):
    if not filename.endswith(".grib2"):
        continue
    step = filename.split("_")[-1].split(".")[0]
    ds = cfgrib.open_dataset(os.path.join(data_dir, filename))
    t2m = ds['t2m'] - 273.15  # K -> °C
    lon = ds['longitude']
    lat = ds['latitude']

    fig, ax = plt.subplots(figsize=(10,10), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([5,16,47,56])  # Deutschland

    # Temperatur Colormap: blau=kalt, rot=heiß
    im = ax.pcolormesh(lon, lat, t2m, cmap='coolwarm', shading='auto', vmin=-20, vmax=35)

    # Bundesländergrenzen
    bundeslaender.boundary.plot(ax=ax, edgecolor='black', linewidth=1)

    # Städte einzeichnen
    for idx, city in cities.iterrows():
        ax.plot(city['lon'], city['lat'], 'ko', markersize=4)
        ax.text(city['lon']+0.1, city['lat']+0.1, city['name'], fontsize=8)

    # Küstenlinien und Ländergrenzen
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.COASTLINE)

    # Legende
    cbar = fig.colorbar(im, ax=ax, orientation='vertical', label='Temperatur 2m [°C]')

    # Titel
    ax.set_title(f"ICON-D2 2m Temperatur - Schritt {step}")

    plt.savefig(os.path.join(output_dir, f"output_{step}.png"), dpi=150)
    plt.close()
