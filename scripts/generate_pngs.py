import sys
import cfgrib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import pandas as pd
import os
from datetime import datetime, timedelta
from zoneinfo import ZoneInfo

data_dir = sys.argv[1]
output_dir = sys.argv[2]

os.makedirs(output_dir, exist_ok=True)

# Bundesländergrenzen
bundeslaender = gpd.read_file("scripts/bundeslaender.geojson")

# Städte
cities = pd.DataFrame({
    'name': ['Berlin','Hamburg','München','Köln','Frankfurt','Dresden','Stuttgart','Düsseldorf',
             'Bremen','Kiel','Rostock','Mainz','Karlsruhe'],
    'lat': [52.52,53.55,48.14,50.94,50.11,51.05,48.78,51.23,53.08,54.32,54.09,50.00,49.01],
    'lon': [13.40,9.99,11.57,6.96,8.68,13.73,9.18,6.78,8.80,10.13,12.10,8.27,8.40]
})

# Liste aller Run-Zeiten in UTC
RUNS_UTC = [0, 3, 6, 9, 12, 15, 18, 21]

# Durchlaufe alle Runs
for run_hour in RUNS_UTC:
    run_dir = os.path.join(data_dir, f"{run_hour:02d}")  # z.B. "21" für 21z
    if not os.path.exists(run_dir):
        continue  # überspringen, falls Ordner noch nicht existiert

    grib_files = sorted([f for f in os.listdir(run_dir) if f.endswith(".grib2")])
    
    for idx, filename in enumerate(grib_files):
        step_hour = idx  # Forecast-Stunde nach Run
        
        ds = cfgrib.open_dataset(os.path.join(run_dir, filename))
        t2m = ds['t2m'] - 273.15
        if len(t2m.shape) == 3:
            t2m = t2m[0,:,:]

        lon = ds['longitude']
        lat = ds['latitude']

        fig, ax = plt.subplots(figsize=(10,10), subplot_kw={'projection': ccrs.PlateCarree()})
        ax.set_extent([5,16,47,56])

        im = ax.pcolormesh(lon, lat, t2m, cmap='RdYlBu_r', shading='auto', vmin=-20, vmax=35)

        bundeslaender.boundary.plot(ax=ax, edgecolor='black', linewidth=1)
        for _, city in cities.iterrows():
            ax.plot(city['lon'], city['lat'], 'ko', markersize=4)
            ax.text(city['lon']+0.1, city['lat']+0.1, city['name'], fontsize=8)

        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.COASTLINE)

        cbar = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, label='Temperatur 2m [°C]')

        # Forecast-Zeit in lokale Zeit berechnen
        today = datetime.utcnow().date()
        forecast_utc = datetime.combine(today, datetime.min.time()) + timedelta(hours=run_hour + step_hour)
        forecast_local = forecast_utc.astimezone(ZoneInfo("Europe/Berlin"))
        forecast_str = forecast_local.strftime("%d.%m.%Y %H:%M %Z")

        ax.set_title(f"ICON-D2 2m Temperatur - {forecast_str}")

        # Speichern
        out_dir = os.path.join(output_dir, f"{run_hour:02d}/t_2m")
        os.makedirs(out_dir, exist_ok=True)
        plt.savefig(os.path.join(out_dir, f"output_{step_hour:03d}.png"), dpi=150)
        plt.close()
