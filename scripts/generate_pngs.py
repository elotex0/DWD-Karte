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

# Ordner sicher erstellen
os.makedirs(output_dir, exist_ok=True)

# Lade Bundesländergrenzen
bundeslaender = gpd.read_file("scripts/bundeslaender.geojson")

# Einige Städte in Deutschland
cities = pd.DataFrame({
    'name': ['Berlin', 'Hamburg', 'München', 'Köln', 'Frankfurt', 'Dresden', 'Stuttgart', 'Düsseldorf'],
    'lat': [52.52, 53.55, 48.14, 50.94, 50.11, 51.05, 48.78, 51.23],
    'lon': [13.40, 9.99, 11.57, 6.96, 8.68, 13.73, 9.18, 6.78]
})

for filename in sorted(os.listdir(data_dir)):
    if not filename.endswith(".grib2"):
        continue
   
    # Lade die GRIB-Datei
    ds = cfgrib.open_dataset(os.path.join(data_dir, filename))
   
    # Extrahiere die gültige Zeit (valid_time) der Vorhersage
    if 'valid_time' in ds:
        valid_time = ds.valid_time.values
        # Konvertiere valid_time (UTC) nach lokaler Zeit (MESZ/MEZ)
        local_time = pd.Timestamp(valid_time).tz_localize('UTC').tz_convert('Europe/Berlin')
        # Bestimme, ob Sommerzeit (MESZ) oder Winterzeit (MEZ)
        tz_name = 'MESZ' if local_time.dst() != pd.Timedelta(0) else 'MEZ'
        # Formatiere die Zeit
        time_str = f"{local_time.strftime('%Y-%m-%d %H:%M')} {tz_name}"
    else:
        # Fallback: Wenn valid_time nicht verfügbar ist
        print(f"Warnung: valid_time nicht in {filename} gefunden. Verwende Dateinamen für Zeit.")
        step = int(filename.split("_")[-1].split(".")[0])  # z. B. '000' -> 0, '003' -> 3
        model_run = filename.split("_")[-2]  # z. B. '00z', '03z'
        model_hour = int(model_run[:-1])  # z. B. '00' -> 0, '03' -> 3
        base_time = pd.Timestamp(f"2025-09-19 {model_hour:02d}:00:00").tz_localize('UTC')
        local_time = (base_time + pd.Timedelta(hours=step)).tz_convert('Europe/Berlin')
        tz_name = 'MESZ' if local_time.dst() != pd.Timedelta(0) else 'MEZ'
        time_str = f"{local_time.strftime('%Y-%m-%d %H:%M')} {tz_name}"
   
    t2m = ds['t2m'] - 273.15  # K -> °C
   
    # Prüfen ob 2D oder 3D
    if len(t2m.shape) == 3:
        t2m = t2m[0,:,:]
   
    lon = ds['longitude']
    lat = ds['latitude']
   
    fig, ax = plt.subplots(figsize=(10,10), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([5, 16, 47, 56])  # Deutschland
   
    # Temperatur Colormap
    im = ax.pcolormesh(lon, lat, t2m, cmap='coolwarm', shading='auto', vmin=-20, vmax=35)
   
    # Bundesländer
    bundeslaender.boundary.plot(ax=ax, edgecolor='black', linewidth=1)
   
    # Städte
    for idx, city in cities.iterrows():
        ax.plot(city['lon'], city['lat'], 'ko', markersize=4)
        ax.text(city['lon']+0.1, city['lat']+0.1, city['name'], fontsize=8)
   
    # Küstenlinien & Ländergrenzen
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.COASTLINE)
   
    # Legende
    cbar = fig.colorbar(im, ax=ax, orientation='vertical', label='Temperatur 2m [°C]')
   
    # Titel mit der berechneten Uhrzeit
    ax.set_title(f"ICON-D2 2m Temperatur - {time_str}")
   
    # Speichern
    plt.savefig(os.path.join(output_dir, f"output_{filename.split('_')[-1].split('.')[0]}.png"), dpi=150)
    plt.close()
