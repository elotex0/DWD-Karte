import sys
import cfgrib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import pandas as pd
import os
import matplotlib.colors as mcolors

data_dir = sys.argv[1]
output_dir = sys.argv[2]

os.makedirs(output_dir, exist_ok=True)

bundeslaender = gpd.read_file("scripts/bundeslaender.geojson")

cities = pd.DataFrame({
    'name': ['Berlin','Hamburg','München','Köln','Frankfurt','Dresden','Stuttgart','Düsseldorf'],
    'lat': [52.52,53.55,48.14,50.94,50.11,51.05,48.78,51.23],
    'lon': [13.40,9.99,11.57,6.96,8.68,13.73,9.18,6.78]
})

# Farbbereiche ähnlich DWD, von -30 bis +40 in 5° Schritten
bounds = list(range(-30, 45, 5))

# Farbpalette, angelehnt an typische DWD Modellkarten
colors = [
    "#002070",  # sehr kalt
    "#0040ff",  # kalt
    "#00a0ff",
    "#00ffff",
    "#80ff80",
    "#c0ff00",
    "#ffff00",  # moderat
    "#ffcc00",
    "#ff8000",
    "#ff4000",
    "#ff0000",  # heiß
    "#cc0000",
    "#990000",  # sehr heiß
]

cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(bounds, cmap.N)

for filename in sorted(os.listdir(data_dir)):
    if not filename.endswith(".grib2"):
        continue

    ds = cfgrib.open_dataset(os.path.join(data_dir, filename))
    t2m = ds['t2m'] - 273.15

    if len(t2m.shape) == 3:
        t2m = t2m[0,:,:]

    lon = ds['longitude']
    lat = ds['latitude']

    valid_time = pd.to_datetime(ds.valid_time.values)

    fig, ax = plt.subplots(figsize=(18,10), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([5,16,47,56])

    im = ax.pcolormesh(lon, lat, t2m, cmap=cmap, norm=norm, shading='auto')

    bundeslaender.boundary.plot(ax=ax, edgecolor='black', linewidth=1)

    for idx, city in cities.iterrows():
        ax.plot(city['lon'], city['lat'], 'ko', markersize=4)
        ax.text(city['lon']+0.1, city['lat']+0.1, city['name'], fontsize=8)

    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.COASTLINE)

    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=60, ticks=bounds)
    cbar.set_label("Temperatur 2m [°C]", color="black")
    cbar.ax.tick_params(colors="black", labelsize=8)
    cbar.outline.set_edgecolor("black")
    cbar.ax.set_facecolor("white")  # oder hellem Grau, damit die Zahlen schön sichtbar sind

    ax.set_title(f"ICON-D2 2m Temperatur - {valid_time:%d.%m.%Y %H:%M} UTC")

    outname = f"output_{valid_time:%Y%m%d_%H%M}.png"
    plt.savefig(os.path.join(output_dir, outname), dpi=150, bbox_inches="tight")
    plt.close()
