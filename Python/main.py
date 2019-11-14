#!/usr/bin/env python

from read_netcdf_general import read_nc_data
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs

data_dir = '/Users/christoph/Desktop/DAMOCLES_training_school/WorkingGroup1/nc_data/'
file_name = 'NH_yield_and_meteovar.nc'
var_list = ['lon', 'lat', 'year', 'yield', 'month', 'vpd', 'tasmax', 'pr']
data = read_nc_data(data_dir, file_name, var_list)

lon = data['lon']
lat = data['lat']

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent([-180, 180, 0, 90], crs=ccrs.PlateCarree())
ax.coastlines()
yd_mean = np.mean(data['yield'], axis=0)
yd_max = np.max(yd_mean)

plt.scatter(lon, lat, c=yd_mean, vmin=0, vmax=yd_max, s=1.5, cmap=cm, marker='s')
cbar = plt.colorbar(orientation='horizontal', fraction=0.09, pad=0.04,)
cbar.ax.tick_params(labelsize=12)
plt.title('Average crop yield [kg] over all years', fontsize=18)
plt.tight_layout()
fig.savefig('av_yield_global.pdf')
plt.show()

