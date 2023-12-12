
###############################################################################
# LICENSE
#Copyright (C) 2018 - INPE - NATIONAL INSTITUTE FOR SPACE RESEARCH
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
###############################################################################
#======================================================================================================
# GNC-A Blog Python Tutorial: Derived Motion Winds
#======================================================================================================

# Required libraries ==================================================================================
import matplotlib.pyplot as plt                         # Import the Matplotlib package
from mpl_toolkits.basemap import Basemap                # Import the Basemap toolkit 
import numpy as np                                      # Import the Numpy package
from utilities import remap                             # Import the Remap function  
# from cpt_convert import loadCPT                       # Import the CPT convert function
from matplotlib.colors import LinearSegmentedColormap   # Linear interpolation for color maps
import datetime                                         # Library to convert julian day to dd-mm-yyyy
from netCDF4 import Dataset                             # Import the NetCDF Python interface
from utilities import download_CMI, download_DMW,loadCPT,reprojectBruno

import cartopy, cartopy.crs as ccrs  # Plot maps
import cartopy.io.shapereader as shpreader  # Import shapefiles
import cartopy.feature as cfeature # features

#======================================================================================================
# Input and output directories
## Criar diretorio dir_in band02, dmw
dir_in = "G:\Meu Drive\Minhas_coisas\Faculdade\ProjetoSAE\GOES\clone_dir_servidor\goes\\"
#Criar diretorio Scripts >> output,shapefiles,colortables,logos
dir_main = "G:\Meu Drive\Minhas_coisas\Faculdade\ProjetoSAE\GOES\clone_dir_servidor\Scripts\goes\\"
dir_out = dir_main + "output\\"
dir_libs = dir_main + "libs\\"
dir_shapefiles = dir_main + "shapefiles\\"
dir_colortables = dir_main + "colortables\\"
dir_logos = dir_main + "logos\\"
# # Desired extent
extent = [-90.0, -40.0, -20.0, 10.0]  # Max Lat, Max lon, Min lat, Min Lon

yyyymmddhhmn = '202312111200'
# Initial time and date
yyyy = datetime.datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%Y')
mm = datetime.datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%m')
dd = datetime.datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%d')
hh = datetime.datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%H')
mn = datetime.datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%M')

date_ini = str(datetime.datetime(int(yyyy), int(mm), int(dd), int(hh), int(mn)))


print(date_ini)
# Load the Data =======================================================================================
# Path to the GOES-16 image file
yyyymmddhhmn = datetime.datetime.strptime(date_ini, '%Y-%m-%d %H:%M:%S').strftime('%Y%m%d%H%M')

path = download_CMI(yyyymmddhhmn,2,dir_in+"band02")
path_reproject = reprojectBruno(f'{dir_in+"band02"}\\{path}','CMI',extent,2,f"{dir_in}band02\\")

# Path to the GOES-16 derived motion winds and file reading: 
path_dmwf = download_DMW(yyyymmddhhmn,2, dir_in+"dmw")

# Getting information from the file name ==============================================================
# Search for the Scan start in the file name
Start = (path[path.find("_s")+2:path.find("_e")])

# Search for the GOES-16 channel in the file name
Band = int((path[path.find("M6")+3:path.find("_G16")]))
# Create a GOES-16 Bands string array
Wavelenghts = ['[]','[0.47 μm]','[0.64 μm]','[0.865 μm]','[1.378 μm]','[1.61 μm]','[2.25 μm]','[3.90 μm]','[6.19 μm]','[6.95 μm]','[7.34 μm]','[8.50 μm]','[9.61 μm]','[10.35 μm]','[11.20 μm]','[12.30 μm]','[13.30 μm]']
                
# Converting from julian day to dd-mm-yyyy
year = int(Start[0:4])
dayjulian = int(Start[4:7]) - 1 # Subtract 1 because the year starts at "0"
dayconventional = datetime.datetime(year,1,1) + datetime.timedelta(dayjulian) # Convert from julian to conventional
date = dayconventional.strftime('%d-%b-%Y')  # Format the date according to the strftime directives

time = Start [7:9] + ":" + Start [9:11] + ":" + Start [11:13] + " UTC" # Time of the Start of the Scan

# Get the unit based on the channel. If channels 1 trough 6 is Albedo. If channels 7 to 16 is BT.
if Band <= 6:
    Unit = "Reflectance"
    data = Dataset(path_reproject).variables['Band1'][:] 
else:
    Unit = "Brightness Temperature [°C]"
    data = Dataset(path_reproject).variables['Band1'][:] - 273.15
    

# Choose a title for the plot
Title = " GOES-16 ABI CMI Band " + str(Band) + "       " +  Wavelenghts[int(Band)] + "       " + Unit + "       " + date + "       " + time
# Insert the institution name
Institution = "GNC-A Blog"
# =====================================================================================================


#======================================================================================================

# Define the size of the saved picture=================================================================
# Define the size of the saved picture=================================================================
#print (data.shape)
DPI = 150 
fig = plt.figure(figsize=(data.shape[1]/float(DPI), data.shape[0]/float(DPI)), frameon=False, dpi=DPI)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax = plt.axis('off')
#======================================================================================================

# Plot the Data =======================================================================================
# Create the basemap reference for the Rectangular Projection
bmap = Basemap(llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3], epsg=4326)

if Band <= 6:
    # Converts a CPT file to be used in Python
    cpt = loadCPT(f'{dir_colortables}\\Square Root Visible Enhancement.cpt')
    # Makes a linear interpolation
    cpt_convert = LinearSegmentedColormap('cpt', cpt)
    # Plot the GOES-16 channel with the converted CPT colors (you may alter the min and max to match your preference)
    bmap.imshow(data, origin='upper', cmap=cpt_convert, vmin=0, vmax=1, alpha = 1.0)  
    shapefile = "white"
    # Insert the colorbar at the bottom
elif Band == 7:
    # Converts a CPT file to be used in Python
    cpt = loadCPT(f'{dir_colortables}\\SVGAIR2_TEMP.cpt')
    # Makes a linear interpolation
    cpt_convert = LinearSegmentedColormap('cpt', cpt) 
    # Plot the GOES-16 channel with the converted CPT colors (you may alter the min and max to match your preference)
    #bmap.imshow(data, origin='upper', cmap=cpt_convert, vmin=-112.15, vmax=56.85, alpha = 1.0) 
    bmap.imshow(data, origin='upper', cmap='gray_r', vmin=-93.15, vmax=46.85, alpha = 1.0)
    shapefile = "white"
    # Insert the colorbar at the bottom
elif Band > 7 and Band < 11:
    # Converts a CPT file to be used in Python
    cpt = loadCPT(f'{dir_colortables}\\SVGAWVX_TEMP.cpt')
    # Makes a linear interpolation
    cpt_convert = LinearSegmentedColormap('cpt', cpt) 
    # Plot the GOES-16 channel with the converted CPT colors (you may alter the min and max to match your preference)
    bmap.imshow(data, origin='upper', cmap=cpt_convert, vmin=-112.15, vmax=56.85, alpha = 1.0)
    shapefile = "cyan"
    # Insert the colorbar at the bottom
elif Band > 10:
    # Converts a CPT file to be used in Python
    cpt = loadCPT(f'{dir_colortables}\\IR4AVHRR6.cpt')   
    # Makes a linear interpolation
    cpt_convert = LinearSegmentedColormap('cpt', cpt) 
    # Plot the GOES-16 channel with the converted CPT colors (you may alter the min and max to match your preference)
    #bmap.imshow(data, origin='upper', cmap=cpt_convert, vmin=-103, vmax=84, alpha = 1.0)
    bmap.imshow(data, origin='upper', cmap='gray_r', vmin=-93.15, vmax=46.85, alpha = 1.0)
    shapefile = "white"
    
# Date as string

# Time (UTC) as string
time_save = Start [7:9] + Start [9:11]

# Required libraries ==========================================================
import matplotlib.pyplot as plt          # Import the Matplotlib package
from netCDF4 import Dataset              # Import the NetCDF Python interface
import numpy as np                       # Import the Numpy package
from mpl_toolkits.basemap import Basemap # Import the Basemap toolkit 
import math                              # Import the Math package
#import time as t                         # Import the Time package
import cartopy, cartopy.crs as ccrs  # Plot maps
import cartopy.io.shapereader as shpreader  # Import shapefiles
import cartopy.feature as cfeature # features

# Opening the NetCDF DMW
nc = Dataset(f'{dir_in+"dmw"}\\{path_dmwf}')

# Read the required variables: ================================================ 
pressure = nc.variables['pressure'][:]
temperature = nc.variables['temperature'][:]
wind_direction = nc.variables['wind_direction'][:]
wind_speed = nc.variables['wind_speed'][:]
lats = nc.variables['lat'][:]
lons = nc.variables['lon'][:]

# Selecting data only from the region of interest: ============================
# Detect Latitude lower and upper index, according to the selected extent: 
latli = np.argmin( np.abs( lats - extent[1] ) ) # Lower index
latui = np.argmin( np.abs( lats - extent[3] ) ) # Upper index

# Detect the Longitude index:
# Store the indexes where the lons are between the selected extent:
lon_ind = np.where(( lons >= extent[0]) & (lons <= extent[2] ))[0]
# Eliminate the lon indexes where we don't have the lat indexes:
lon_ind = lon_ind[(lon_ind >= latui) & (lon_ind <= latli)]

# Create the variables lists ==================================================
pressure_a = []
temperature_a = []
wind_direction_a = []
wind_speed_a = []
lats_a = []
lons_a = []

# For each item, append the values to the respective variables ================
for item in lon_ind:
    lons_a.append(lons[item])
    lats_a.append(lats[item])
    pressure_a.append(pressure[item])
    temperature_a.append(temperature[item])
    wind_direction_a.append(wind_direction[item])
    wind_speed_a.append(wind_speed[item])

# Read the variables as numpy arrays
temperature = np.asarray(temperature_a)
wind_direction = np.asarray(wind_direction_a)
wind_speed = np.asarray(wind_speed_a)
lons = np.asarray(lons_a)
lats = np.asarray(lats_a)

# If you want to measure time, uncomment the next line
#start = t.time()
        
for x in range(1, 7):
        
    # Read the pressures as python arrays
    pressure = np.asarray(pressure_a)
    
    # Plot the wind vectors divided in 6 pressure ranges, separated by color
    if (x == 1): 
        print ("Plotting Pressure Range 1: (249-100 hPa)")
        pressure_index = np.where(( pressure >= 100 ) & ( pressure <= 249 ))[0]
        color = '#0000FF' # Blue 
    elif (x == 2):
        print ("Plotting Range 2: (399-250 hPa)")
        pressure_index = np.where(( pressure >= 250 ) & ( pressure <= 399 ))[0]
        color = '#309AFF' # Light Blue
    elif (x == 3):
        print ("Plotting Range 3: (400-549 hPa)")
        pressure_index = np.where(( pressure >= 400 ) & ( pressure <= 549 ))[0]
        color = '#00FF00' # Green
    elif (x == 4):
        print ("Plotting Range 4: (699-550 hPa)")
        pressure_index = np.where(( pressure >= 550 ) & ( pressure <= 699 ))[0]
        color = '#FFFF00' # Yellow
    elif (x == 5):
        print ("Plotting Range 5: (849-700 hPa)")
        pressure_index = np.where(( pressure >= 700 ) & ( pressure <= 849 ))[0]
        color = '#FF0000' # Red
    elif (x == 6):
        print ("Plotting Range 6: (1000-850 hPa)")
        pressure_index = np.where(( pressure >= 850 ) & ( pressure <= 1000 ))[0]
        color = '#FF2FCD' # Violet   
    
    # Create the variables lists (considerign only the given pressure range)
    pressure_b = []
    temperature_b = []
    wind_direction_b = []
    wind_speed_b = []
    lats_b = []
    lons_b = []

    # For each item, append the values to the respective variables 
    for item in pressure_index:
        lons_b.append(lons_a[item])
        lats_b.append(lats_a[item])
        pressure_b.append(pressure_a[item])
        temperature_b.append(temperature_a[item])
        wind_direction_b.append(wind_direction_a[item])
        wind_speed_b.append(wind_speed_a[item])
        
    # Final variables for the given pressure range
    # Read the variables as numpy arrays
    pressure = np.asarray(pressure_b)
    temperature = np.asarray(temperature_b)
    wind_direction = np.asarray(wind_direction_b)
    wind_speed = np.asarray(wind_speed_b)
    lons = np.asarray(lons_b)
    lats = np.asarray(lats_b)
        
    # Calculating the u and v components using the wind_speed and wind direction
    # in order to plot the barbs. Reference:
    # https://earthscience.stackexchange.com/questions/11982/plotting-wind-barbs-in-python-no-u-v-component
    u = []
    v = []
    for item in range(lons.shape[0]):
        u.append(-(wind_speed[item]) * math.sin((math.pi / 180) * wind_direction[item]))
        v.append(-(wind_speed[item]) * math.cos((math.pi / 180) * wind_direction[item]))
    
    # Read the u and v components as numpy arrays
    u_comp = np.asarray(u) 
    v_comp = np.asarray(v)
    
    # Create the Basemap
    bmap = Basemap(llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3], epsg=4326)
    x,y = bmap(lons, lats)
    # Make the barb plot  
    bmap.barbs(x, y, u_comp, v_comp, length=7, pivot='middle', barbcolor=color)

# If you want to measure the time required to process, uncomment the next line
#print ('- Finished! Time required:', t.time() - start, 'seconds')
    
# Add the brazilian states shapefile
#bmap.readshapefile('E:\\VLAB\\Python\\Shapefiles\\estados_2010','estados_2010',linewidth=1.00,color='khaki')
bmap.readshapefile(f'{dir_shapefiles}divisao_estados\\gadm40_BRA_1','gadm40_BRA_1', linewidth=0.65, color=shapefile)
# Add the countries shapefile

# bmap.readshapefile(f'{dir_shapefiles}america_latina\\ne_110m_admin_0_countries','ne_10m_admin_0_countries',linewidth=0.4,color=shapefile)

# Draw parallels and meridians
bmap.drawparallels(np.arange(-90.0, 90.0, 5.0), linewidth=0.25, dashes=[5, 5], color='white', labels=[False,False,False,False], fmt='%g', labelstyle="+/-", xoffset=-0.80, yoffset=-1.00, size=7)
bmap.drawmeridians(np.arange(0.0, 360.0, 5.0), linewidth=0.25, dashes=[5, 5], color='white', labels=[False,False,False,False], fmt='%g', labelstyle="+/-", xoffset=-0.80, yoffset=-1.00, size=7)


# Insert the legend
plt.text(extent[0] + 0.5,extent[1] + 6.6,'249-100 hPa', fontsize = 25,color='#0000FF')
plt.text(extent[0] + 0.5,extent[1] + 5.5,'399-250 hPa', fontsize = 25,color='#309AFF')
plt.text(extent[0] + 0.5,extent[1] + 4.5,'400-549 hPa', fontsize = 25,color='#00FF00')
plt.text(extent[0] + 0.5,extent[1] + 3.5,'699-550 hPa', fontsize = 25,color='#FFFF00')
plt.text(extent[0] + 0.5,extent[1] + 2.5,'849-700 hPa', fontsize = 25,color='#FF0000')
plt.text(extent[0] + 0.5,extent[1] + 1.5,'1000-850 hPa', fontsize = 25,color='#FF2FCD')

Title = "GOES-16 ABI CMI Band " + str(Band) + " " + Wavelenghts[int(Band)] + " " + Unit + " " + date + " " + time
plt.text(extent[0] + 0.5,extent[1] + 0.5,Title, fontsize = 25,color='black')

# Save the result

date_save = datetime.datetime.strptime(date_ini, '%Y-%m-%d %H:%M:%S').strftime('%Y%m%d%H%M')
plt.savefig(f'{dir_out}dmw\\'+ 'dmw' + str(Band) + '_' + date_save + '.png', dpi=DPI, pad_inches=0,bbox_inches='tight')
# plt.show()
plt.close()

