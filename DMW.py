# Required libraries ==================================================================================
import matplotlib.pyplot as plt                         # Import the Matplotlib package
from mpl_toolkits.basemap import Basemap                # Import the Basemap toolkit 
import numpy as np                                      # Import the Numpy package
from utilities import remap                             # Import the Remap fudmwtion  
# from cpt_convert import loadCPT                       # Import the CPT convert fudmwtion
from matplotlib.colors import LinearSegmentedColormap   # Linear interpolation for color maps
import datetime                                         # Library to convert julian day to dd-mm-yyyy
from netCDF4 import Dataset                             # Import the NetCDF Python interface
from utilities import download_CMI, download_DMW,loadCPT,reprojectBruno

import cartopy, cartopy.crs as ccrs  # Plot maps
import cartopy.io.shapereader as shpreader  # Import color_shapefiles
import cartopy.feature as cfeature # features


def process_dmw(file_dmw,file_fundo,v_extent):
    global dir_maps
    start = time.time()  
    #---------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------
    # Area de interesse para recorte
    extent, resolution = area_para_recorte(v_extent)
    variable = "CMI"

    # Getting information from the file name ==============================================================
    # Search for the Scan start in the file name
    Start = (file_fundo[file_fundo.find("_s")+2:file_fundo.find("_e")])

    # Search for the GOES-16 channel in the file name
    Band = int((file_fundo[file_fundo.find("M6")+3:file_fundo.find("_G16")]))
    # Create a GOES-16 file_fundos string array
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
        # reprojetando file fundo
        grid = remap(file_fundo, variable, extent, resolution)
        # Lê o retorno da função
        data_fundo = grid.ReadAsArray()
    else:
        Unit = "Brightness Temperature [°C]"
        # reprojetando band 01
        grid = remap(file_fundo, variable, extent, resolution)
        # Lê o retorno da função
        data_fundo = grid.ReadAsArray()- 273.15
        
    # Choose a title for the plot
    description = f"GOES-16 Derived Motion Winds - Band {str(Band)} {Wavelenghts[int(Band)]} {Unit} {date} {time}"
    # Insert the institution name
    institution = "CEPAGRI-UNICAMP"

    # Define the size of the saved picture=================================================================
    d_p_i = 150
    fig = plt.figure(figsize=(2000/ float(d_p_i), 2000/ float(d_p_i)), frameon=True, 
    dpi=d_p_i, edgecolor='black', facecolor='black')
    
    # Define the projection
    proj = ccrs.PlateCarree()
    
    # Use the PlateCarree projection in cartopy
    ax = plt.axes([0, 0, 1, 1], projection=proj)
    
    # Define the image extent
    img_extent = [extent[0], extent[2], extent[1], extent[3]]
    ax.set_extent(img_extent, ccrs.PlateCarree())
  
    # Insert the legend
    # ax.text(extent[2] - 11,extent[1] + 8.5,'249-100 hPa', fontsize = 20,color='#0000FF')
    # ax.text(extent[2] - 11,extent[1] + 7,'399-250 hPa', fontsize = 20,color='#309AFF')
    # ax.text(extent[2] - 11,extent[1] + 5.5,'400-549 hPa', fontsize = 20,color='#00FF00')
    # ax.text(extent[2] - 11,extent[1] + 4,'699-550 hPa', fontsize = 20,color='#FFFF00')
    # ax.text(extent[2] - 11,extent[1] + 2.5,'849-700 hPa', fontsize = 20,color='#FF0000')
    # ax.text(extent[2] - 11,extent[1] + 1,'1000-850 hPa', fontsize = 20,color='#FF2FCD') 
    legenda = plt.imread(dir_logos + 'dmw_legend.png')  # Lendo o arquivo do logo
    fig.figimage(legenda, 1750, 70, zorder=1, alpha=0.8, origin='upper')  # Plotando legenda

    # Plot the Data =======================================================================================
    if Band <= 6:
        ## Checar esses diretórios da colortables e ve se tem essas palhetas de cores, acredito que tenha
        # Converts a CPT file to be used in Python
        cpt = loadCPT(f'{dir_colortables}\\Square Root Visible Enhancement.cpt')
        # Makes a linear interpolation
        cpt_convert = LinearSegmentedColormap('cpt', cpt)
      
        # Plot the GOES-16 channel with the converted CPT colors (you may alter the min and max to match your preference)
        ax.imshow(data_fundo, origin='upper', cmap=cpt_convert, vmin=0, vmax=1, alpha = 1.0,extent=img_extent)  
        color_shapefile = "white"
      
        # Adicionando o shapefile dos estados brasileiros
        adicionando_shapefile(v_extent, ax,color_shapefile)

        # Adicionando  linhas dos litorais
        adicionando_linhas(ax)

        # Adicionando descricao da imagem
        adicionando_descricao_imagem(description, institution, ax, fig)

        # Adicionando os logos
        adicionando_logos(fig)
    elif Band == 7:
        # Converts a CPT file to be used in Python
        cpt = loadCPT(f'{dir_colortables}\\SVGAIR2_TEMP.cpt')
      
        # Makes a linear interpolation
        cpt_convert = LinearSegmentedColormap('cpt', cpt) 
      
        # Plot the GOES-16 channel with the converted CPT colors (you may alter the min and max to match your preference)
        ax.imshow(data_fundo, origin='upper', cmap='gray_r', vmin=-93.15, vmax=46.85, alpha = 1.0,extent=img_extent)
        color_shapefile = "white"
      
        # Adicionando o shapefile dos estados brasileiros
        adicionando_shapefile(v_extent, ax,color_shapefile)
        # Adicionando  linhas dos litorais
        adicionando_linhas(ax)
        # Adicionando descricao da imagem
        adicionando_descricao_imagem(description, institution, ax, fig)
        # Adicionando os logos
        adicionando_logos(fig)
      
    elif Band > 7 and Band < 11:
        # Converts a CPT file to be used in Python
        cpt = loadCPT(f'{dir_colortables}\\SVGAWVX_TEMP.cpt')
      
        # Makes a linear interpolation
        cpt_convert = LinearSegmentedColormap('cpt', cpt) 
      
        # Plot the GOES-16 channel with the converted CPT colors (you may alter the min and max to match your preference)
        ax.imshow(data_fundo, origin='upper', cmap=cpt_convert, vmin=-112.15, vmax=56.85, alpha = 1.0,extent=img_extent)
        color_shapefile = "cyan"
      
        # Adicionando o shapefile dos estados brasileiros
        adicionando_shapefile(v_extent, ax,color_shapefile)
        # Adicionando  linhas dos litorais
        adicionando_linhas(ax)
        # Adicionando descricao da imagem
        adicionando_descricao_imagem(description, institution, ax, fig)
        # Adicionando os logos
        adicionando_logos(fig)
      
    elif Band > 10:
        # Converts a CPT file to be used in Python
        cpt = loadCPT(f'{dir_colortables}\\IR4AVHRR6.cpt')  
      
        # Makes a linear interpolation
        cpt_convert = LinearSegmentedColormap('cpt', cpt) 
      
        # Plot the GOES-16 channel with the converted CPT colors (you may alter the min and max to match your preference)
        #ax.imshow(data, origin='upper', cmap=cpt_convert, vmin=-103, vmax=84, alpha = 1.0)
        ax.imshow(data_fundo, origin='upper', cmap='gray_r', vmin=-93.15, vmax=46.85, alpha = 1.0,extent=img_extent)
        color_shapefile = "white"
      
        # Adicionando o shapefile dos estados brasileiros
        adicionando_shapefile(v_extent, ax,color_shapefile)
        # Adicionando  linhas dos litorais
        adicionando_linhas(ax)
        # Adicionando descricao da imagem
        adicionando_descricao_imagem(description, institution, ax, fig)
        # Adicionando os logos
        adicionando_logos(fig)
    
    # Time (UTC) as string
    time_save = Start [7:9] + Start [9:11]

    # Opening the NetCDF DMW
    dmw = Dataset(file_dmw)
    # Read the required variables: ================================================ 
    pressure = dmw.variables['pressure'][:]
    temperature = dmw.variables['temperature'][:]
    wind_direction = dmw.variables['wind_direction'][:]
    wind_speed = dmw.variables['wind_speed'][:]
    lats = dmw.variables['lat'][:]
    lons = dmw.variables['lon'][:]
    
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
    pressure_a = np.take(pressure,lon_ind)
    wind_direction_a = np.take(wind_direction,lon_ind)
    wind_speed_a = np.take(wind_speed,lon_ind)
    lats_a = np.take(lats,lon_ind)
    lons_a = np.take(lons,lon_ind)
    
    # Read the variables as numpy arrays
    temperature = np.asarray(temperature_a)
    wind_direction = np.asarray(wind_direction_a)
    wind_speed = np.asarray(wind_speed_a)
    lons = np.asarray(lons_a)
    lats = np.asarray(lats_a)

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
        wind_direction_b = np.take(wind_direction_a,pressure_index)
        wind_speed_b = np.take(wind_speed_a,pressure_index)
        lats_b = np.take(lats_a,pressure_index)
        lons_b = np.take(lons_a,pressure_index)
    
            
        # Final variables for the given pressure range
        # Read the variables as numpy arrays
        wind_direction = np.asarray(wind_direction_b)
        wind_direction[wind_direction== -999]=np.nan
        #Converte a direção de graus para radianos
        wind_direction = np.deg2rad(wind_direction)
        wind_speed = np.asarray(wind_speed_b)
        wind_speed[wind_speed== -999]=np.nan
        lons = np.asarray(lons_b)
        lons[lons== -999]=np.nan
        lats = np.asarray(lats_b)
        lats[lats== -999]=np.nan
            
        # Calculating the u and v components using the wind_speed and wind direction
        # in order to plot the barbs. Referedmwe:
        # https://earthsciedmwe.stackexchange.com/questions/11982/plotting-wind-barbs-in-python-no-u-v-component
        u = []
        v = []
        for item in range(lons.shape[0]):
            u.append(-(wind_speed[item]) * np.sin(wind_direction[item]))
            v.append(-(wind_speed[item]) * np.cos(wind_direction[item]))
        
        # Read the u and v components as numpy arrays
        u_comp = np.asarray(u) 
        v_comp = np.asarray(v)
    
        ax.barbs(lons,lats,u_comp,v_comp,length=5, barbcolor=color,pivot='middle',transform=ccrs.PlateCarree())
    
    # Salvando a imagem de saida -- Alterar diretório pra salvar
    plt.savefig(f'{dir_out}{rgb_type}/{rgb_type}_{date_file}_{v_extent}.png', bbox_inches='tight', pad_inches=0, dpi=d_p_i)

    # Fecha a janela para limpar a memoria
    plt.close()
    
    logging.info(f'Derived Motion Winds - Tempo de Processamento: {round((time.time() - start), 2)} segundos.')
    
