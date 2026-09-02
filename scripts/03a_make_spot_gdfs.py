import os
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from tifffile import imread
import matplotlib.pyplot as plt
import gc
import json
os.chdir("/gpfs/projects/b1169/boles/als_cns_visium")

output_dir = "/gpfs/projects/b1169/boles/als_cns_visium/data/03a_make_spot_gdfs/"

# Define path to Space Ranger output 
spaceranger_dir = "/gpfs/projects/b1042/Gate_Lab/boles/als_motor_circuit_visium/spaceranger/"

# Define path to images 
img_dir = "/gpfs/projects/b1042/Gate_Lab/boles/als_motor_circuit_visium/images_for_alignment/"

sample_ids = [
"AN67-1",   "AN67-2",   "AN67-3",   "AN67-4",   "AN67-5",   "AN67-6",   "AN67-7",   "AN67-8",   "AN72-1",   "AN72-2",   "AN72-3",   "AN72-4",   "AN72-5",  
"AN72-6",   "AN72-7",   "AN72-8",   "JSB146-1", "JSB146-2", "JSB146-3", "JSB146-4", "JSB146-5", "JSB146-6", "JSB146-7", "JSB146-8", "JSB147-1", "JSB147-2",
"JSB147-3", "JSB147-4", "JSB147-5", "JSB147-6", "JSB164-1", "JSB164-2", "JSB164-3", "JSB164-4", "JSB164-5", "JSB164-6", "JSB164-7", "JSB164-8", "JSB171-1",
"JSB171-2", "JSB171-3", "JSB171-4", "JSB171-5", "JSB171-6", "JSB171-7", "JSB171-8", "JSB171-9"
]

pxl_per_µm = 0.5

# Define ST spot radius in pixels 
radius = 55

for sample in sample_ids: 
  
  # Define output folder 
  sample_dir = f"{output_dir}{sample}/"
  os.makedirs(sample_dir, exist_ok = True)

  # Load Space Ranger coordinates
  coords = pd.read_csv(f"{spaceranger_dir}{sample}/outs/spatial/tissue_positions.csv")

  # Extract x,y coordinates 
  spot_gdf = gpd.GeoDataFrame({"geometry": [Point(x,y) for x,y in zip(coords['pxl_col_in_fullres'], coords['pxl_row_in_fullres'])], 
                               "barcode": [sample + "_" + bc for bc in coords['barcode']]})

  # Specify index column
  spot_gdf = spot_gdf.set_index("barcode")

  # Construct ST spots
  spot_gdf['geometry'] = spot_gdf.geometry.buffer(radius)
    
  # Save ST GDF
  os.makedirs(f"{sample_dir}gdf_spots", exist_ok = True)
  spot_gdf.to_file(f"{sample_dir}gdf_spots/gdf.shp")

  # Load image    
  img = imread(f"{img_dir}{sample}.tiff")
    
  # Calculate bounds of Space Ranger window
  minx, miny, maxx, maxy = spot_gdf.geometry.total_bounds

  # Define x, y boundaries ~100 pixels beyond bounds of GDF
  minx = max(0, int(minx)-100)
  maxx = min(img.shape[1], int(maxx)+100)   
  miny = max(0, int(miny)-100)
  maxy = min(img.shape[0], int(maxy)+100) 

  # Crop image (x pixels: minx to maxx-1; y pixels: miny to maxy-1 pixel)
  img_window = img[miny:maxy, minx:maxx]

  # Initialize plot
  fig, ax = plt.subplots()

  # Boundaries of min pixels: miny-0.5 to miny+0.5; minx-0.5 to minx+0.5
  # Boundaries of max pixels: maxy-1.5 to maxy-0.5; maxx-1.5 to maxx-0.5
  ax.imshow(img_window, extent=[minx-0.5, maxx-0.5, maxy-0.5, miny-0.5]) # left, right, bottom, top

  # Plot ST spots
  spot_gdf.plot(ax=ax, facecolor = "#8DD3C7")

  # Crop axes 
  ax.set_xlim(minx, maxx)
  ax.set_ylim(miny, maxy)
  ax.axis("off")

  # Save plot
  fig.savefig(f"{sample_dir}ST_spots.pdf", format = "pdf", bbox_inches = "tight")
    
  # Free memory
  plt.close(fig)
  del fig, spot_gdf, img
  gc.collect()
  
