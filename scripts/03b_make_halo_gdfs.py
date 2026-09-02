import os
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, Point, MultiPolygon
import numpy as np
from tifffile import imread, imwrite
import matplotlib.pyplot as plt
from shapely.geometry import box
import math 
import matplotlib.colors as mcolors
import gc
from shapely.ops import unary_union

output_dir = "/gpfs/projects/b1169/boles/als_cns_visium/data/03b_make_halo_gdfs/"

# Define path to halo output 
roi_dir = "/gpfs/projects/b1169/boles/als_cns_visium/data/halo_annotations/mcx_meninges"

# Define path to images 
img_dir = "/gpfs/projects/b1042/Gate_Lab/boles/als_motor_circuit_visium/images_for_alignment/"

def convert_to_polygon(geom):
  
  # Skip if not at least 3 points
  unique_coords = list(dict.fromkeys(geom.coords))
  if len(unique_coords) < 3: 
    return None

  # Convert to Polygon 
  poly = Polygon(geom)
    
  # Fix invalid polygons
  if not poly.is_valid:
    poly = poly.buffer(0)
  if not poly.is_valid:
    print("Polygon still not valid")
    return None
  else:
    return poly
  
sample_ids = [
    "JSB146-1", "JSB146-8", "AN67-3", "AN67-7", "AN68-1"
]

for sample in sample_ids: 
  
  # Define output folder 
  sample_dir = f"{output_dir}{sample}/"
    
  # Load ST spot GDF 
  spot_gdf = gpd.read_file(f"{sample_dir}gdf_spots/gdf.shp")
    
  roi_name = sample
        
  roi_file = f"{roi_dir}{roi_name}_Scan1.unmixed.geojson"
    
  # Load ROIs
  rois = gpd.read_file(roi_file)

  # Check initial data type 
  if not (rois.geometry.type == 'LineString').all(): 
    print("Not all LineString before transformation")
    continue

  # Remove CRS 
  rois.geometry.crs = None
    
  # Delete unused columns 
  del rois['object_type']
  del rois['classification']
  del rois['isLocked']
    
  # Convert ROIs to Polygons 
  rois['geometry'] = rois.geometry.apply(convert_to_polygon)

  # Remove empty geometries 
  rois = rois[~rois.geometry.is_empty]
  rois = rois[rois.geometry.notnull()]
        
  # Save ROI GDF
  os.makedirs(f"{sample_dir}gdf_halo", exist_ok = True)
  rois.to_file(f"{sample_dir}gdf_halo/gdf.shp")
  
  # Identify ST spots covered by gray matter ROIs 
  spot_gdf_temp = spot_gdf[['geometry']].copy()
  rois_temp = rois[['geometry']].copy()
  joined = gpd.sjoin(spot_gdf_temp, rois_temp, how="inner", predicate="within")
  spot_gdf['in_roi'] = False
  spot_gdf.loc[joined.index, 'in_roi'] = True

  # Save meta data
  os.makedirs(f"{sample_dir}results", exist_ok = True)
  spot_gdf[['barcode', 'in_roi']].to_csv(f"{sample_dir}results/in_roi_mng.csv")
  
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
  
  # Plot ROIs
  rois.plot(ax=ax, facecolor="#FFDBBB", clip_on = True)

  # Plot ST spots in gray matter
  spot_gdf[spot_gdf["in_roi"]].plot(ax=ax, facecolor="#003151")

  # Crop axes 
  ax.set_xlim(minx, maxx)
  ax.set_ylim(miny, maxy)
  ax.axis("off")

  # Save plot
  fig.savefig(f"{sample_dir}results/plt_in_roi_mng.pdf", format = "pdf", bbox_inches = "tight")
    
  # Free memory
  plt.close(fig)
  del fig, spot_gdf, img
  gc.collect()
