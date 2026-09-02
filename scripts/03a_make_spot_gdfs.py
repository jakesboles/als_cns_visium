import os
import sys
import math
import gc
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from tifffile import imread
import matplotlib.pyplot as plt

# Print output is redirected to a SLURM log file rather than a terminal, so
# it's block-buffered by default and won't show up until the buffer fills or
# the job ends -- this forces line buffering so `tail -f` on the log reflects
# progress in real time.
sys.stdout.reconfigure(line_buffering=True)

os.chdir("/gpfs/projects/b1169/boles/als_cns_visium")

data_dir = "/gpfs/projects/b1169/boles/als_cns_visium/data/03a_make_spot_gdfs/"
results_dir = "/gpfs/projects/b1169/boles/als_cns_visium/results/03a_make_spot_gdfs/"
os.makedirs(data_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

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

# Define ST spot radius in pixels
radius = 55

# The combined QC figure below stacks one panel per sample into a single
# file; each panel's crop is downsampled by this factor before plotting so
# ~47 full-res image crops don't all have to be held/rendered at full
# resolution at once. Fine for a by-eye overview, not meant for pixel-level
# inspection -- the full-res crop is still recomputed fresh here each time
# from the source tiff if that's ever needed.
DOWNSAMPLE = 5

# Saved as PNG rather than PDF: a PDF keeps every spot as a separate vector
# path, which balloons file size fast once a figure has dozens of panels
# each with thousands of spots. A rasterized PNG at this DPI stays small
# regardless of spot count.
PLOT_DPI = 150

print(f"Processing {len(sample_ids)} samples")

ncols = math.ceil(math.sqrt(len(sample_ids)))
nrows = math.ceil(len(sample_ids) / ncols)
fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows))
axes = axes.flatten()

for i, sample in enumerate(sample_ids):

  print(f"[{i + 1}/{len(sample_ids)}] {sample}")

  # Define output folder
  sample_dir = f"{data_dir}{sample}/"
  os.makedirs(sample_dir, exist_ok = True)

  # Load Space Ranger coordinates, keeping only spots actually called as
  # over tissue -- carrying the rest along is dead weight in every step
  # downstream (the GDF, the spatial join in 03b, the plots).
  coords = pd.read_csv(f"{spaceranger_dir}{sample}/outs/spatial/tissue_positions.csv")
  coords = coords[coords["in_tissue"] == 1]

  # Extract x,y coordinates. `barcode` is kept as a plain column rather than
  # the index, so it survives a round-trip through any file format
  # unambiguously instead of depending on that format's index handling.
  spot_gdf = gpd.GeoDataFrame({"geometry": [Point(x,y) for x,y in zip(coords['pxl_col_in_fullres'], coords['pxl_row_in_fullres'])],
                               "barcode": [sample + "_" + bc for bc in coords['barcode']]})

  # Construct ST spots
  spot_gdf['geometry'] = spot_gdf.geometry.buffer(radius)

  # Save ST GDF as GeoParquet -- a single file, no attribute field name
  # truncation, and no ambiguity about whether the index survives the
  # round-trip (it doesn't matter here since barcode is a real column).
  os.makedirs(f"{sample_dir}gdf_spots", exist_ok = True)
  spot_gdf.to_parquet(f"{sample_dir}gdf_spots/gdf.parquet")

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

  ax = axes[i]

  # Boundaries of min pixels: miny-0.5 to miny+0.5; minx-0.5 to minx+0.5
  # Boundaries of max pixels: maxy-1.5 to maxy-0.5; maxx-1.5 to maxx-0.5
  ax.imshow(img_window[::DOWNSAMPLE, ::DOWNSAMPLE], extent=[minx-0.5, maxx-0.5, maxy-0.5, miny-0.5]) # left, right, bottom, top

  # Plot ST spots
  spot_gdf.plot(ax=ax, facecolor = "#8DD3C7", alpha = 0.4)

  # Crop axes
  ax.set_xlim(minx, maxx)
  ax.set_ylim(miny, maxy)
  ax.axis("off")
  ax.set_title(sample, fontsize = 8)

  # Free memory
  del img, img_window, spot_gdf
  gc.collect()

# Turn off any unused panels (the grid is sized to the nearest full rows/cols)
for ax in axes[len(sample_ids):]:
  ax.axis("off")

print("Saving combined QC figure")
fig.suptitle("ST spots")
fig.savefig(f"{results_dir}ST_spots.png", format = "png", dpi = PLOT_DPI, bbox_inches = "tight")
plt.close(fig)

print("Done")
