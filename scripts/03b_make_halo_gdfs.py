import os
import sys
import re
import glob
import math
import gc
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
from shapely.ops import unary_union
from tifffile import imread
import matplotlib.pyplot as plt

# Print output is redirected to a SLURM log file rather than a terminal, so
# it's block-buffered by default and won't show up until the buffer fills or
# the job ends -- this forces line buffering so `tail -f` on the log reflects
# progress in real time.
sys.stdout.reconfigure(line_buffering=True)

os.chdir("/gpfs/projects/b1169/boles/als_cns_visium")

# 03a_make_spot_gdfs.py's output (one GDF of ST spots per sample) -- read
# from, never written to, here.
spot_gdf_dir = "/gpfs/projects/b1169/boles/als_cns_visium/data/03a_make_spot_gdfs/"

# Halo annotations, organized as one subfolder per annotated region/feature
# (e.g. mcx_meninges, sc_gm), each containing one geoJSON per annotated
# sample. Every category under this directory is processed in one run.
annotations_dir = "/gpfs/projects/b1169/boles/als_cns_visium/data/halo_annotations/"

data_dir = "/gpfs/projects/b1169/boles/als_cns_visium/data/03b_make_halo_gdfs/"
results_dir = "/gpfs/projects/b1169/boles/als_cns_visium/results/03b_make_halo_gdfs/"
os.makedirs(data_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

# Define path to images
img_dir = "/gpfs/projects/b1042/Gate_Lab/boles/als_motor_circuit_visium/images_for_alignment/"

# Samples excluded from the cohort in 01_obj_creation.R. Checked again here
# so a stray annotation for an excluded sample (as happened before with
# AN68-1) never silently gets processed.
exclude = ["137-1", "137-2", "AN16-1", "AN68-1", "AN68-2", "AN69-7", "AN69-8"]

# A spot is called "in" a region if more than this fraction of its area is
# covered by that region's annotated polygon(s), rather than requiring the
# whole spot to fall inside it. Spots straddling a region boundary are left
# unclassified for that region either way; any tissue not covered by a
# majority in any annotated category is meant to be picked up downstream as
# whichever region wasn't explicitly annotated (e.g. gray matter, once
# meninges/white matter are accounted for).
OVERLAP_THRESHOLD = 0.5

# Each combined per-category figure is downsampled by this factor before
# plotting, same rationale as 03a_make_spot_gdfs.py -- keeps memory and file
# size sane when a category's figure holds ~20+ full-res image crops at once.
DOWNSAMPLE = 5

# Saved as PNG rather than PDF: a PDF keeps every spot/ROI as a separate
# vector path, which balloons file size fast once a figure has dozens of
# panels each with thousands of spots. A rasterized PNG at this DPI stays
# small regardless of spot count.
PLOT_DPI = 150

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

categories = sorted(
  d for d in os.listdir(annotations_dir)
  if os.path.isdir(os.path.join(annotations_dir, d))
)

print(f"Found {len(categories)} annotation categories: {', '.join(categories)}")

# Each sample's spot GDF is read at most once no matter how many categories
# reference it (e.g. AN67-3 shows up in mcx_meninges, mcx_ptdp, and mcx_wm).
spot_gdf_cache = {}

def load_spot_gdf(sample):
  if sample not in spot_gdf_cache:
    path = f"{spot_gdf_dir}{sample}/gdf_spots/gdf.parquet"
    if not os.path.exists(path):
      print(f"No spot GDF found for {sample} at {path} -- has 03a_make_spot_gdfs.py been run for it?")
      spot_gdf_cache[sample] = None
    else:
      spot_gdf_cache[sample] = gpd.read_parquet(path)
  return spot_gdf_cache[sample]

# One combined results table per sample, one boolean column per annotated
# category, built up as categories are processed below.
sample_results = {}

for category in categories:

  category_dir = f"{annotations_dir}{category}/"
  geojson_files = sorted(glob.glob(f"{category_dir}*.geojson"))

  print(f"== {category} ({len(geojson_files)} samples) ==")

  # Collected here, then plotted together into one combined figure per
  # category once the loop below is done.
  plot_entries = []

  for j, geojson_file in enumerate(geojson_files):

    # Halo filenames aren't consistently "<sample>_Scan1..." -- some are
    # "_Scan2" -- so the sample id is recovered by splitting on the scan
    # number rather than assuming a fixed suffix.
    filename = os.path.basename(geojson_file)
    sample = re.split(r"_Scan[0-9]+", filename)[0]

    if sample in exclude:
      print(f"  [{j + 1}/{len(geojson_files)}] Skipping {sample} -- excluded from cohort")
      continue

    print(f"  [{j + 1}/{len(geojson_files)}] {sample}")

    spot_gdf = load_spot_gdf(sample)
    if spot_gdf is None:
      continue

    # Load ROIs
    rois = gpd.read_file(geojson_file)

    # Check initial data type
    if not (rois.geometry.type == 'LineString').all():
      print(f"  [{j + 1}/{len(geojson_files)}] {sample}: not all annotations are LineString before transformation, skipping")
      continue

    # Remove CRS
    rois = rois.set_crs(None, allow_override=True)

    # Halo's export carries bookkeeping columns alongside the geometry
    # (object type, lock state, a classification label), and the exact set
    # of columns isn't guaranteed to be the same across annotation types --
    # e.g. an AI-model-derived category like mcx_ptdp need not match a
    # manually-drawn region like mcx_meninges. Since every polygon in a
    # given category's file is unioned into one region below regardless of
    # any per-polygon label, none of those columns are actually needed, so
    # they're dropped unconditionally instead of by name. NOTE: this assumes
    # every polygon in one file belongs to the single class implied by its
    # category folder -- if any category's geoJSON actually mixes multiple
    # classifications that should be treated differently, filter on
    # `classification` before this line rather than dropping it.
    rois = rois[["geometry"]]

    # Convert ROIs to Polygons
    rois["geometry"] = rois.geometry.apply(convert_to_polygon)

    # Remove empty/invalid geometries
    rois = rois[rois.geometry.notnull() & ~rois.geometry.is_empty]

    if len(rois) == 0:
      print(f"  [{j + 1}/{len(geojson_files)}] {sample}: no valid polygons after conversion, skipping")
      continue

    # Save ROI GDF as GeoParquet (see 03a_make_spot_gdfs.py for why)
    sample_dir = f"{data_dir}{sample}/"
    os.makedirs(f"{sample_dir}gdf_halo", exist_ok=True)
    rois.to_parquet(f"{sample_dir}gdf_halo/{category}.parquet")

    # Identify ST spots with >OVERLAP_THRESHOLD of their area covered by
    # this category's (dissolved) annotated region.
    roi_union = unary_union(rois.geometry)
    spot_area = spot_gdf.geometry.area
    overlap_area = spot_gdf.geometry.intersection(roi_union).area
    in_roi = (overlap_area / spot_area) > OVERLAP_THRESHOLD

    if sample not in sample_results:
      sample_results[sample] = pd.DataFrame(index=spot_gdf["barcode"].values)
      sample_results[sample].index.name = "barcode"
    sample_results[sample][f"in_{category}"] = in_roi.values

    plot_entries.append((sample, spot_gdf, rois, in_roi))

  # Combined diagnostic figure for this category: one panel per sample
  # annotated for it, instead of one file per sample to click through.
  if len(plot_entries) == 0:
    print(f"  No samples produced valid ROIs for {category}, skipping its plot")
    continue

  print(f"  Saving combined plot for {category} ({len(plot_entries)} panels)")

  ncols = math.ceil(math.sqrt(len(plot_entries)))
  nrows = math.ceil(len(plot_entries) / ncols)
  fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows), squeeze=False)
  axes = axes.flatten()

  for ax, (sample, spot_gdf, rois, in_roi) in zip(axes, plot_entries):

    img = imread(f"{img_dir}{sample}.tiff")

    minx, miny, maxx, maxy = spot_gdf.geometry.total_bounds
    minx = max(0, int(minx) - 100)
    maxx = min(img.shape[1], int(maxx) + 100)
    miny = max(0, int(miny) - 100)
    maxy = min(img.shape[0], int(maxy) + 100)

    img_window = img[miny:maxy, minx:maxx]

    ax.imshow(img_window[::DOWNSAMPLE, ::DOWNSAMPLE], extent=[minx-0.5, maxx-0.5, maxy-0.5, miny-0.5])

    # Plot ROIs
    rois.plot(ax=ax, facecolor="#FFDBBB", alpha=0.4)

    # Plot ST spots in this category's region
    spot_gdf[in_roi].plot(ax=ax, facecolor="#003151", alpha=0.4)

    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)
    ax.axis("off")
    ax.set_title(sample, fontsize=8)

    del img, img_window
    gc.collect()

  for ax in axes[len(plot_entries):]:
    ax.axis("off")

  fig.suptitle(category)
  fig.savefig(f"{results_dir}{category}.png", format="png", dpi=PLOT_DPI, bbox_inches="tight")
  plt.close(fig)
  del fig
  gc.collect()

print(f"Saving combined results tables for {len(sample_results)} samples")

# Save one combined results table per sample (one boolean column per
# category that sample was annotated for -- a category a sample was never
# annotated for is simply absent as a column, not filled with a false/NaN).
for sample, df in sample_results.items():
  sample_dir = f"{data_dir}{sample}/"
  os.makedirs(f"{sample_dir}results", exist_ok=True)
  df.to_csv(f"{sample_dir}results/in_roi.csv")

print("Done")
