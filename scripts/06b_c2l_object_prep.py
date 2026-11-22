# Assembles 06a_c2l_object_prep.R's per-tissue .mtx/.csv triplets into the
# .h5ad reference/query objects 07a/07b actually read. See 06a's own
# header for why the R-to-Python handoff goes through plain .mtx/.csv
# rather than zellkonverter/reticulate.
#
# Builds adata.obsm['spatial'] for each Visium query object from whichever
# pixel-coordinate columns actually survived into the metadata CSV --
# needed by cell2location's run_colocation() (08_nmf.py) for its
# per-sample spatial plots, which fall back to a coordinate-only scatter
# (cell2location.plt.plot_factor_spatial()) when adata.uns['spatial']
# (the full tissue image/scalefactors) isn't present, rather than the
# full histology-backed scanpy sc.pl.spatial(). Wiring up actual tissue
# images/scalefactors here would be a much bigger lift than this task
# asked for, and the core Cell2Location deconvolution math itself never
# needs spatial coordinates at all -- only this one downstream plotting
# path does. Checks a couple of candidate column name pairs rather than
# assuming one (this session can't inspect the real metadata.rds columns
# directly), and warns rather than crashing if neither is found.

import os
import sys
import scanpy as sc
import pandas as pd

sys.stdout.reconfigure(line_buffering=True)

os.chdir("/projects/b1169/boles/als_cns_visium")

data_dir = "data/06_c2l_object_prep/"

TISSUES = ["mcx", "sc"]

SPATIAL_COORD_CANDIDATES = [
    ("pxl_col_in_fullres", "pxl_row_in_fullres"),
    ("array_col", "array_row"),
]


def add_spatial_coords(adata):
    for col_x, col_y in SPATIAL_COORD_CANDIDATES:
        if col_x in adata.obs.columns and col_y in adata.obs.columns:
            adata.obsm["spatial"] = adata.obs[[col_x, col_y]].to_numpy(dtype=float)
            return adata
    print(f"WARNING: no recognized spatial coordinate columns found "
          f"(checked {SPATIAL_COORD_CANDIDATES}) -- adata.obsm['spatial'] "
          f"will not be set. Coordinate-based spatial plots downstream "
          f"(08_nmf.py's run_colocation() call) will be skipped or fail.")
    return adata


for tissue in TISSUES:

    print(f"INFO: Assembling {tissue} scRNAseq reference")

    adata = sc.read_mtx(f"{data_dir}{tissue}_sc_counts.mtx").T
    cell_meta = pd.read_csv(f"{data_dir}{tissue}_sc_cell_metadata.csv", index_col=0)
    gene_meta = pd.read_csv(f"{data_dir}{tissue}_sc_gene_metadata.csv")

    adata.obs = cell_meta
    adata.var = gene_meta
    adata.var.index = gene_meta["gene"]

    adata.write_h5ad(f"{data_dir}{tissue}_ref.h5ad")

    print(f"INFO: Assembling {tissue} Visium query")

    adata = sc.read_mtx(f"{data_dir}{tissue}_spatial_counts.mtx").T
    cell_meta = pd.read_csv(f"{data_dir}{tissue}_spatial_cell_metadata.csv", index_col=0)
    gene_meta = pd.read_csv(f"{data_dir}{tissue}_spatial_gene_metadata.csv")

    adata.obs = cell_meta
    adata.var = gene_meta
    adata.var.index = gene_meta["gene"]

    # Fix batch class -- written out as a factor level by R's write.csv()
    adata.obs['batch'] = adata.obs['batch'].astype(str)

    adata = add_spatial_coords(adata)

    adata.write_h5ad(f"{data_dir}{tissue}_query.h5ad")

print("INFO: Done")
