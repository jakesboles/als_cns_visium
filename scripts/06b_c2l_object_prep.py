import os
import warnings
warnings.filterwarnings('ignore')
import torch
import pickle
import sklearn
import scipy as sp
import numpy as np
import scanpy as sc
import pandas as pd
import squidpy as sq
import anndata as ad
import matplotlib as mpl
import cell2location as c2l
import matplotlib.pyplot as plt
from matplotlib import rcParams
from cell2location.plt import plot_spatial
from cell2location import run_colocation
from cell2location.utils import select_slide
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from sklearn.model_selection import train_test_split
rcParams['pdf.fonttype'] = 42
os.chdir("/gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/CreateObjects")

# MCX single cell
adata = sc.read_mtx("mcx_sc_counts.mtx").T 
cell_meta = pd.read_csv("mcx_sc_cell_metadata.csv", index_col=0)
gene_meta = pd.read_csv("mcx_sc_gene_metadata.csv")

# Assign metadata
adata.obs = cell_meta
adata.var = gene_meta
adata.var.index = gene_meta["gene"]

# Save object
adata.write_h5ad("mcx_ref.h5ad")


# SC single cell
adata = sc.read_mtx("sc_sc_counts.mtx").T 
cell_meta = pd.read_csv("sc_sc_cell_metadata.csv", index_col=0)
gene_meta = pd.read_csv("sc_sc_gene_metadata.csv")

# Assign metadata
adata.obs = cell_meta
adata.var = gene_meta
adata.var.index = gene_meta["gene"]

adata.write_h5ad("sc_ref.h5ad")


# MCX spatial
adata = sc.read_mtx("mcx_spatial_counts.mtx").T 
cell_meta = pd.read_csv("mcx_spatial_cell_metadata.csv", index_col=0)
gene_meta = pd.read_csv("mcx_spatial_gene_metadata.csv")

# Assign metadata
adata.obs = cell_meta
adata.var = gene_meta
adata.var.index = gene_meta["gene"]

# Fix batch class
adata.obs['batch'] = adata.obs['batch'].astype(str)

adata.write_h5ad("mcx_query.h5ad")


# SC spatial
adata = sc.read_mtx("sc_spatial_counts.mtx").T 
cell_meta = pd.read_csv("sc_spatial_cell_metadata.csv", index_col=0)
gene_meta = pd.read_csv("sc_spatial_gene_metadata.csv")

# Assign metadata
adata.obs = cell_meta
adata.var = gene_meta
adata.var.index = gene_meta["gene"]

# Fix batch class
adata.obs['batch'] = adata.obs['batch'].astype(str)

adata.write_h5ad("sc_query.h5ad")
