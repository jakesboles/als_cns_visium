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
import flax
import cell2location as c2l
import matplotlib.pyplot as plt
from matplotlib import rcParams
from cell2location.plt import plot_spatial
from cell2location import run_colocation
from cell2location.utils import select_slide
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from sklearn.model_selection import train_test_split
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs
os.chdir("/gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/mcx/")

adata_vis = sc.read_h5ad("/gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/mcx/adata_objects/adata_vis_after.h5ad")

adata_vis = adata_vis[adata_vis.obs['region'] == "WM"]

from cell2location import run_colocation
# res_dict, adata_vis = run_colocation(
#     adata_vis,
#     model_name='CoLocatedGroupsSklearnNMF',
#     train_args={
#       'n_fact': np.arange(5,20), # IMPORTANT: use a wider range of the number of factors (5-30)
#       'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
#       'n_restarts': 3 # number of training restarts
#     },
#     # the hyperparameters of NMF can be also adjusted:
#     model_kwargs={'alpha': 0.01, 'init': 'random', "nmf_kwd_args": {"tol": 0.000001}},
#     export_args={'path': 'NMF_WM/'}
# )

res_dict, adata_vis = run_colocation(
    adata_vis,
    model_name='ArchetypalAnalysis', 
    train_args={
        'n_fact': np.arange(5, 20),  
        'sample_id': 'sample_id',    # Note: ArchetypalAnalysis uses 'sample_id', not 'sample_name_col'
        'n_iter': 5000,              # Replaces n_restarts / sklearn training args
        'random_state': 0
    },
    export_args={'path': 'Archetypes_WM/'}
)

adata_vis.write_h5ad("AfterNMF_WM.h5ad")

with open("nmf_res_WM.pkl", "wb") as f:
    pickle.dump(res_dict, f)
