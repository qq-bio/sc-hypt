import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys
import cell2location
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42

os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:64"

os.chdir('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cell2location')

results_folder = '/xdisk/mliang1/qqiu/project/multiomics-hypertension/cell2location'
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

# load data & process
sp_data_file = sys.argv[1]
sn_data_file = sys.argv[2]
model_folder = sys.argv[3]
output_file = sys.argv[4]

adata_sp = sc.read(sp_data_file)
#adata_sn = sc.read(sn_data_file)

## remove mitochondria genes in visium data
adata_sp.var['MT_gene'] = [gene.startswith('Mt-') for gene in adata_sp.var.index]
adata_sp.obsm['MT'] = adata_sp[:, adata_sp.var['MT_gene'].values].X.toarray()
adata_sp = adata_sp[:, ~adata_sp.var['MT_gene'].values]


# load model and anndata
adata_file = f"{ref_run_name}/{sn_data_file}"
adata_sn = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}/{model_folder}", adata_sn)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_sn.varm.keys():
    inf_aver = adata_sn.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_sn.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_sn.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_sn.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_sn.uns['mod']['factor_names']
#inf_aver.iloc[0:5, 0:5]


# spatial mapping
## find shared genes
intersect = np.intersect1d(adata_sp.var_names, inf_aver.index)
adata_sp = adata_sp[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

## prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_sp, batch_key="orig.ident")
mod = cell2location.models.Cell2location(
    adata_sp, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=2,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection; test detection_alpha=200
    detection_alpha=200
)

## training cell2location
mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True
         )


# export the estimated cell abundance (summary of the posterior distribution).
adata_sp = mod.export_posterior(
    adata_sp, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

## save model & output
mod.save(f"{run_name}/{model_folder}", overwrite=True)
# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_sp)
adata_file = f"{run_name}/{output_file}"
adata_sp.__dict__['_raw'].__dict__['_var'] = adata_sp.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
# del(adata_sp.raw)
# del(adata_sp.var['_index'])
adata_sp.write(adata_file)

## load model & output
#adata_file = f"{run_name}/{output_file}"
#adata_sp = sc.read_h5ad(adata_file)
#mod = cell2location.models.Cell2location.load(f"{run_name}/{model_folder}", adata_sp)

adata_sp.obsm['means_cell_abundance_w_sf'].to_csv(f"{model_folder}.means_cell_abundance_w_sf.csv")
adata_sp.obsm['stds_cell_abundance_w_sf'].to_csv(f"{model_folder}.stds_cell_abundance_w_sf.csv")
adata_sp.obsm['q05_cell_abundance_w_sf'].to_csv(f"{model_folder}.q05_cell_abundance_w_sf.csv")
adata_sp.obsm['q95_cell_abundance_w_sf'].to_csv(f"{model_folder}.q95_cell_abundance_w_sf.csv")

#from cell2location import run_colocation
#res_dict, adata_sp = run_colocation(
#    adata_sp,
#    model_name='CoLocatedGroupsSklearnNMF',
#    train_args={
#      'n_fact': np.arange(11, 13), # IMPORTANT: use a wider range of the number of factors (5-30)
#      'sample_name_col': 'Sample_ID', # columns in adata_sp.obs that identifies sample
#      'n_restarts': 3 # number of training restarts
#    },
#    # the hyperparameters of NMF can be also adjusted:
#    model_kwargs={'alpha': 0.01, 'init': 'random', "nmf_kwd_args": {"tol": 0.000001}},
#    export_args={'path': f'{run_name}/CoLocatedComb/'}
#)


# Compute expected expression per cell type
#expected_dict = mod.module.model.compute_expected_per_cell_type(
#    mod.samples["post_sample_q05"], mod.adata_manager
#)

# Add to anndata layers
#for i, n in enumerate(mod.factor_names_):
#    adata_sp.layers[n] = expected_dict['mu'][i]

# Save anndata object with results
#adata_file = f"{run_name}/sp.h5ad"
#adata_sp.write(adata_file)
