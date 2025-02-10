import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

import cell2location

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42

os.chdir('/scratch/g/mliang/multiomics-CKD/cell2location')

results_folder = '/scratch/g/mliang/multiomics-CKD/cell2location'

# create result folder
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

# load data
sp_data_file = "/scratch/g/mliang/multiomics-CKD/spaceranger/multiomics-CKD.space_merge.h5ad"
sn_data_file="/scratch/g/mliang/multiomics-CKD/cluster/multiomics-CKD.rna.cluster.RNA-0.9.anno.h5ad"

adata_sp = sc.read(sp_data_file)
adata_sn = sc.read(sn_data_file)

adata_sn = adata_sn.raw.to_adata()
adata_sn.layers["counts"] = adata_sn.X.copy()

from cell2location.utils.filtering import filter_genes

# increase cut-off to exclude more genes, to select between 8k-16k genes
selected = filter_genes(adata_sn, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_sn = adata_sn[:, selected].copy()

cell2location.models.RegressionModel.setup_anndata(adata=adata_sn,
                        # 10X reaction / sample / batch
                        batch_key='orig.ident',
                        # cell type, covariate used for constructing signatures
                        labels_key='new.cluster.ids.v2'
                       )

from cell2location.models import RegressionModel
mod = RegressionModel(adata_sn)

# view anndata_setup as a sanity check
# mod.view_anndata_setup()

mod.train(max_epochs=250, use_gpu=True)

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_sn = mod.export_posterior(
    adata_sn, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
# adata_file = f"{ref_run_name}/sc.h5ad"
# adata_sn.write(adata_file)
# adata_file

if 'means_per_cluster_mu_fg' in adata_sn.varm.keys():
    inf_aver = adata_sn.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_sn.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_sn.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_sn.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_sn.uns['mod']['factor_names']

intersect = np.intersect1d(adata_sp.var_names, inf_aver.index)
adata_sp = adata_sp[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_sp, batch_key="orig.ident")

mod = cell2location.models.Cell2location(
    adata_sp, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=8,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)

mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True,
         )


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_sp = mod.export_posterior(
    adata_sp, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

mod.save(f"{run_name}", overwrite=True)

adata_sp.obsm['means_cell_abundance_w_sf'].to_csv("multiomics-CKD.cell2location.means_cell_abundance_w_sf.csv")
adata_sp.obsm['stds_cell_abundance_w_sf'].to_csv("multiomics-CKD.cell2location.stds_cell_abundance_w_sf.csv")
adata_sp.obsm['q05_cell_abundance_w_sf'].to_csv("multiomics-CKD.cell2location.q05_cell_abundance_w_sf.csv")
adata_sp.obsm['q95_cell_abundance_w_sf'].to_csv("multiomics-CKD.cell2location.q95_cell_abundance_w_sf.csv")



