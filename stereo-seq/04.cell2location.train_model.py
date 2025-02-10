import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys
import cell2location
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42

os.chdir('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cell2location')

results_folder = '/xdisk/mliang1/qqiu/project/multiomics-hypertension/cell2location'
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

# load data & process
sn_data_file = sys.argv[1]
count_data_file = sys.argv[2]
output_file = sys.argv[3]


adata_sn = sc.read(sn_data_file)
count_sn = sc.read_10x_h5(count_data_file)

## transfer count data to adata
adata_sn = adata_sn.raw.to_adata()
adata_sn.layers["counts"] = count_sn.X.astype(int)
adata_sn.X = adata_sn.layers["counts"]


# filter gene
from cell2location.utils.filtering import filter_genes
## increase cut-off to exclude more genes, to select between 8k-16k genes
selected = filter_genes(adata_sn, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.05)
adata_sn = adata_sn[:, selected].copy()


# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_sn,
                        # 10X reaction / sample / batch
                        batch_key='orig.ident',
                        # cell type, covariate used for constructing signatures
                        labels_key='subclass_level2'
                       )

# create and train the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_sn)
# os.environ["SLURM_NTASKS_PER_NODE"]="10"
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)
#mod.plot_history(20)

## export the estimated cell abundance (summary of the posterior distribution).
adata_sn = mod.export_posterior(
    adata_sn, 
    sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

## save model and anndata object with results
mod.save(f"{ref_run_name}", overwrite=True)
adata_file = f"{ref_run_name}/{output_file}"
adata_sn.write(adata_file)

# examine QC plots
#mod.plot_QC()
