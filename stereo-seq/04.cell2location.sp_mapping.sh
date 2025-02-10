#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:4
#SBATCH --mem=128gb
#SBATCH --job-name=cell2location_gpu
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/cell2location_gpu\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

module load cuda11/11.2
module load anaconda
conda activate cell2location

#export 'PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:512'

cd /xdisk/mliang1/qqiu/project/multiomics-hypertension/src/stereo-seq
/home/u1/qqiu/.conda/envs/cell2location/bin/python 04.cell2location.sp_mapping.py /xdisk/mliang1/qqiu/project/multiomics-hypertension/stereo_seq/result/849BL_A4.bin50.h5ad mouse.LV.h5ad mouse.LV 849BL_A4.bin50.mouse.LV.mapping.h5ad

