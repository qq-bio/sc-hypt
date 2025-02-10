#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --job-name=cell2location_gpu
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/cell2location_gpu\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

export SLURM_CLUSTER_NAME=Ocelote

module load cuda11/11.2
module load anaconda
conda activate cell2location

cd /xdisk/mliang1/qqiu/project/multiomics-hypertension/src/stereo-seq/
/home/u1/qqiu/.conda/envs/cell2location/bin/python 04.cell2location.GPU.py /xdisk/mliang1/qqiu/project/multiomics-hypertension/stereo_seq/result/952BR_B4-L.bin50.QC.h5ad /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.h5ad /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.count.10x.h5 952BR_B4-L

