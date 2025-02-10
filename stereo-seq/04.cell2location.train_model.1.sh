#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --job-name=cell2location_gpu
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/cell2location_gpu\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

module load cuda11/11.2
module load anaconda
conda activate cell2location

cd /xdisk/mliang1/qqiu/project/multiomics-hypertension/src/stereo-seq
/home/u1/qqiu/.conda/envs/cell2location/bin/python 04.cell2location.train_model.py /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.multiomics.anno.v2.h5ad /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.multiomics.anno.v2.count.10x.h5 mouse.LV.h5ad

