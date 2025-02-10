#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=gpu_standard
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5gb
#SBATCH --gres=gpu:1
#SBATCH --job-name=03.x1.EC.merge
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/03.x1.EC.merge\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

module load cuda11/11.8
source /opt/ohpc/pub/apps/anaconda/2022.05/etc/profile.d/conda.sh
conda activate cell2location
module load R/4.2.2

#which python
#python --version
#python -c "import scanpy as sc; print(sc.__version__)"

export LD_LIBRARY_PATH=/home/u1/qqiu/my_libs:$LD_LIBRARY_PATH
export PATH="/home/u1/qqiu/.conda/envs/cell2location/bin:$PATH"
#export RETICULATE_PYTHON="/home/u1/qqiu/.conda/envs/cell2location/bin/python3.9"

cd /xdisk/mliang1/qqiu/project/multiomics-hypertension/src
Rscript 03.x1.EC.merge.R

