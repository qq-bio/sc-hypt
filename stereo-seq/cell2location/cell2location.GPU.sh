#!/bin/bash
#SBATCH --job-name=cell2location
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:1
#SBATCH --mem-per-cpu=12gb
#SBATCH --time=99:00:00
#SBATCH --account=mliang
#SBATCH --output=/home/qqiu/src/stereo-seq/cell2location/cell2location.GPU.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@mcw.edu

export PYTHONNOUSERSITE="literallyanyletters"
module load miniconda3
conda activate cell2loc_env

cd /home/qqiu/src/stereo-seq/cell2location/
python cell2location.GPU.test.py
