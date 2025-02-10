#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --job-name=cellbender_RLK3
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/cellbender_RLK3\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

### this script is to run cellbender with GPU and CUDA. 
## make sure in the folder, there is a file from the cell ranger count with name of raw_feature_bc_matrix.h5
## the parameter for the --expected-cells  and --total-droplets-included need to be estimated from the cell ranger count report plot
## --fpr 0.01 --epochs 150 

module load cuda11/11.2

[[ -d /xdisk/mliang1/qqiu/project/multiomics-hypertension/cellbender/RLK3 ]] || mkdir /xdisk/mliang1/qqiu/project/multiomics-hypertension/cellbender/RLK3
cd /xdisk/mliang1/qqiu/project/multiomics-hypertension/cellbender/RLK3

apptainer exec --nv /home/u1/qqiu/software/cellbender_latest.sif cellbender remove-background --input /xdisk/mliang1/qqiu/project/multiomics-hypertension/cellranger/RLK3/outs/raw_feature_bc_matrix.h5 --output RLK3\_cellbender_output.h5 --cuda --expected-cells 9041 --total-droplets-included 16000 --fpr 0.01 --epochs 150 --posterior-batch-size 5 --cells-posterior-reg-calc 20

