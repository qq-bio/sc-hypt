#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5gb
#SBATCH --job-name=cellchat
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/cellchat\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

module load R/4.2.2

cd /xdisk/mliang1/qqiu/project/multiomics-hypertension/src
Rscript 03.x6.EC.cellchat.R

