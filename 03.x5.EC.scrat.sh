#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=32gb
#SBATCH --constraint=hi_mem
#SBATCH --job-name=scrat
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/scrat\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

module load R/4.2.2

cd /xdisk/mliang1/qqiu/project/multiomics-hypertension/src
Rscript 03.x5.EC.scrat.R

