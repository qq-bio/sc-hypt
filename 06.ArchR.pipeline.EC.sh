#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4gb
#SBATCH --job-name=archR
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/archR-rat\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

module load R/4.2.2

export HDF5_USE_FILE_LOCKING=FALSE
export RHDF5_USE_FILE_LOCKING=FALSE

cd /xdisk/mliang1/qqiu/project/multiomics-hypertension/src
Rscript 06.ArchR.pipeline.EC.R

