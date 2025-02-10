#!/bin/bash
#SBATCH --job-name=952BR_E4.annh5ad2rds
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=48:00:00
#SBATCH --account=mliang
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@mcw.edu

/scratch/u/qqiu/project/multiomics/stereo_seq/result

Rscript /home/qqiu/src/stereo-seq/annh5ad2rds/R --infile 952BR_E4.bin50.h5ad --outfile 952BR_E4.bin50.rds
