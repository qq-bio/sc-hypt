#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5gb
#SBATCH --job-name=ss.cluster
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/ss.cluster\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

module load R/4.2.2

cd /xdisk/mliang1/qqiu/project/multiomics-hypertension/src
Rscript 02.cluster.wnn.harmony.ss.R

