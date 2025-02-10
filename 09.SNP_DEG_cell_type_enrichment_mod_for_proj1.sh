#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4gb
#SBATCH --job-name=archR
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/09.SNP_DEG_cell_type_enrichment_mod_for_proj1\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

module load R/4.2.2


cd /xdisk/mliang1/qqiu/project/multiomics-hypertension/src
Rscript 09.SNP_DEG_cell_type_enrichment_mod_for_proj1.R

