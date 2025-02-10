#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4gb
#SBATCH --job-name=GTEx.vcf.processing
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/08.0.GTEx.vcf.processing\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

module load bcftools

bcftools view -i 'ID=@/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/eqtl.variant_id.txt' /xdisk/mliang1/qqiu/data/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz -Oz -o /xdisk/mliang1/qqiu/project/multiomics-hypertension/data/GTEx.eqtl.filtered.vcf.gz

