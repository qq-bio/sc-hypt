#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --job-name=saw_build_index
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/saw_build_index\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

# module load singularity
SAW=/xdisk/mliang1/qqiu/software/SAW_6.1.sif


# cd /xdisk/mliang1/qqiu/reference/SAW/rn7/

singularity exec $SAW mapping --runMode genomeGenerate \
    --genomeDir /xdisk/mliang1/qqiu/reference/SAW/mm10/STAR_SJ100 \
    --genomeFastaFiles /xdisk/mliang1/qqiu/reference/SAW/mm10/genome/genome.fa \
    --sjdbGTFfile /xdisk/mliang1/qqiu/reference/SAW/mm10/genes/genes.gtf \
    --sjdbOverhang 99 \
    --runThreadN 8

