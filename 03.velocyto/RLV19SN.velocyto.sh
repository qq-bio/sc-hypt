#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5gb
#SBATCH --job-name=velocyto_RLV19SN
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/velocyto_RLV19SN\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu


module load samtools
module load anaconda/2020
source /home/u1/qqiu/.bashrc && conda activate velocyto

cd /xdisk/mliang1/qqiu/project/multiomics-hypertension/velocyto/RLV19SN
velocyto run10x -m /xdisk/mliang1/qqiu/reference/rn7_rmsk.gtf /xdisk/mliang1/qqiu/project/multiomics-hypertension/cellranger/RLV19SN /xdisk/mliang1/qqiu/reference/cellranger/refdata-gex-rat7-2022-sb-r108/mRatBN7/genes/genes.gtf


