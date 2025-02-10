#!/bin/bash
#SBATCH --partition=standard
#SBATCH --account=mliang1
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5gb
#SBATCH --job-name=A02095A3.stereoPipeline_6.1
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/A02095A3.stereoPipeline\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

# module load singularity
module load python
# module load hdf5/openmpi/gcc/64

# dataDir=/xdisk/mliang1/qqiu/data/multiomics-hypertension/stereo_seq_062623/A02095A3/00.Rawdata


apptainer exec /xdisk/mliang1/qqiu/software/SAW_6.1.sif mapping \
--outSAMattributes spatial --outSAMtype BAM SortedByCoordinate \
--genomeDir /xdisk/mliang1/qqiu/reference/SAW/rn7/STAR_SJ100/ \
--runThreadN 16 \
--outFileNamePrefix /xdisk/mliang1/qqiu/project/multiomics-hypertension/stereo_seq/result/952BR_B4/00.mapping/E100067636_L01_103_1. \
--sysShell /bin/bash \
--stParaFile /xdisk/mliang1/qqiu/project/multiomics-hypertension/stereo_seq/result/952BR_B4/00.mapping/E100067636_L01_103_1.bcPara \
--readNameSeparator \" \" --outBAMsortingBinsN 50
