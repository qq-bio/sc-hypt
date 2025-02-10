para_file=/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/velocyto_parameter.txt

while read -r line
do
sample_id=$(echo $line | awk '{print $1}')
ref=$(echo $line | awk '{print $2}')
data_folder=/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellranger/$sample_id
out_folder=/xdisk/mliang1/qqiu/project/multiomics-hypertension/velocyto/$sample_id
log_folder=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log

if [ $ref = "hg38" ]
then
gtf=/xdisk/mliang1/qqiu/reference/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf
msk=/xdisk/mliang1/qqiu/reference/hg38_rmsk.gtf
elif [ $ref = "mm10" ]
then
gtf=/xdisk/mliang1/qqiu/reference/cellranger/refdata-gex-mm10-2020-A/genes/genes.gtf
msk=/xdisk/mliang1/qqiu/reference/mm10_rmsk.gtf
elif [ $ref = "rn7" ]
then
gtf=/xdisk/mliang1/qqiu/reference/cellranger/refdata-gex-rat7-2022-sb-r108/mRatBN7/genes/genes.gtf
msk=/xdisk/mliang1/qqiu/reference/rn7_rmsk.gtf
else
echo "Argument 2: Reference must be one of hg38, mm10 or rn7" 1>&2; exit 1
fi

[[ -d $out_folder ]] || mkdir $out_folder
# [[ -d $data_folder ]] || (mkdir $data_folder && cp -R $sample_folder/*fastq.gz $data_folder)

echo "#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5gb
#SBATCH --job-name=velocyto_$sample_id
#SBATCH --output=$log_folder/velocyto_$sample_id\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu


module load samtools
module load anaconda/2020
source /home/u1/qqiu/.bashrc && conda activate velocyto

cd $out_folder
velocyto run10x -m $msk $data_folder $gtf

" > /xdisk/mliang1/qqiu/project/multiomics-hypertension/src/03.velocyto/$sample_id\.velocyto.sh
done < $para_file
