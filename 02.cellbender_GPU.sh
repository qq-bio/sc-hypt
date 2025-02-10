para_file=/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/cellbender_parameter.txt
log_folder=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log

while read -r line
do
sample_id=$(echo $line | awk '{print $1}')
ec=$(echo $line | awk '{print $2}')
tdi=$(echo $line | awk '{print $3}')
input_folder=/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellranger/$sample_id/outs
out_folder=/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellbender/$sample_id

echo "#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --job-name=cellbender_$sample_id
#SBATCH --output=$log_folder/cellbender_$sample_id\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

### this script is to run cellbender with GPU and CUDA. 
## make sure in the folder, there is a file from the cell ranger count with name of raw_feature_bc_matrix.h5
## the parameter for the --expected-cells  and --total-droplets-included need to be estimated from the cell ranger count report plot
## --fpr 0.01 --epochs 150 

module load cuda11/11.2

[[ -d $out_folder ]] || mkdir $out_folder
cd $out_folder

apptainer exec --nv /home/u1/qqiu/software/cellbender_latest.sif cellbender remove-background --input $input_folder/raw_feature_bc_matrix.h5 --output $sample_id\_cellbender_output.h5 --cuda --expected-cells $ec --total-droplets-included $tdi --fpr 0.01 --epochs 150 --posterior-batch-size 5 --cells-posterior-reg-calc 20
" > /xdisk/mliang1/qqiu/project/multiomics-hypertension/src/02.cellbender_GPU/$sample_id\.cellbender_GPU.sh
done < $para_file


for file in `ls /xdisk/mliang1/qqiu/project/multiomics-hypertension/src/02.cellbender_GPU/*.sh`
do
sample_id=$(basename $file | cut -d . -f 1)
if [ ! -f $log_folder/cellbender_$sample_id\-*.out ]
then
echo sbatch $file
fi
done

