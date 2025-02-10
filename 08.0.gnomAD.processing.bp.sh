#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4gb
#SBATCH --job-name=gnomAD.processing
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/08.0.gnomAD.processing\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

module load bcftools

# cd /xdisk/mliang1/qqiu/project/multiomics-hypertension/gnomAD
# Rscript 08.0.gnomAD.processing.R


cd /xdisk/mliang1/qqiu/project/multiomics-hypertension/data

# Randomly extract 30k SNPs
# shuf -n 300000 /xdisk/mliang1/qqiu/reference/dbSNP/dbSNP_list.txt > random_30k_snps.txt
# Remove BP SNPs
# grep -v -F -f /xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_blood_pressure.snp.txt random_30k_snps.txt > random_30k_snps_filtered.txt


# Define the directory containing your VCF files
vcf_dir="/xdisk/mliang1/qqiu/reference/gnomAD"
# Define the output file for the combined random SNPs
combined_output="gwas_catalog_bp_relevant.snp.vcf"
# List of chromosomes
chromosomes=$(echo {1..22} X Y)  # Modify if you have X, Y chromosomes
chromosomes=$(echo 12)

# Extract SNPs from each chromosome and combine them
for chr in $chromosomes; do
 vcf_file="${vcf_dir}/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz"
 output_file="chr${chr}_bp_snps.vcf"
  
  # Extract the random SNPs for each chromosome
 bcftools view --include ID==@/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_bp_relevant.snp.txt ${vcf_file} > ${output_file}
  
done


ls chr*vcf > vcf_list.txt 
bcftools concat -f vcf_list.txt -o gwas_catalog_bp_relevant.snp.vcf

# Remove temporary files
# rm chr*_random_snps.vcf



# awk 'NR==FNR { pos[$2] = $1; next } FNR==1 { print $0, "POS_dbSNP"; next } { print $0, pos[$1] }' /xdisk/mliang1/qqiu/reference/dbSNP/dbSNP.pos.chr.reformat.txt gwas_catalog_bp_relevant.snp.txt > gwas_catalog_bp_relevant.snp.pos.txt
