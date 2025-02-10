#!/bin/bash
#SBATCH --account=mliang1
#SBATCH --partition=standard
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4gb
#SBATCH --job-name=eqtl_processing
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/08.eqtl_process\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu


cd /xdisk/mliang1/qqiu/data/GTEx/GTEx_Analysis_v8_eQTL

input=Artery_Tibial.v8.signif_variant_gene_pairs.txt
sample=$(echo $input | cut -d. -f1)

awk 'NR > 1 { print $1 }' $input | sort | uniq > $sample.snp_ids.txt
awk 'NR==FNR { snp_ids[$1]; next } $1 in snp_ids' $sample.snp_ids.txt /xdisk/mliang1/qqiu/reference/dbSNP/dbSNP.pos.chr.reformat.txt > $sample.filtered_dbSNP.txt


awk '
BEGIN {
    FS = OFS = "\t"
}
FILENAME == ARGV[1] {
    if ($1 ~ /^chr/) {
        rs_ids[$1] = $2
    }
    next
}
FILENAME == ARGV[2] {
    if ($1 ~ /^ENSG/) {
        gene_names[$1] = $7
    }
    next
}
FILENAME == ARGV[3] {
    if (FNR == 1) {
        print $0, "rs_id", "gene_name"
    } else {
        gene_id = $2
        sub(/\..*/, "", gene_id)  # Remove the dot and everything after it
        print $0, rs_ids[$1], gene_names[gene_id]
    }
}
' $sample.filtered_dbSNP.txt /xdisk/mliang1/qqiu/reference/biomaRt.gene.GRCH38.out $input > $sample.v8.signif_variant_gene_pairs.processed.txt