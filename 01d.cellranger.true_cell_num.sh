folder=/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellranger/

cd $folder
for i in `ls *LK*/outs/web_summary.html`
do
sample=$(echo $i | cut -d / -f 1)

if ! grep -q $sample /xdisk/mliang1/qqiu/project/multiomics-hypertension/cellbender/cellbender_parameter.txt; then
  num=$(cat $i | sed 's/[{}]/\n/g' | grep "100% Cells" | sed 's/.*(//' | sed 's/\/.*//')
  echo $sample $num >> /xdisk/mliang1/qqiu/project/multiomics-hypertension/cellbender/cellbender_parameter.txt
fi
done
