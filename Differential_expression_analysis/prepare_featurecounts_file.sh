
for i in `cat sample_prefix.txt`
do
    cat  <(sed 1,2d ${i}_featurecounts.txt | awk '{print $NF}') > ./featureCounts/${i}_featureCounts.txt
done
sed 1,2d Ctr3_featurecounts.txt  | awk '{print $1}' > gene_names.txt
