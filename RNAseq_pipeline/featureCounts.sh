cd $4

# sort by readnames
samtools sort -@ $3 -n -o $1_STARAligned.filtered.SortByRead.bam $1_STARAligned.filtered.bam

# index
samtools index -@ $3 $1_STARAligned.filtered.SortByRead.bam
# index
#  -B only count read pairs with both ends aligned
featureCounts -p -t exon -g gene_id -a $2 -o ./counts/${1}_featurecounts.txt ${1}_STARAligned.filtered.SortByRead.bam