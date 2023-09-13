cd $4
# use flag to keep properly paired reads
# read paired (0x1)
# read mapped in proper pair (0x2)
samtools view -h -q 255 -f 3 $1 > $2_STARAligned.filtered.bam

# sort by coordinates
samtools sort -@ $3 -o $2_STARAligned.filtered.SortByCoord.bam $2_STARAligned.filtered.bam

# index
samtools index -@ $3 $2_STARAligned.filtered.SortByCoord.bam

# bam to bigwig
# forward strand
#bamCoverage -b $2_STARAligned.filtered.SortByCoord.bam -o $2_STARAligned.filtered.SortByCoord.forward.bw --filterRNAstrand forward
# reverse strand
#bamCoverage -b $2_STARAligned.filtered.SortByCoord.bam -o $2_STARAligned.filtered.SortByCoord.reverse.bw --filterRNAstrand reverse
#bamCoverage -b $2_STARAligned.filtered.SortByCoord.bam -o $2_STARAligned.filtered.SortByCoord.bw  --normalizeUsing RPKM