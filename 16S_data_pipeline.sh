##################
# qiime2 #
##################
flash-merged single reads
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path [input_path] \
  --input-format 'CasavaOneEightSingleLanePerSampleDirFmt' \
  --output-path demux-paired-end-v34.qza 
wait

qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux-paired-end-v34.qza \
  --p-trim-left 0 \
  --p-trunc-len 0 \
  --o-table table-v34.qza \
  --o-representative-sequences rep-seqs-v34.qza \
  --p-n-threads 16 \
  --output-dir dada2-v34

wait

# checking statistics output
qiime tools export \
  --input-path dada2-v34/denoising_stats.qza \
  --output-path dada2-v34/
wait
# less dada2-v34/stats.tsv

qiime tools export --input-path rep-seqs-v34.qza --output-path exported-rep-seqs-v34 

# assigning taxonomy silva
qiime feature-classifier classify-sklearn \
  --i-classifier [qiime2_path]/qiime2/qiime2-2023.5/silva138/silva138_v34_classifier.qza \
  --i-reads rep-seqs-v34.qza \
  --o-classification taxonomy-v34.qza \
  --p-read-orientation 'same' \
  --p-n-jobs 32
wait

# assigning taxonomy homd
qiime feature-classifier classify-sklearn \
  --i-classifier [qiime2_path]/qiime2-2023.5/homd1523/homd_v34_classifier.qza \
  --i-reads rep-seqs-v34.qza \
  --o-classification homd-v34.qza \
  --p-read-orientation 'same' \
  --p-n-jobs 32
wait

qiime tools export --input-path homd-v34.qza --output-path exported-homd-v34 
less -S exported-homd-v34/taxonomy.tsv 

# Generate a phylogenetic tree for beta diversity analysis
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs-v34.qza --p-n-threads 32 --output-dir fasttree-tree 
wait
qiime tools export --input-path fasttree-tree/rooted_tree.qza --output-path exported-fasttree-tree  
wait

# make rep-seqs tree using pplacer algorithem
qiime fragment-insertion sepp --p-threads 32 --i-representative-sequences rep-seqs-v34.qza --i-reference-database [qiime2_path]/sepp/sepp-refs-silva-128.qza --output-dir insertion-tree 
wait
qiime tools export --input-path insertion-tree/tree.qza --output-path exported-insertion-tree  
wait



##################
# subset dataset #
##################
# use pplacer to assign taxonomy against 99 database
cd [input_path]
screen -S pplacer_flash1
conda activate metagenomic

mkdir pplacer_fasttree.100; cd pplacer_fasttree.100
pkg_fold=[pkg_fold] # ok

[pplacer_path]/hmmalign --dna -o asv.sto --mapali ${pkg_fold}/reference.sto ${pkg_fold}/reference.hmm ../exported-rep-seqs-v34/dna-sequences.fasta
wait
[pplacer_path]/pplacer --out-dir ./ -c ${pkg_fold} asv.sto
wait
rm classifications.db; [pplacer_path]/rppr prep_db -c ${pkg_fold} --sqlite classifications.db
wait
[pplacer_path]/guppy classify --sqlite classifications.db -c ${pkg_fold} asv.jplace
wait
[pplacer_path]/guppy tog -o asv.tre asv.jplace
wait
rm asv.csv; [pplacer_path]/sqlite3_script_v1.1.alpha16.sh -i classifications.db -o asv.csv
wait
[pplacer_path]/raw_results_table_RDPconvert.v4.R asv.csv asv.txt
wait
awk '{print $0}' ${pkg_fold}/seqinfo.csv | cut -d',' -f2,4 | awk '{gsub(",","\t",$0); print;}' >> asv.txt
wait
ulimit -s 65536 # enlarge stack limit to 64M, max value
# R --slave -e 'Cstack_info()["size"]'
[pplacer_path]/ncbi_blastn_otus2_tree_rename_exe.R asv.tre asv.txt
wait
[pplacer_path]/tree_plot.R asv.tre ${pkg_fold} asv.csv asv_tree
wait
# rm asv.sto



##########
# make a subset cohort including the surveyed samples only
##########
screen -S qiime2
cd [input_path]
mkdir subset; cd subset
conda activate qiime2-2023.5

# prepare a new taxonomy file based on pplacer, blastn, silva & homd alignment
qiime tools import --type 'FeatureData[Taxonomy]' --input-path ../exported-taxonomy-v34/taxonomy_pplacer.tsv --output-path taxonomy-pplacer.qza

# Filtering "non-Bacteria" or "mitochondria" or "chloroplast" ASV
qiime taxa filter-table --i-table ../table-v34.qza --i-taxonomy taxonomy-pplacer.qza --p-include "D_0__Bacteria;D_1__" --p-exclude "Mitochondria","Chloroplast" --o-filtered-table table-v34.qza
# qiime taxa filter-seqs --i-sequences ../rep-seqs-v34.qza --i-taxonomy taxonomy-pplacer.qza --p-include "D_0__Bacteria;D_1__" --p-exclude "Mitochondria","Chloroplast" --o-filtered-sequences rep-seqs-v34.qza

# prepare sample_subset.txt with header of sampleid, and remove rep sequences with <10 reads across the samples
# qiime feature-table filter-samples --i-table table-v34.qza --m-metadata-file ../sample_subset.txt --o-filtered-table template.qza; mv template.qza table-v34.qza
qiime feature-table filter-features --i-table table-v34.qza --p-min-frequency 10 --o-filtered-table template.qza; mv template.qza table-v34.qza

# filter rep sequences based on filtered features in the asv table
qiime feature-table filter-seqs --i-data ../rep-seqs-v34.qza --i-table table-v34.qza  --o-filtered-data rep-seqs-v34.qza

# exporting OTU feature table
qiime tools export --input-path table-v34.qza --output-path exported-otu-table-v34 
biom convert -i exported-otu-table-v34/feature-table.biom -o exported-otu-table-v34/feature-table-v34.tsv --to-tsv --table-type "OTU table" 
less -S exported-otu-table-v34/feature-table-v34.tsv

# exporting representative OTU sequences
qiime tools export --input-path rep-seqs-v34.qza --output-path exported-rep-seqs-v34 
less -S exported-rep-seqs-v34/dna-sequences.fasta 
# wait
qiime tools export --input-path taxonomy-pplacer.qza --output-path exported-taxonomy-pplacer 
less -S exported-taxonomy-pplacer/taxonomy.tsv 

##########
# collapsing feature tables by the taxonomy at the specific level (	--p-level 6==genus, 5==family, 4==order, 3==class, 2==phylum)
##########
nnn=2 
wait
qiime taxa collapse --i-table table-v34.qza --i-taxonomy taxonomy-pplacer.qza --p-level ${nnn} --o-collapsed-table collapsed-table
wait
qiime tools export --input-path collapsed-table.qza --output-path exported-collapsed-pplacer-v34
wait
biom convert -i exported-collapsed-pplacer-v34/feature-table.biom -o exported-collapsed-pplacer-v34/otu-table-v34-l${nnn}.tsv --to-tsv --table-type "OTU table" 
wait
mv exported-collapsed-pplacer-v34/feature-table.biom exported-collapsed-pplacer-v34/feature-table-v34-l${nnn}.biom
wait
qiime feature-table relative-frequency --i-table collapsed-table.qza --o-relative-frequency-table exported-collapsed-pplacer-v34/otu-table-v34-l${nnn}-frequency
wait
qiime tools export --input-path exported-collapsed-pplacer-v34/otu-table-v34-l${nnn}-frequency.qza --output-path exported-collapsed-pplacer-v34
wait
biom convert -i exported-collapsed-pplacer-v34/feature-table.biom -o exported-collapsed-pplacer-v34/otu-table-v34-l${nnn}-frequency.tsv --to-tsv --table-type "OTU table"; rm exported-collapsed-pplacer-v34/feature-table.biom; rm exported-collapsed-pplacer-v34/otu-table-v34-l${nnn}-frequency.qza
wait
rm collapsed-table.qza

