#!/bin/bash -e
# 5. This script imports the collapsed and renamed taxonomies for each thing
# Then, it merges the tables and taxonomies
#

cd /nesi/nobackup/uoo02328/clare/Munida/taxonomy-Oct/merged/
module load QIIME2/2020.8

## COI

biom convert \
-i coi-output1_collapsed_frequency.tsv \
 -o coi.OTU.biom --to-hdf5

qiime tools import \
  --input-path coi.OTU.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path coi_freq_table.qza

biom convert \
-i coi-output1_collapsed_taxonomy.tsv \
 -o coi.tax.biom --to-hdf5

qiime tools import \
  --input-path coi-output1_collapsed_taxonomy.tsv \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --output-path coi_feat_table.qza

## Crust

biom convert \
-i crust_collapsed_frequency.tsv \
 -o crust.OTU.biom --to-hdf5

qiime tools import \
  --input-path crust.OTU.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path crust_freq_table.qza

biom convert \
-i crust_collapsed_taxonomy.tsv \
 -o crust.tax.biom --to-hdf5

qiime tools import \
  --input-path crust_collapsed_taxonomy.tsv \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --output-path crust_feat_table.qza

## Fish

biom convert \
-i fish_collapsed_frequency.tsv \
 -o fish.OTU.biom --to-hdf5

qiime tools import \
  --input-path fish.OTU.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path fish_freq_table.qza

biom convert \
-i fish_collapsed_taxonomy.tsv \
 -o fish.tax.biom --to-hdf5

qiime tools import \
  --input-path fish_collapsed_taxonomy.tsv \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --output-path fish_feat_table.qza


## Merge tables 

qiime feature-table merge \
  --i-tables coi_freq_table.qza \
  --i-tables crust_freq_table.qza \
  --i-tables fish_freq_table.qza \
  --p-overlap-method sum \
  --o-merged-table full_freq_table.qza

qiime feature-table merge-taxa \
  --i-data coi_feat_table.qza \
  --i-data crust_feat_table.qza \
  --i-data fish_feat_table.qza \
  --o-merged-data full_feat_table.qza

## Visualize just to check because like ... that's a good thing to do

qiime taxa barplot \
	--i-table full_freq_table.qza \
	--i-taxonomy full_feat_table.qza \
	--m-metadata-file Fixedmetadata.txt \
	--o-visualization test

# In this instance it actually looks OK

# make a tree so that you can export it according to the QIIME2_to_Phyloseq pdf

#Statistically, this is likely wrong because the sequences are of different lengths and also
# the sequeneces are from different markers as well

#this combines the rep seqs (sequences)
qiime feature-table merge-seqs \
   --i-data coi_rep_seqs.qza \
   --i-data Crust_rep_seqs.qza \
   --i-data fish_rep_seqs.qza \
   --o-merged-data merged-rep-seqs.qza

#this creates 
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences merged-rep-seqs.qza \
  --o-alignment merged-aligned-rep-seqs.qza \
  --o-masked-alignment merged-masked-aligned-rep-seqs.qza \
  --o-tree merged-unrooted-tree.qza \
  --o-rooted-tree merged-rooted-tree.qza


#####
## And now we export for use with Phyloseq in R! 

qiime tools export \
--input-path full_feat_table.qza \
--output-path full_feat_table

qiime tools export \
--input-path full_freq_table.qza \
--output-path full_freq_table

biom convert --input-fp full_freq_table/feature-table.biom --output-fp full_freq_table/otu_table.txt --to-tsv

#open up otu_table and change #OTUID to OTUID without the #

qiime tools export \
--input-path merged-unrooted-tree.qza \
--output-path merged-unrooted-tree

qiime tools export \
--input-path merged-rooted-tree.qza \
--output-path merged-rooted-tree

qiime tools export \
--input-path merged-rep-seqs.qza \
--output-path full_req_seqs