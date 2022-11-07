#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J Qiime2-taxonomy
#SBATCH --time 3:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH -n 1
#SBATCH --mem=48G
#SBATCH --partition=large
#SBATCH --mail-user=clare.adams@postgrad.otago.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --output QiimeTax.%j.out # CHANGE map1 part each run
#SBATCH --error QiimeTax.%j.err # CHANGE map1 part each run

datadir="/nesi/nobackup/uoo02328/clare/Munida/taxonomy-Oct/crust/"
cd ${datadir} || exit

#This command will cut off the numbers from the rep seques file if still there. Please go check.
#it cuts stuff, with the ";" as the delimiter keeping only the first field (before the deimilter)
#Best to run manually
#cut -d ';' -f1 allsample_OTU_reps.fasta > Fish-Reps.fasta
#we did this at the end of the swarm script but be sure that you actually did it


datadir="/nesi/nobackup/uoo02328/clare/Munida/taxonomy-Oct/crust/"  #change to whatever directory you're using
OTUs="Crust_zotus.fasta"    #change to whatever your OTUs are, likely "${assay}"_zotus.fasta
OTUtab="Crust_OTU_table.txt"            #change to whatever your OTU table is, likely "${assay}"_OTU_table.txt
name="Crust"                            #change to whatever you're analyzing: COI, Fish, or Crustaceans
metadata="Fixedmetadata.txt"            # Fixedmetadata.txt should be the right file
classifier="MIDORI_UNIQ_GB239_lrRNA_crustQiiExt_classifier.qza" ## download database from Midori; Classifier script can run for a long time and take a lot of memory

module load QIIME2/2020.8

# Import the rep seqs 

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ${datadir}/${OTUs} \
  --output-path ${name}_rep_seqs.qza


echo "Classifying sequences; I think this is matching our sequences to the sequences in the database"

# classify sequences
qiime feature-classifier classify-sklearn \
 --i-classifier ${datadir}/${classifier} \
 --i-reads ${datadir}/${name}_rep_seqs.qza \
 --p-n-jobs 16 \
 --o-classification ${name}_NB_taxonomy.qza

echo "Create a visiualization of our own sequences with the taxonomy available"

qiime metadata tabulate \
  --m-input-file ${name}_NB_taxonomy.qza \
  --o-visualization ${name}_NB_taxonomy.qzv


echo "Importing our OTU table using biom"

biom convert -i \
 ${datadir}/${OTUtab} \
 -o ${name}.biom --to-hdf5

qiime tools import \
  --input-path ${name}.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path ${name}_freq_table.qza


echo "Creating a bar plot graph of the taxa in each sample according to the frequency table"

#This step may fuck up if the names are not 100% correct
qiime taxa barplot \
  --i-table ${name}_freq_table.qza \
  --i-taxonomy ${name}_NB_taxonomy.qza \
  --m-metadata-file ${metadata} \
  --o-visualization ${name}-bar-plots.qzv

echo "Generating a phylogeny using your rep seqs"

# Crust_rep_seqs with the underscore is the sequences
# Crust-rep-seqs is the alignment.... The alignment is needed to do phylogeny

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${name}_rep_seqs.qza \
  --o-alignment ${name}-aligned-rep-seqs.qza \
  --o-masked-alignment ${name}-masked-aligned-rep-seqs.qza \
  --o-tree ${name}-unrooted-tree.qza \
  --o-rooted-tree ${name}-rooted-tree.qza

echo "and based on your alingment/phylogeny, get some basic diversity data to play with"

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ${name}-rooted-tree.qza \
  --i-table ${name}_freq_table.qza \
  --p-sampling-depth 10 \
  --m-metadata-file ${metadata} \
  --output-dir ${name}-diversity-results
