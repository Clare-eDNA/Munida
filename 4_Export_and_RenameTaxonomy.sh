#!/bin/bash -e
# This script exports all the different taxonomies
# Then, it renames them and collapses the taxonomies

cd /nesi/nobackup/uoo02328/clare/Munida/taxonomy-Oct/coi

module load QIIME2/2020.8

qiime tools export \
--input-path coi_NB_taxonomy.qza \
--output-path coi_taxonomy

qiime tools export \
--input-path coi_rep_seqs.qza \
--output-path coi_rep_seqs

qiime tools export \
--input-path coi-unrooted-tree.qza \
--output-path coi_unrooted-tree

/nesi/project/uoo02328/programs/munida_transect/scripts/collapse_tables_by_taxonomy.py -f OTU-coi-Swarm-labeled.txt -t ./coi_taxonomy/taxonomy.tsv -p coi-output1



cd /nesi/nobackup/uoo02328/clare/Munida/taxonomy-Oct/crust

module load QIIME2/2020.8

qiime tools export \
--input-path Crust_NB_taxonomy.qza \
--output-path crust_taxonomy

qiime tools export \
--input-path Crust_rep_seqs.qza \
--output-path crust_rep_seqs

qiime tools export \
--input-path Crust-unrooted-tree.qza \
--output-path crust_unrooted-tree

/nesi/project/uoo02328/programs/munida_transect/scripts/collapse_tables_by_taxonomy.py -f OTU-Crust-Swarm-labeled.txt -t ./crust_taxonomy/taxonomy.tsv -p crust



module load QIIME2/2020.8

qiime tools export \
--input-path fish_NB_taxonomy.qza \
--output-path fish_taxonomy

qiime tools export \
--input-path fish_rep_seqs.qza \
--output-path fish_rep_seqs

qiime tools export \
--input-path fish-unrooted-tree.qza \
--output-path fish_unrooted-tree

/nesi/project/uoo02328/programs/munida_transect/scripts/collapse_tables_by_taxonomy.py -f OTU-fish-Swarm-labeled.txt -t ./fish_taxonomy/taxonomy.tsv -p fish


#####

The script is ready. I tested it on a few files, though not your new ones. It should be fine, but let me know. 

/nesi/project/uoo02328/programs/munida_transect/scripts/collapse_tables_by_taxonomy.py 

You can use it from there. Run help to see options:

/nesi/project/uoo02328/programs/munida_transect/scripts/collapse_tables_by_taxonomy.py -h

Once you have exported your tables, you just have to enter those in the command, with an optional prefix for output collapsed tables:

/nesi/project/uoo02328/programs/munida_transect/scripts/collapse_tables_by_taxonomy.py -f frequency_table.tsv -t taxonomy_table.tsv -p output1

It is also on the Munida transect GitHub.

Cheers,

Hugh

Hugh Cross