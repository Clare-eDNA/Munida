#!/usr/bin/env bash

#SBATCH -A uoo02328
#SBATCH -J fishBLAST
#SBATCH --time 1:30:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH -n 1
#SBATCH --mem=48G
#SBATCH --partition=large
#SBATCH --mail-user=clare.adams@postgrad.otago.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --output FishBLAST.%j.out # CHANGE map1 part each run
#SBATCH --error FishBLAST.%j.err # CHANGE map1 part each run


# a small BASH script
# by clare
# Note: Qiime2 was used to assign taxonomy, but BLAST can be useful for looking at the data.

module load BLAST/2.10.0-GCC-9.2.0
module load BLASTDB/2020-07-v5

cd /nesi/nobackup/uoo02328/clare/Munida/Fish || exit

blastn -query Fish-Reps.fasta -evalue 1e-10 -max_target_seqs 5 -qcov_hsp_perc 90 -perc_identity 50 -out BlastTry1.txt -outfmt 6 -num_threads 24 -db nt
