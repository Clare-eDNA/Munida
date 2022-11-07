#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J swarm
#SBATCH --time 1:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH -n 1
#SBATCH --mem=48G
#SBATCH --partition=large
#SBATCH --mail-user=clare.adams@postgrad.otago.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --output swarm.%j.out # CHANGE 
#SBATCH --error swarm.%j.err # CHANGE 

# enter the directory
cd /nesi/nobackup/uoo02328/clare/Munida/LongCOI/filtering/Swarm || exit

echo "starting script"

module load USEARCH/11.0.667-i86linux32
module load VSEARCH/2.14.2-GCC-9.2.0
module load swarm/2.2.2-GCC-9.2.0

echo "Loaded modules"


## run swarm on dereplicated seqs

# the /nesi/project is where Hugh installed it
# Requires uppercase stuff
# OTU representatives is the new fasta file that is denoised
# amplicons-COI.fasta is the un-denoised file that we want to denoise
# /dev/null is the junk stuff that is discarded
# takes about ~5m

swarm -f -z -t 4 -w allsample_OTU_reps.fasta Amplicons_COI.fasta > /dev/null

## create otu frequency table using usearch

# requires uppercase stuff
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' allsample_OTU_reps.fasta > allsample_OTU_reps_upper.fasta


usearch -fastx_relabel Amplicons-COI.fasta -prefix coiOTU. -fastaout coi_swarmOTUs_12oct20.fasta -keep_annots
# run a zotu output table
usearch -otutab allsamples_COI.fasta -otus coi_swarmOTUs_12oct20.fasta -otutabout OTU-coi-Swarm-labeled.txt


#old code

#usearch -otutab allsamples_COI.fasta -otus allsample_OTU_reps_upper.fasta -otutabout OTU-Swarm.txt
# this command should have output a zotu table
# Takes about 15m