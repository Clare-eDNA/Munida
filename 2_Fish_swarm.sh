#!/bin/bash
#SBATCH -A uoo02328
#SBATCH -J swarm
#SBATCH --time 2:00:00
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
cd /nesi/nobackup/uoo02328/clare/Munida/LongFish/filtering/ACTUALFiltering || exit

echo "starting script"

module load USEARCH/11.0.667-i86linux32
module load VSEARCH/2.14.2-GCC-9.2.0
module load swarm/2.2.2-GCC-9.2.0

echo "Loaded modules"

#This bit of code filters for a max ee (error rate) of 0.1 

for samp in *.fastq
do
	name=$(basename "${samp}" .fastq)
	echo "${name}"
	usearch -fastq_filter "${name}".fastq -fastq_maxee 0.5 -fastqout "${name}".filt.fastq -fastqout_discarded "${name}".dis.fastq
done

#This relabels things

echo "done with filtering"

for samp in *.filt.fastq
do
	name=$(basename "${samp}" .filt.fastq)
	echo "${name}"
	usearch -fastx_relabel "${samp}" -prefix "${name}". -fastqout relabeled_"${samp}"
done


# we are assuming that samples are labeled T1_1_1.fastq; the name command pulls out the basename of the sample
# the echo calls it out so we know what samples have been processed
# the usearch fastx_relabel relabels are sample fastq headings inside of the file with the sample name; each sequence has its own number

#put all the relabled into one file
cat relab*fastq >> allsamples.fq

echo "Now booting up bbmap in order to change the fastq file into a fasta file"

module load BBMap/38.81-gimkl-2020a
reformat.sh in=allsamples.fq out=allsamples.fa

echo "Swarm wants things linearized so here's a bit of code from the website to linearize things"

#swarm calls for the linearization of samples
awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' allsamples.fa > allsamps_linearized.fasta

#swarm wants thing de-repped with vsearch but honestly I think the names of everything are getting in the way?

echo "Swarm wants things to be dereplicated so we are using vsearch to derepelicate stuff"

module load VSEARCH/2.14.2-GCC-9.2.0
vsearch \
    --derep_fulllength allsamps_linearized.fasta \
    --sizeout \
    --relabel_sha1 \
    --fasta_width 0 \
    --output amplicons_dereped.fasta

echo "OK now all samples are linearized and de-repped"

echo "OK now we are going to just cut everything off that doesn't appear at least 5x with usearch"

module load USEARCH/11.0.667-i86linux32

usearch -sortbysize amplicons_dereped.fasta -fastaout Ampli-Full-Derep-5.fasta -minsize 5

echo "and because usearch is a complete BUTT and puts everyting in lowercase, we are changing things to uppercase"

awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' Ampli-Full-Derep-5.fasta > Amplicons-Fish.fasta
#here's a full, uppercase file to SWARM
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' allsamples.fa > allsamples_Fish.fasta

echo "Now files allsamples_Fish.fasta (all sequences with the proper names) and \
Amplicons-Fish.fasta (derepped sequences) are ready to be imported into SWARM and usearch (to make an OTU table)"



#### Swarm

module load swarm/2.2.2-GCC-9.2.0
module load USEARCH/11.0.667-i86linux32

mkdir swarm
cp *Fish.fasta swarm/
cd swarm/

## run swarm on dereplicated seqs

# the /nesi/project is where Hugh installed it
# Requires uppercase stuff
# OTU representatives is the new fasta file that is denoised
# amplicons-COI.fasta is the un-denoised file that we want to denoise
# /dev/null is the junk stuff that is discarded
# takes about ~5m

swarm -f -z -t 4 -w allsample_OTU_reps.fasta Amplicons-Fish.fasta > /dev/null

## create otu frequency table using usearch

# requires uppercase stuff
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' allsample_OTU_reps.fasta > allsample_OTU_reps_upper.fasta


usearch -fastx_relabel Amplicons-Fish.fasta -prefix fishOTU. -fastaout fish_swarmOTUs_12oct20.fasta -keep_annots
# run a zotu output table
usearch -otutab allsamples_Fish.fasta -otus fish_swarmOTUs_12oct20.fasta -otutabout OTU-fish-Swarm-labeled.txt

# old code
#usearch -otutab allsamples_Fish.fasta -otus allsample_OTU_reps_upper.fasta -otutabout OTU-Fish-Swarm.txt
# this command should have output a zotu table
# Takes about 15m
