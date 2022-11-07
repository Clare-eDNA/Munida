#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J Fish
#SBATCH --time 2:00:00
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -n 1
#SBATCH --mem=48G
#SBATCH --partition=large
#SBATCH --mail-user=clare.adams@postgrad.otago.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --output COI.%j.out # CHANGE map1 part each run
#SBATCH --error COI.%j.err # CHANGE map1 part each run

#This script is for use on NeSI, and is a general pipeline
#This script takes barcode-assigned (demultiplexed) and paired COI sequences for the Munida Transect
#This script starts looking at them, size selecting them, and then filters them
#This script will also prepare for SWARM, although it is reccomended to use USEARCH/VSEARCH instead of SWARM
assay="Fish" ## change to COI, Fish, or Crustaceans (Crust)
datadir="/nesi/nobackup/uoo02328/clare/Munida/FINAL-GO/Fish/" #Change to whatever directory


cd "${datadir}"

module load Miniconda3
source activate /nesi/project/uoo02328/programs/miniconda_envs/OBItools

###################################
## Demultiplex based on Barcodes ##
###################################

cd /home/clare/MunidaTransect/RAW/

#you will want to run a fastqc of the raw files beforehand!

echo "Ensure you do a FASTQC"

# Use Obitools to demux using the barcode files
ngsfilter -t /home/clare/MunidaTransect/barcodes/Transect_FISH_FTP154_Barcode_Info.txt -u Unidentified_sequences.fasta -e 1 /home/clare/MunidaTransect/RAW/FTP154_R1 > /home/clare/MunidaTransect/25Jun/merged_assigned.fastq
ngsfilter -t /home/clare/MunidaTransect/barcodes/Transect_FISH_FTP154_Barcode_Info_1.txt -u Unidentified_sequences1.fasta -e 1 /home/clare/MunidaTransect/RAW/FTP154_R1 > /home/clare/MunidaTransect/25Jun/merged_assigned1.fastq
ngsfilter -t /home/clare/MunidaTransect/barcodes/Transect_FISH_FTP154_Barcode_Info_2.txt -u Unidentified_sequences2.fasta -e 1 /home/clare/MunidaTransect/RAW/FTP154_R1 > /home/clare/MunidaTransect/25Jun/merged_assigned2.fastq

#Put them all into one file
cat merged_assigned.fastq merged_assigned1.fastq merged_assigned2.fastq >> Fish_154.fastq

file="Fish_154.fastq"   #Change to whatever file you have

#Take a look at the head
obihead -n 1 "${file}"
#Take a look at the stats per sample
obistat -c sample "${file}" > stats.txt


##############################
## Obitools Size Control ##
##############################

# Trim by max length - there's a little bump at 340
# Filtered the bump at 340 out by filtering for 305 to 320 

# ALWAYS CHANGE THE LENGTH TO WHAT YOU NEED, this may vary depending on what 

obigrep -L 230 -l 190 "${file}" > Length.fastq

# large L for top length you'll accept and little l for smallest length you'll accept

echo "Length filter"

#if you want to look at the weird long lengths:
#obigrep -L 355 -l 320 COI_150.fastq > WeirdLongLength.fastq
#Weird stuff doesn't blast to anything

module load FastQC/0.11.9
fastqc Length.fastq

obiannotate -k sample Length.fastq > LengthAnnotated.fastq

obisplit -t sample LengthAnnotated.fastq

echo "spliting into groups"

mkdir filtering
cp T*.fastq filtering/
cp N*.fastq filtering/

cd filtering/

#############################
## USEARCH Qualtiy Control ##
#############################

module load USEARCH/11.0.667-i86linux32

echo "Starting filtering ee 1"

for samp in *.fastq
do
	name=$(basename "${samp}" .fastq)
	echo "${name}"
	usearch -fastq_filter "${name}".fastq -fastq_maxee 1 -fastq_maxns 0 -relabel "${name}". -fastqout "${name}".filt.fastq -fastaout "${name}".filt.fasta -fastqout_discarded "${name}".discard.fastq
done

echo "done filtering based on ee 1"

# put it all back into one file and change it to uppercase

cat *.filt.fasta >> pooled-"${assay}".fasta
tr ‘[:lower:]’ ‘[:upper:]’ < pooled-"${assay}".fasta > pooled_upper_"${assay}".fasta

###################################
## Dereplicate, considering size ##
###################################

#take all samples and make the minimum size of each OTU to be 10
module load VSEARCH/2.4.3-gimkl-2017a
vsearch --derep_fulllength pooled_upper_"${assay}".fasta -sizeout -relabel OTU. --output "${assay}"_derep.fasta

###################################
## Denoise the dataset w USEARCH ##
###################################

# denoise the dataset using unoise 3
#first sort the data
usearch -sortbysize "${assay}"_derep.fasta -fastaout uniques_sorted.fasta
#then get a fasta file with OTUs and zOTUs using the unoise3 alogrithm
#minimum OTU size is 10
usearch -unoise3 uniques_sorted.fasta -zotus "${assay}"_zotus.fasta -tabbedout unoise3_"${assay}"_report.txt -minsize 10
# gnerate an OTU table based on the thing
# For COI this takes like an hour
usearch -otutab pooled_upper_"${assay}".fasta -zotus "${assay}"_zotus.fasta -otutabout "${assay}"_OTU_table.txt

echo "done with filtering and re-naming"

## Your OTU table is "${assay}"_OTU_table.txt
## Your zOTUs (dereped sequences) are "${assay}"_zotus.fasta









