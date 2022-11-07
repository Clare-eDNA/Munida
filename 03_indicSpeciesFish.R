# Phyloseq looking at data for eDNA for the Munida Transect
# Script for R
## 2020 Oct 29

setwd("~/Documents/MunidaTransect/data/merged/phyloseq")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("readxl")
library("dplyr")
library("vegan")
theme_set(theme_bw())

# read in OTU table
otu<-read.table(file="otu_table.txt",sep="\t", header=TRUE)
head(otu)

#read in taxonomy
tax<-read.csv(file="taxonomy2.csv", header=TRUE)
head(tax)

#read in metadata
metadata<-read.csv("Fixedmetadata2.csv")
head(metadata)

# define row names from the otu column
row.names(otu) <- otu$OTU.ID
# remove row names
otu <- otu %>% select (-OTU.ID)
# define row names for taxonomies
row.names(tax) <- tax$Feature.ID
tax <- tax %>% select (-Feature.ID)
# define row names for metadata (samples)
row.names(metadata) <- metadata$SampleID
metadata <- metadata %>% select (-SampleID)

#Transform into matricies
otu_mat <- as.matrix(otu)
tax_mat <- as.matrix(tax)

#Tranform and read into phyloseq
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(metadata)

#Import into phyloseq
Munida <- phyloseq(OTU, TAX, samples)
Munida
### THIS "MUNIDA" is your phyloseq object! ###
sample_names(Munida)
rank_names(Munida)
sample_variables(Munida)


## How many taxa per family
get_taxa_unique(Munida, "Family")
get_taxa_unique(Munida, "Genus")
get_taxa_unique(Munida, "Species")
get_taxa_unique(Munida, "Phylum")

# remove negatives
Munida <- prune_samples(sample_names(Munida) != "Neg_1", Munida)
Munida <- prune_samples(sample_names(Munida) != "Neg_2", Munida)


#####################
## Filter for Species
#####################

## Species ##
#filter out what you don't want
phylaSFilter = c("", "s__")
Munida10 <- prune_taxa(taxa_sums(Munida) > 10, Munida)  ## everything has to have more than 10 reads
# and then subset the dataset
MunidaNR.Spe = subset_taxa(Munida10, !Species %in% phylaSFilter)
MunidaNR.Spe
tax_table(MunidaNR.Spe)
OTU_Tab_IndVal<-otu_table(MunidaNR.Spe)

#####################
## Filter for only chordates
#####################

# sequences that have not been rarified but have been assigned down to species-level taxonomy
# with no negatives
MunidaNR.Spe
sample_names(MunidaNR.Spe)


# Filter for vertebrata
# for crustaceans do :  Phylum=="p__Arthropoda_6656")
FishSpecies = subset_taxa(MunidaNR.Spe, Phylum=="p__Chordata_7711")

FishSpecies                     #the phyloseq object
FishSpe<-FishSpecies@otu_table  # get the otu table out of the phyloseq object
FishSpe2<-t(FishSpe)            # transpose (rotate) the otu table
FishSpe2[1:4,]                  # double check
metadata<-read.csv("Fixedmetadata2.csv")  # read in metadata
metadata2<-subset(metadata, JulianDate!=0)
metadata2[85:90,]               # double check

library(tibble)                 # change row to an actual column and name it
FishSpe3<-data.frame(FishSpe2)
FishSpe4<-rownames_to_column(FishSpe3, var="SampleID")

# Join based on sample ID
FishSpe5 <- left_join(metadata2, FishSpe4, by = c("SampleID" = "SampleID"))
FishSpe5<-subset(FishSpe5, JulianDate!=0) # get rid of negatives

# separate out days
FishD1<-subset(FishSpe5, JulianDate!=54)
FishD2<-subset(FishSpe5, JulianDate!=33)

# Load indicspecies package
library(indicspecies); packageVersion('indicspecies')

# set indicipieces factors
# this was before I used IndVal, so change r.g to IndVal!!

# all days
abund = FishSpe5[,20:44] #select just the ASV table part of the file
wat = FishSpe5$WaterType
FishSpeIn = multipatt(abund, wat, func = "r.g", control = how(nperm=9999))
summary(FishSpeIn)

#day 1
abundD1 = FishD1[,20:44] #select just the ASV table part of the file
watD1 = FishD1$WaterType
FishSpeInD1 = multipatt(abundD1, watD1, func = "r.g", control = how(nperm=9999))
summary(FishSpeInD1)

#day 2
abundD2 = FishD2[,20:44] #select just the ASV table part of the file
watD2 = FishD2$WaterType
FishSpeInD2 = multipatt(abundD2, watD2, func = "r.g", control = how(nperm=9999))
summary(FishSpeInD2)