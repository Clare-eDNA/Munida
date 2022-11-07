# Phyloseq looking at data for eDNA for the Munida Transect
## 2021 Feb 24

#setwd("~/Documents/MunidaTransect/data/FinalAnalyses/Munida-usearch")
setwd("/Volumes/CLARE_PHD/Documents/MunidaTransect/data/FinalAnalyses/Munida-usearch")

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("readxl")
library("dplyr")
library("vegan")
theme_set(theme_bw())

# read in OTU table
# need to remove the top bit of the otu_table where it has biom stuff
otu<-read.table(file="otu_table.txt",sep="\t", header=TRUE)
head(otu)

#read in taxonomy
tax<-read.csv(file="taxonomy2.csv", header=TRUE)
head(tax)

#read in metadata
metadata<-read.csv("Fixedmetadata2.csv")
head(metadata)

# define row names from the otu column
row.names(otu) <- otu$OTU_ID
# remove row names
otu <- otu %>% select (-OTU_ID)
# define row names for taxonomies
row.names(tax) <- tax$OTU_ID
tax <- tax %>% select (-OTU_ID)
# define row names for metadata (samples)
row.names(metadata) <- metadata$SampleID
metadata <- metadata %>% select (-SampleID)

#phy_tree <- read_tree("phyloseq-tree-21Dec20.nwk")
#taxa_names(phy_tree)
#taxa_names(OTU)

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
sample_names(Munida)
rank_names(Munida)
colnames(tax_table(Munida)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
sample_variables(Munida)

## Normalize number of reads in each sample using median sequenceing depth

#plot out the OTUs vs sample
#this code below takes a long time
#rarecurve(t(otu_table(Munida)), step=50, cex=0.5)

#total=median(sample_sums(Munida))
#standf=function(x,t=total) round(t * (x/sum(x)))
#Munida2=transform_sample_counts(Munida,standf)
# The number of reads used for normalization is 62,018

## Take off OTUs that have less than 10 sequences alltogether (1 taxa has 10 sequences)
Munida_1 <- prune_taxa(taxa_sums(Munida) > 9, Munida)
Munida_1 # samples taken off @ filtering point in bioinfo protocol

###
##
# rarefy without replacement
# grab 90%

#Munida.rarefied = rarefy_even_depth(Munida, rngseed=1, sample.size=0.9*min(sample_sums(Munida)), replace=F)
sum(otu_table(Munida))
## How many taxa per family (remember to subtract the empties)
get_taxa_unique(Munida, "Family") #66
#genus
get_taxa_unique(Munida, "Genus") #75 not counting humans, 76 counting humans
# etc...
get_taxa_unique(Munida, "Species") # ~about 77 depending if you count strains of bacteria or not
get_taxa_unique(Munida, "Phylum") #14

####### Separate out by day ########

#re-level the factors
sample_data(Munida)$Samp <- factor(sample_data(Munida)$Samp, levels = c("One","Two","Three","Four","Five","Six","Seven","Eight", "Nine"))
sample_data(Munida)$WaterType <- factor(sample_data(Munida)$WaterType, levels = c("Neritic","Mixed","SubTropical","Front","SubAntarctic","Depth"))

MunidaD1<-subset_samples(Munida, Day=="Feb-02")
MunidaD2<-subset_samples(Munida, Day=="Feb-23")

MunidaD1
MunidaD2

# go to "RestructureStats.R" script










####### Get an OTU table of only genus ########

Munida

# filter by genus only
Munida_genus_only <- tax_glom(Munida, "Genus")

#Get variables you want
variable1 = as.character(get_variable(Munida_genus_only, "JulianDate"))
variable2 = as.character(get_variable(Munida_genus_only, "WaterType"))

# paste
sample_data(Munida_genus_only)$NewPastedVar <- mapply(paste0,variable1, variable2, collapse = "_")

#merge by only the desired data columns
Munida_Merged=merge_samples(Munida_genus_only, "NewPastedVar")

# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(Munida_Merged), "matrix")
# transpose if necessary
if(taxa_are_rows(Munida)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

head(OTUdf)
OTUdf[OTUdf > 0] <- "Present"

write.csv(OTUdf,'Munida_Genus_Collapsed_OTU_table_present.csv')



####### And filter only by phyla

# filter by genus only
Munida_Phylum_only <- tax_glom(Munida, taxrank="Phylum")

#Get variables you want
variable1 = as.character(get_variable(Munida_genus_only, "JulianDate"))
variable2 = as.character(get_variable(Munida_genus_only, "WaterType"))

# paste
sample_data(Munida_Phylum_only)$NewPastedVar <- mapply(paste0,variable1, variable2, collapse = "_")

#merge by only the desired data columns
Munida_phylum_merged=merge_samples(Munida_Phylum_only, "NewPastedVar")
get_taxa_unique(Munida_phylum_merged, "Phylum")
# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(Munida_phylum_merged), "matrix")
# transpose if necessary
if(taxa_are_rows(Munida_phylum_merged)){OTU1 <- t(OTU1)}
get_taxa_unique(Munida_phylum_merged, "Phylum")
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

head(OTUdf)
OTUdf[OTUdf > 0] <- "Present"

write.csv(OTUdf,'Munida_Phylum_Collapsed_OTU_table_present.csv')
# will have to replace names
get_taxa_unique(Munida_phylum_merged, "Phylum")
