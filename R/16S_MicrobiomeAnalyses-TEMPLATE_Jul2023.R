###### 16S Microbiome Analyses - Molecular Microbiology & Genomics Lab ######
### This R-Script includes routine microbiome analyses conducted using
### the 16S rRNA gene (V4 region (primers 515F/806R)).
### Initial data cleaning starts from QIIME2, so files needed are from QIIME2
### resulting files.
### Feel free to change anything in this script to fit your project analyses!
### Please make sure to keep this original copy and to save any edited versions
### with your project data!
### Created by Paisley Samuel - Summer 2023

###### List of Packages Used ######
## Here is a list of all the packages used throughout this R script.
## Take the necessary steps to install any packages that you do not have
## installed already on your computer (may have to use install.packages()
## or Bioconductor depending on the package)

library(ggplot2)
library(microbiome)
library(MMUPHin)
library(multcompView)
library(pgirmess)
library(phyloseq)
library(psych)
library(tidyverse)
library(vegan)

###### Set Working Directory and Seed - ALWAYS RUN THIS FIRST! ######
setwd_seed <- function(file_path,seed_number) {
  setwd(file_path) #include the file path to the folder where all your files are stored
  set.seed(seed_number) #choose any random number
}
## insert your files into the function
# "file_path" = insert the file path of your desired folder
# seed_number = insert a random number
setwd_seed("file_path", seed_number)

###### Producing Relative Abundance Data ######
### NOTE: If you ignore metadata, you will NOT be able to perform statistical
###       analyses on your data, only taxonomic analyses.

## RUN FUNCTIONS FIRST!
# use "calculate_abund.metadta" if you have your metadata file
# use "calculate_abund" if ignoring metadata (NO metadata file)
calculate_abund.metadata <- function(feature_csv,metadata_csv){
  library(vegan)
  dat<-t(data.matrix(read.csv(feature_csv, header=TRUE, row.names = 1)))
  #transposed so that the row names are now the sample names and the ASVs are the columns
  #BUT it is still a MATRIX
  metadata <- read.csv(metadata_csv, header = TRUE, row.names = 1)
  dat <- as.data.frame(dat)
  typeof(dat) #list
  common.rownames <- intersect(rownames(dat), rownames(metadata)) #returns the same row names found in each dataframe
  dat <- dat[common.rownames,] #subsets dat to only include the same row names as metadata
  metadata <- metadata[common.rownames,] #subsets metadata to only include the same row names as dat
  otu.abund<-which(colSums(dat)>2) #removes singletons and doubletons
  dat.dom<-dat[,otu.abund] #include dominant taxa
  dat.pa<-decostand(dat.dom, method ="pa") #turns dat.dom into presence/absence data (1/0)
  dat.otus.01per<-which(colSums(dat.pa) > (0.01*nrow(dat.pa)))
  dat.01per<-dat.dom[,dat.otus.01per] #removed ASVs that occur less than 0.1%
  dat.otus.001per<-which(colSums(dat.pa) > (0.001*nrow(dat.pa)))
  dat.001per<-dat.dom[,dat.otus.001per] #removed ASVs that occur less than 0.01%; increases the number of ASVs - includes more "microdiversity"
  dat.ra<-decostand(dat.01per, method = "total") #relative abundance of >1% taxa
  dfs_to_return <- list(as.data.frame(dat),as.data.frame(metadata),
                        as.data.frame(dat.dom),as.data.frame(dat.pa),
                        as.data.frame(dat.01per),as.data.frame(dat.001per),
                        as.data.frame(dat.ra))
  names(dfs_to_return) <- c("dat", "metadata","dat.dom","dat.pa","dat.01per","dat.001per","dat.ra")
  return(dfs_to_return)
}
calculate_abund <- function(feature_csv){
  library(vegan)
  dat<-t(data.matrix(read.csv(feature_csv, header=TRUE, row.names = 1)))
  #transposed so that the row names are now the sample names and the ASVs are the columns
  dat <- as.data.frame(dat)
  otu.abund<-which(colSums(dat)>2) #removes singletons and doubletons
  dat.dom<-dat[,otu.abund] #include dominant taxa
  dat.pa<-decostand(dat.dom, method ="pa") #turns dat.dom into presence/absence data (1/0)
  dat.otus.01per<-which(colSums(dat.pa) > (0.01*nrow(dat.pa)))
  dat.01per<-dat.dom[,dat.otus.01per] #removed ASVs that occur less than 0.1%
  dat.otus.001per<-which(colSums(dat.pa) > (0.001*nrow(dat.pa)))
  dat.001per<-dat.dom[,dat.otus.001per] #removed ASVs that occur less than 0.01%; increases the number of ASVs - includes more "microdiversity"
  dat.ra<-decostand(dat.01per, method = "total") #relative abundance of >1% taxa
  dfs_to_return <- list(as.data.frame(dat),as.data.frame(dat.dom),as.data.frame(dat.pa),
                        as.data.frame(dat.01per),as.data.frame(dat.001per),
                        as.data.frame(dat.ra))
  names(dfs_to_return) <- c("dat","dat.dom","dat.pa","dat.01per","dat.001per","dat.ra")
  return(dfs_to_return)
}

## insert your files into the function
# "feature.csv" = insert the name of your feature table .csv file
# "metadata.csv" = insert the name of your metadata .csv file

abund_metadata <- calculate_abund.metadata("feature.csv","metadata.csv") #if you have metadata
abund <- calculate_abund("feature.csv") #if you have NO metadata

# functions will return a list of data frames:
# dat = your original feature table that is transposed (or flipped)
# metadata = your metadata table             **dat and metadata should have the same samples!**
# dat.dom = the dominant ASVs found in dat
# dat.pa = the dat.dom table transformed into presence/absence data
# dat.01per = ASVs from dat that have an abundance of 0.1% or higher
# dat.001per = ASVs from dat that have an abundance of 0.01% or higher
# dat.ra = ASVs from dat.01per normalized into relative abundance

# To save any dataframe (df) within the function result list as a single df in
# your workspace, simply assign it to its own variable
# Example saving dat and metadata
feature <- abund_metadata$dat
metadta <- abund_metadata$metadata

###### Check for Batch Correction ######
## make sure you run the abundance function first so you have your abundance
## tables to input into the following functions (again RUN FUNCTIONS FIRST!)

## Batch correction
## first, create Bray-Curtis dissimilarity distance matrix
## you can either use the dat.01per or dat.ra tables as input,
## they give the same result

ra.bc.dist <- vegdist(dat.ra, method = "bray")
#OR dat01.bc.dist <- vegdist(dat.01per, method = "bray") if using dat01.per abundance

## next, determine if the batch effect is present and significant
# NOTE: if checking for the batch effect, you must have a column named
# "Batch" in your metadata how the batch is defined is dependent on your project.
# For example, if you had multiple individuals aid in numerous sequence runs
# then it will be good to base your "batches" on the sequence run that
# each sample was a part of.

## function
batch_test <-function(bray_dist_matrix, metadata){
  library(vegan)
  dis.Batch <- betadisper(bray_dist_matrix,metadata$Batch) # betadisper calculates dispersion (variances) within each group
  test <- permutest(dis.Batch, pairwise=TRUE, permutations=999) #determines if the variances differ by groups
  if (test$tab$`Pr(>F)`[1] <= 0.05){    #differences are SIGNIFICANT - use ANOSIM
    ano_sim <- anosim(bray_dist_matrix, metadata$Batch, permutations = 999)
    return(ano_sim)
  }
  else{            #differences are NOT SIGNIFICANT - use PERMANOVA (adonis))
    p_anova <- adonis2(bray_dist_matrix~metadata$Batch, permutations = 999)
    return(p_anova)
  }
}

## insert your input into the function
# "bray_dist_matrix" = insert the distance matrix you just created
# "metadata" = insert your metadata variable (may be already named 'metadata')
batch_test(ra.bc.dist, metadata)

## If p <= 0.05, then the batch effect was found to be significant and you MUST
## correct the batch effect BEFORE moving on to further analyses

## CONTINUE HERE IF YOUR DATA IS SIGNIFICANT FOR BATCH EFFECT!! SKIP IF NOT!
# Using a package called MMUPHin, we will have to adjust the data so that the
# batch effect is no longer affecting the statistical outcome of the data.

batch_correct <- function(feature_csv,metadata_csv){
  library(vegan)
  library(MMUPHin)
  dat<-t(data.matrix(read.csv(feature_csv, header=TRUE, row.names = 1)))
  metadata <- read.csv(metadata_csv, header = TRUE, row.names = 1)
  dat <- as.data.frame(dat)
  typeof(dat)
  common.rownames <- intersect(rownames(dat), rownames(metadata))
  dat <- dat[common.rownames,]
  metadata <- metadata[common.rownames,]
  #Adjusting (removing) batch effect
  fit_adjust_batch <- adjust_batch(feature_abd = t(dat), # ASVs should be rows in feature table (MATRIX)
                                   batch = "Batch",
                                   data = metadata)   # samples should be rows in metadata (DATAFRAME)
  feat_abd_adj <- fit_adjust_batch$feature_abd_adj #now adjusted feature table MATRIX
  feat_abd_adj <- as.data.frame(feat_abd_adj) #converting to data frame
  write.csv(feat_abd_adj, "feature_ADJUSTED.csv") #saving as csv
  dfs_to_return <- list(as.data.frame(dat),as.data.frame(metadata),
                        as.data.frame(feat_abd_adj))
  names(dfs_to_return) <- c("dat", "metadata","adj-feature")
  return(dfs_to_return)
}

## insert your files into the function
# "feature.csv" = insert the name of your feature table .csv file
# "metadata.csv" = insert the name of your metadata .csv file
# The adjusted feature table is saved as a csv under your working directory as "feature_ADJUSTED.csv"

batchadj <- batch_correct("feature.csv","metadata.csv")

# function will return a list of data frames:
# dat = your original feature table that is transposed (or flipped)
# metadata = your metadata table             **dat and metadata should have the same samples!**
# adj-feature = the adjusted feature table that is now batch corrected

# YOU NEED TO USE THE BATCH CORRECTED FEATURE TABLE FOR THE REST OF THE
# ANALYSES IF YOUR DATA WAS BATCH CORRECTED!!

###### Generating Rarefaction Curve #######
library(vegan)
## assigning the batch corrected feature table to its own variable; MAKE SURE TO READ ROWNAMES
rardat<-read.csv("feature_ADJUSTED.csv", header=TRUE, row.names=1, sep=',')

# samples are in columns and need to be in the rows so we need to flip or transpose the file
# transpose the data to rows; transposing changes the data frame to a matrix!
trans.rardat <- t(rardat)
# check file to make sure it worked
trans.rardat[1:5,1:5] #shows rows 1 through 5 and the samples should now be the rows
# assign the transformed data matrix into main data variable
rardat <- trans.rardat
# change back into data frame instead of matrix
rardat <-as.data.frame(rardat)
#check data file to make sure it looks okay
View(rardat)

rowSums(rardat) #sums the value of each row in the data frame;
# this shows the total sequencing reads for each sample

## Creating the rarefaction curve
# count the number of species within each sample
S <- specnumber(rardat)
raremax <- min(rowSums(rardat)) # takes the sample with the lowest number of sequencing reads

## Plotting the rarefaction curves
# ** auto removes samples that have no reads **

# creating color palette
col <- c("darkred", "forestgreen", "hotpink", "blue")  # feel free to edit the colors
                                                      # keep it between 3 and 5 colors
grp <- factor(sample(seq_along(col), nrow(rardat), replace = TRUE))
cols <- col[grp]

# creating rarefaction curve
# create the curve to estimate where the inflection point (the point at which
# most of the lines being to plateau) lies, then assign that value to the
# variable "inflection" below
rarecurve(rardat, step = 500, sample=raremax, col = cols, label = TRUE,
          main="Title", cex= 0.35, cex.axis= 0.95, cex.lab= 1, xlim=c(0,200000),
          xlab = "# of Sequencing Reads", ylab = "# of ASVs")
inflection <- 10000 # insert your estimated inflection point value
abline(0,0) # creates a vertical line at 0,0
abline(v = inflection, col="black", lwd=1.4) # creates the inflection line at specified value
                                             # feel free to change the line color and thickness (col and lwd)

###### Reproducing Abundance Data using Batch-corrected Data ######
### NOTE: If you ignore metadata, you will NOT be able to perform statistical
###       analyses on your data, only taxonomic analyses.

## RUN FUNCTIONS FIRST!
# use "calculate_abund.metadta" if you have your metadata file
# use "calculate_abund" if ignoring metadata (NO metadata file)
calculate_abund.metadata <- function(feature_csv,metadata_csv){
  library(vegan)
  dat<-t(data.matrix(read.csv(feature_csv, header=TRUE, row.names = 1)))
  #transposed so that the row names are now the sample names and the ASVs are the columns
  #BUT it is still a MATRIX
  metadata <- read.csv(metadata_csv, header = TRUE, row.names = 1)
  dat <- as.data.frame(dat)
  typeof(dat) #list
  common.rownames <- intersect(rownames(dat), rownames(metadata)) #returns the same row names found in each dataframe
  dat <- dat[common.rownames,] #subsets dat to only include the same row names as metadata
  write.csv(dat, "feature_ADJUSTED-Transposed_matched.csv")
  metadata <- metadata[common.rownames,] #subsets metadata to only include the same row names as dat
  write.csv(dat, "metadata_matched.csv")
  otu.abund<-which(colSums(dat)>2) #removes singletons and doubletons
  dat.dom<-dat[,otu.abund] #include dominant taxa
  dat.pa<-decostand(dat.dom, method ="pa") #turns dat.dom into presence/absence data (1/0)
  dat.otus.01per<-which(colSums(dat.pa) > (0.01*nrow(dat.pa)))
  dat.01per<-dat.dom[,dat.otus.01per] #removed ASVs that occur less than 0.1%
  write.csv(dat.01per, "feature_01percent.csv")
  dat.otus.001per<-which(colSums(dat.pa) > (0.001*nrow(dat.pa)))
  dat.001per<-dat.dom[,dat.otus.001per] #removed ASVs that occur less than 0.01%; increases the number of ASVs - includes more "microdiversity"
  write.csv(dat.001per, "feature_001percent.csv")
  dat.ra<-decostand(dat.01per, method = "total") #relative abundance of >1% taxa
  write.csv(dat.ra, "relative-abundance.csv")
  dfs_to_return <- list(as.data.frame(dat),as.data.frame(metadata),
                        as.data.frame(dat.dom),as.data.frame(dat.pa),
                        as.data.frame(dat.01per),as.data.frame(dat.001per),
                        as.data.frame(dat.ra))
  names(dfs_to_return) <- c("dat", "metadata","dat.dom","dat.pa","dat.01per","dat.001per","dat.ra")
  return(dfs_to_return)
}
calculate_abund <- function(feature_csv){
  library(vegan)
  dat<-t(data.matrix(read.csv(feature_csv, header=TRUE, row.names = 1)))
  #transposed so that the row names are now the sample names and the ASVs are the columns
  dat <- as.data.frame(dat)
  write.csv(dat, "feature_ADJUSTED-Transposed.csv")
  otu.abund<-which(colSums(dat)>2) #removes singletons and doubletons
  dat.dom<-dat[,otu.abund] #include dominant taxa
  dat.pa<-decostand(dat.dom, method ="pa") #turns dat.dom into presence/absence data (1/0)
  dat.otus.01per<-which(colSums(dat.pa) > (0.01*nrow(dat.pa)))
  dat.01per<-dat.dom[,dat.otus.01per] #removed ASVs that occur less than 0.1%
  write.csv(dat.01per, "feature_01percent.csv")
  dat.otus.001per<-which(colSums(dat.pa) > (0.001*nrow(dat.pa)))
  dat.001per<-dat.dom[,dat.otus.001per] #removed ASVs that occur less than 0.01%; increases the number of ASVs - includes more "microdiversity"
  write.csv(dat.001per, "feature_001percent.csv")
  dat.ra<-decostand(dat.01per, method = "total") #relative abundance of >1% taxa
  write.csv(dat.ra, "relative-abundance.csv")
  dfs_to_return <- list(as.data.frame(dat),as.data.frame(dat.dom),as.data.frame(dat.pa),
                        as.data.frame(dat.01per),as.data.frame(dat.001per),
                        as.data.frame(dat.ra))
  names(dfs_to_return) <- c("dat","dat.dom","dat.pa","dat.01per","dat.001per","dat.ra")
  return(dfs_to_return)
}

## insert your files into the function
# "feature.csv" = insert the name of your ADJUSTED feature table .csv file
# "metadata.csv" = insert the name of your metadata .csv file

abund_metadata <- calculate_abund.metadata("feature_ADJUSTED.csv","metadata.csv") #if you have metadata
abund <- calculate_abund("feature_ADJUSTED.csv") #if you have NO metadata

# functions will return a list of data frames (must assign to a variable in order to access them):
# dat = your original adjusted feature table that is transposed (or flipped)
# metadata = your metadata table             **dat and metadata should have the same samples!**
# dat.dom = the dominant ASVs found in dat
# dat.pa = the dat.dom table transformed into presence/absence data
# dat.01per = ASVs from dat that have an abundance of 0.1% or higher
# dat.001per = ASVs from dat that have an abundance of 0.01% or higher
# dat.ra = ASVs from dat.01per normalized into relative abundance

## NOTE: The dat, dat.01per, dat.001per, and dat.ra tables are saved as csvs into your working directory

# To save any dataframe (df) within the function result list as a single df in
# your workspace, simply assign it to its own variable
# Example saving dat and metadata
feature <- abund_metadata$dat
metadta <- abund_metadata$metadata
###### Taxonomy Analyses ######
## Creating a phyloseq object
# RUN FUNCTION FIRST!!
create_phyloseq <- function(abund, taxonomy, metadata) {
  library(phyloseq)
  library(microbiome)
  asvdat <- as.data.frame(t(abund)) #ASV/taxa should be the rows ("# of obs.")
  taxdat <- read.csv(taxonomy, header = TRUE, row.names = 1)
  meta <- read.csv(metadata, header = TRUE, row.names = 1)
  asvmat <- data.matrix(asvdat)
  taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
  ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
  TAX <- tax_table(taxmat)
  META <- sample_data(meta)
  #Merging metadata, taxonomy, and ASV tables into one phyloseq object
  physeq <- phyloseq(ASV,TAX,META)
  #Use transform functions from microbiome package
  transform <- microbiome::transform #converts into relative abundances
  #Merge rare taxa in to "Other"
  physeq_transform <- transform(physeq, "compositional")
  return(physeq_transform)
}
## insert your input into the function
# abundance = insert the abundance table you created (dat.ra or dat.01per)
# "taxonomy.csv" = insert your taxonomy .csv file
# "metadata.csv" = insert your metadata .csv file
physeq_transform <- create_phyloseq(abundance, "taxonomy.csv", "metadata.csv")


## Basic sequencing statistics
# Check number of sequencing reads observed in each sample
sample_sums(physeq)

# Calculate the total sequencing reads of your data
sum(sample_sums(physeq))

# Calculate the average number sequencing reads
mean(sample_sums(physeq))

# Find the lowest number of sequencing reads in your data
min(sample_sums(physeq))

# Find the highest number of sequencing reads in your data
max(sample_sums(physeq))

# Calculate the standard deviation of sequencing reads between samples
sd(sample_sums(physeq))

# Find the total amount of ASVs within your data
ntaxa(physeq)

## Retrieves the unique taxonomic ranks observed in the data set
## [#] = rank (starting from Domain and onward DPCOFGS)
get_taxa_unique(physeq, taxonomic.rank=rank_names(physeq)[7], errorIfNULL=TRUE)
#Unique Domains =
#Unique Phyla =
#Unique Classes =
#Unique Orders =
#Unique Families =
#Unique Genera =
#Unique Species =

## Aggregating by Taxonomic level
# This function allows you to aggregate the taxonomy based on Taxonomic
# level, gives you the counts of each level, and saves as a CSV file

# RUN FUNCTION FIRST!!
taxon_aggregate <- function(phyloseq_object,tax_level) {
  aggregation <- aggregate_taxa(phyloseq_object, tax_level) # aggregates the phyloseq object based on the taxonomic level provided
  counts <- as.data.frame(taxa_sums(aggregation)) # changes the aggregation to a table that can be exported and/or viewed
  colnames(counts)[1] <- paste(tax_level, "Count", sep="") # changes the name of first column in the table
  write.csv(counts, paste(tax_level, "Counts.csv", sep="_")) # exports the table as a csv file
  return(counts) # function will return the count table
}
## insert your input into the function
# phyloseq_object = insert the phyloseq object you created (physeq_transform)
# "tax_level" = insert the taxonomic level you will like to aggregate on
agg.tax.level <- taxon_aggregate(phyloseq_object,"tax_level")
# For example for Class, the code would read:
# Class <- taxon_aggregate(physeq_transform,"Class")


## Calculating the top taxonomic levels
# First, find the top taxa names and abundances
# This function allows you to specify the top taxonomy of your data and
# exports the resulting vector as a csv file

# RUN FUNCTION FIRST!!
top_taxa <- function(phyloseq_object, tax_level, n){
  top_taxanames <- sort(tapply(taxa_sums(phyloseq_object), tax_table(phyloseq_object)[, tax_level], sum), TRUE)[1:n]
  top_names <- as.data.frame(top_taxanames)
  colnames(top_names) <- paste(tax_level, "Abundance", sep=".")
  write.csv(top_names, paste("Top",as.character(n),tax_level, ".csv", sep="")) # exports the table as a csv file
  return(top_taxanames)
}
## insert your input into the function
# phyloseq_object = insert the phyloseq object you created (physeq_transform)
# "tax_level" = insert the taxonomic level you will like to aggregate on
# n = insert the number of taxa you want (e.g. top 5, top 10, etc.)
top.n.names <- top_taxa(phyloseq_object, "tax_level", n)
# For example for top 5 Classes, the code would read:
# top.class.names <- top_taxa(physeq_transform,"Class",5)


# Next, subset your phyloseq object to only include the top taxa you just specified
# cuts down the phyloseq object to only the top n
# Replace the "taxonomic level" with the level you want WITHOUT QUOTES!
top_tax <- subset_taxa(physeq_transform, "taxonomic level" %in% names(top.n.names))
# For example for Class, the code would read:
# top_Class <- subset_taxa(physeq_transform, Class %in% names(top.class.names))

## Creating stacked bar plots from your phyloseq object subset
# Change any taxonomic names to the level that you need
library(ggplot2)
topTaxplot <- plot_bar(top_tax, x="Sample", y="Abundance", fill="Class")
topTaxplot <- topTaxplot +
  geom_bar(aes(fill=Class, colour=Class), stat="identity", position="fill", width = 0.9) +   #width=0.96 removes any space between bars
  ggtitle("Top 5 Classes") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+   #vjust= moves the x-axis text labels
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5)) +   #hjust= 0.5 centers the title
  theme(legend.title = element_text(face="italic"))
topTaxplot # populates the graph into the "plot" panel to the right

###### Alpha Diversity - Measures ######
## Alpha diversity: the species richness that occurs within a given area within a region
## that is smaller than the entire distribution of the species (Moore, 2013)
## You can use the relative abundance or dat.01per data
adivmeasures <- function(abundance_data) {
  library(vegan)
  # Species richness: the number of species within a region
  S <- as.data.frame(specnumber(abundance_data))
  colnames(S)[1] ="S"
  #No. individuals:
  N <- as.data.frame(rowSums(abundance_data))
  colnames(N)[1] ="N"
  #Shannon-Weiner Diversity:
  ## Shannon index: a measure of the information content of a community rather than of the particular species
  ##               that is present (Moore, 2013) [species richness index]
  ## strongly influenced by species richness and by rare species (so sample size is negligible)
  H <- as.data.frame(diversity(abundance_data), index="shannon")
  colnames(H)[1] ="H"
  #Pielou's Evenness:
  ## Pielou's evenness: an index that measures diversity along with the species richness
  ## Formula - J = H/log(S) (aka Shannon evenness index)
  ## evenness = the count of individuals of each species in an area; 0 is no evenness & 1 is complete evenness
  J = H/log(S)
  colnames(J)[1] ="J"
  #Simpson's Diversity (1/D) (inverse):
  ## gives the Simpson index the property of increasing as diversity increases (the dominance of
  ## a few species decreases)
  inv.D <- as.data.frame(diversity(abundance_data, index="inv"))
  colnames(inv.D)[1] ="inv.D"
  #Combine data together into a single new data frame, export as CSV
  diversitybysample <- cbind(S, N, H, J,inv.D)
  write.csv(diversitybysample, "AlphaDiversity.csv")
  return(diversitybysample)
}
## This function returns one dataframe that contains the alpha diversity measure metrics for each sample


# Merge diversity with metadata table and export as csv
# DONE IN EXCEL BEFORE RUNNING NEXT CODE: rename the first column of both data
# tables as the SAME NAME. e.g. both should have the sample column labelled as "Sample"
diversitybysample <- read.csv("AlphaDiversity.csv", row.names = 1)
met <- read.csv("metadata.csv", row.names = 1)
adivmet <- cbind(diversitybysample,met)
write.csv(adivmet,"Metadata-Diversity.csv")
###### Alpha Diversity Statistics ######
library(vegan)

# load in your metadata that has the alpha diversity indices included
metadata <- read.csv("Metadata-Diversity.csv", header = TRUE, row.names = 1)

#### Testing Statistical Significance
## Normality - Shapiro Test (only done on NUMERIC data)
## p <= 0.05 = H0 REJECTED -> DATA IS NOT NORMAL
## p > 0.05 = H0 ACCEPTED -> DATA IS NORMAL

## If data is NOT normal the first time, try transforming the data using log
## and sqrt and retest for normality.

#Alpha Diversity Variables - Test for Normality
shapiro.test(metadata$S)
shapiro.test(metadata$N)
shapiro.test(metadata$H)
shapiro.test(metadata$J)
shapiro.test(metadata$inv.D)

### CONTINUE HERE IF DATA IS NORMAL (OR TRANSFORMATIONS NORMALIZED THE DATA)
# ANOVA: Parametric Data (normal)
# Tukey Test - calculate pairwise comparisons between group levels

# RUN FUNCTION FIRST!!
p_anova <- function(adiv_var, test_var) {
  library(pgirmess)
  library(multcompView)
  anova_res <- aov(adiv_var ~ test_var) # performs an ANOVA
  summary <- summary(anova_res) # provides an ANOVA table with p-values
  Tukey <- TukeyHSD(anova_res) # pairwise comparison
  return_items <- list(summary, Tukey)
  names(return_items) <- c("ANOVA summary", "TukeyTest")
  return(return_items)
}
## insert your files into the function
# "adiv_var" = insert the alpha diversity metric you want to test from the metadata file
#              (metadata$S, metadata$N, metadata$J, metadata$H, or metadata$inv.D)
# "test_var" = insert the variable(s) you want to test against the alpha diversity
#              metric (e.g. year (metadata$year), month, year and month (metadata$year+month), etc.)
#         ** MAKE SURE YOUR VARIABLE IS A FACTOR! IF NOT THEN SURROUND IT WITH as.factor() **
p_anova(adiv_var, test_var) # or p_anova(adiv_var, as.factor(test_var))

# Example inputs:
S_anova <- p_anova(metadata$S, as.factor(metadata$Year))
N_anova <- p_anova(metadata$N, metadata$Year)
J_anova <- p_anova(metadata$J, metadata$Year)
H_anova <- p_anova(metadata$H, metadata$Year)
inv.D_anova <- p_anova(metadata$inv.D, metadata$Year)

### CONTINUE HERE IF DATA IS NOT NORMAL AND TRANSFORMATIONS DID NOT WORK
library(pgirmess)
library(multcompView)

# Kruskal Wallis: Nonparametric Data (not normal)
# Pairwise Wilcox Test - calculate pairwise comparisons between group levels
#                        with corrections for multiple testing (non-parametric)

# RUN FUNCTION FIRST!!
nonp_kruskal <- function(adiv_var, test_var) {
  library(pgirmess)
  library(multcompView)
  kruskal.test(adiv_var ~ test_var)
  pair_WilTest <- pairwise.wilcox.test(adiv_var, test_var, p.adjust.method = "fdr") #pairwise comparisons between the variable levels
  kmc <- kruskalmc(adiv_var ~ test_var) # multiple-comparison test
  # comparisons TRUE= significantly different or FALSE= not significantly different
  # To look for homogeneous groups, and give each group a code (letter):
  test <- kmc$dif.com$stat.signif # select logical vector
  names(test) <- row.names(kmc$dif.com)# add comparison names
  # create a list with "homogeneous groups" coded by letter
  let <- multcompLetters(test, compare="<", threshold=0.05,
                         Letters=c(letters, LETTERS, "."),reversed = FALSE)
  # significant letters for the multiple comparison test
  # if the letter are the SAME, then no significant differences were found
  # between those variables
  returned_items <- list(pair_WilTest,kmc,let)
  names(returned_items) <- c("pairwise", "multiComp","letter-comparisons")
  return(returned_items)
}
## insert your files into the function
# "adiv_var" = insert the alpha diversity metric you want to test from the metadata file
#              (metadata$S, metadata$N, metadata$J, metadata$H, or metadata$inv.D)
# "test_var" = insert the variable(s) you want to test against the alpha diversity
#              metric (e.g. year (metadata$year), month, year and month (metadata$year+month), etc.)
#         ** MAKE SURE YOUR VARIABLE IS A FACTOR! IF NOT THEN SURROUND IT WITH as.factor() **
nonp_kruskal(adiv_var, test_var) # or nonp_kruskal(adiv_var, as.factor(test_var))

# Example inputs:
S_krusk <- nonp_kruskal(metadata$S, metadata$Year)
N_krusk <- nonp_kruskal(metadata$N, metadata$Year)
J_krusk <- nonp_kruskal(metadata$J, metadata$Year)
H_krusk <- nonp_kruskal(metadata$H, metadata$Year)
inv.D_krusk <- nonp_kruskal(metadata$inv.D, metadata$Year)


### Plotting boxplots of alpha diversity by specified variable
## NOTES: You can replace "Year" with your specified variable
##        Adding text to your graph is OPTIONAL but if you are adding it then
##        you'll have to play around with their coordinates, what they say, size, and color

# Creating pdf for the plots to populate
pdf("AlphaDiverisityPlots.pdf")
# plot each boxplot on its own page
par(mar=c(5,6,2,2)+0.1)
boxplot(S~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="Species Richness (S)", line=4.25, cex.lab=1.15)
text(y=1500, x=3, labels="b", col="blue", cex=1.2)
text(y=1420, x=2, labels="a", col="red", cex=1.2)        # labeling which groups are significantly different than the other
text(y=1585, x=1, labels="a", col="red", cex=1.2)

par(mar=c(5,4.5,2,2)+0.1)
boxplot(H~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="Shannon Diversity Index (H)", line=2.8, cex.lab=1.15)
text(y=4, x=3, labels="b", col="blue", cex=1.2)
text(y=3.4, x=2, labels="a", col="red", cex=1.2)
text(y=3.6, x=1, labels="ab", col="purple", cex=1.2)

boxplot(J~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="Species Evenness (J)", line=3, cex.lab=1.15)
text(y=0.73, x=3, labels="b", col="blue", cex=1.2)
text(y=0.685, x=2, labels="ab", col="purple", cex=1.2)
text(y=0.73, x=1, labels="a", col="red", cex=1.2)

par(mar=c(5,6,2,2)+0.1)
boxplot(inv.D~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="inverse Simpson Diversity Index (inv.D)", line=3.6, cex.lab=1.15)
text(y=440, x=3, labels="a", col="red", cex=1.2)
text(y=420, x=2, labels="b", col="blue", cex=1.2)
text(y=340, x=1, labels="a", col="red", cex=1.2)

boxplot(N~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="No. of Individuals (N)", line=4.25, cex.lab=1.15)
text(y=90000, x=3, labels="b", col="blue", cex=1.2)
text(y=130000, x=2, labels="a", col="red", cex=1.2)
text(y=160000, x=1, labels="a", col="red", cex=1.2)

# stop saving to pdf
dev.off()
###### Beta Diversity - Creating Distance Matrix ######
library(vegan)
bc.dist <- vegdist(dat.ra, method = "bray") # you can use either dat.ra or dat.01per to produce the distance matrix
###### Beta Diversity - Statistics ######
# RUN FUNCTION FIRST!!
betadiv_stats <-function(bray_dist_matrix, metadata){
  library(vegan)
  dis.Variable <- betadisper(bray_dist_matrix,metadata) # betadisper calculates dispersion (variances) within each group
  test <- permutest(dis.Variable, pairwise=TRUE, permutations=999) #determines if the variances differ by groups
  if (test$tab$`Pr(>F)`[1] <= 0.05){    #differences are SIGNIFICANT - use ANOSIM
    ano_sim <- anosim(bray_dist_matrix, metadata, permutations = 999)
    return(ano_sim)
  }
  else{            #differences are NOT SIGNIFICANT - use PERMANOVA (adonis))
    p_anova <- adonis2(bray_dist_matrix~metadata, permutations = 999)
    return(p_anova)
  }
}
## insert your input into the function
# "bray_dist_matrix" = insert the distance matrix you created
# "metadata$Variable" = insert your metadata variable with the Variable you want to test (may be already named 'metadata')
betadiv_stats(bc.dist,metadata$Variable)

###### Beta Diversity - nMDS plots ######
library(vegan)
# Creating nMDS plot - 2D
nmds2d <- metaMDS(bc.dist,k=2,autotransform = F,trymax=20)
nmds2d
#Dimensions = 2
#Stress =
stressplot(nmds2d)
#Shepard plot shows scatter around the regression between the inter-point
#distances in the final configuration (i.e., the distances between each pair of communities)
#against their original dissimilarities
nmds.plot <- ordiplot(nmds2d,display="sites") #populates the plot into a variable
## Adding ellipses to group years
ordihull(nmds.plot,groups=metadata$Variable,draw="lines",col=c("tomato3","steelblue3","springgreen3")) # adds ellipses around the point groups (OPTIONAL!)
##adjust colors to match each year, pch=20 makes it bullet points
points(nmds.plot,"sites", pch=20, col= "tomato4", select = metadata$Variable == "Level 1")     # set the metadata to the variable you need it to be
points(nmds.plot,"sites", pch=20, col= "steelblue4", select = metadata$Variable == "Level 2")  # and set a different color to each level of the variable
points(nmds.plot,"sites", pch=20, col= "springgreen4", select = metadata$Variable == "Level 3")
##Add Stress Value
text(1.2,1.5,"2D Stress: ", cex=0.9) # make sure you add the stress value in an empty portion of the graph
##Adding legend
legend("topleft",legend= c("Level 1","Level 2", "Level 3"),   # customize the legend to match the colors and variables
       title = "Year",
       col=c("tomato4","steelblue4","springgreen4"),
       pch=19, cex=1)
##Adding title
title(main="nMDS of Relative Abundances by Variable") # adds a title to the graph

## If making plots for multiple variables, you'll have to redo the point
## customization and legend customization to match the variable
###### CCA analysis & plots ######
library(vegan)
# you can use the dat.ra table or the dat.01per table
# input the quantitative variables (columns that are numeric) from the metadata
ccamodel <- cca(dat.ra~., metadata[,c(7:37)])
finalmodel<- ordistep(ccamodel, scope=formula(ccamodel))
vif.cca(finalmodel) ## values should be under 10
# If VIF>10, the variable presents colinearity with another or other variables.
# In that case, delete the variable from initial dataset and redo the analysis.
# VIF = 1 for completely independent variables,and values above 10 or 20
# (depending on your taste) are regarded as highly multicollinear (dependent on other variables).

# Test the significance of the entire CCA model
anova.cca(finalmodel) #should be significant! (p < 0.05)
#           Df ChiSquare      F Pr(>F)
# Model     18    1.6290 5.7307  0.001 ***
# Residual 522    8.2434
# ---

finalmodel
## Note that "Total Inertia" is the total variance in species (observations matrix) distributions.
## "Constrained Inertia" is the variance explained by the environmental variables (gradients matrix).
## The "Proportion" values represent the percentages of variance of species distributions explained
## by Constrained (environmental) and Unconstrained variables. Eigenvalues of constrained and
## unconstrained axes represent the amount of variance explained by each CCA axis (graphs usually
## present the first two constrained axes, so take a look at their values).

#Total Inertia = total variance in species (observed distributions)
#Unconstrained Inertia = the variance explained by the environmental variables

#               Inertia Proportion Rank
# Total           9.872      1.000
# Constrained     1.629      0.165   18
# Unconstrained   8.243      0.835  522
# Inertia is scaled Chi-square

R2.adj.cca <- RsquareAdj(finalmodel)
# adjusting the R-squared value: The adjusted R2 tells you the percentage of
# variation explained by only the independent variables that actually affect
# the dependent variable; also indicates how well terms fit a curve or line,
# but adjusts for the number of terms in a model
R2.adj.cca
# r.squared: 0.173352
# adj.r.squared: 0.1446893  -> the R2 value you document

summary(finalmodel)

## Correlation between the variables within the model
# Creates pairs plot to see the correlation statistics between each variable
library(psych)
pairs.panels(metadata[,c(7:19,21:24,33)]) # use the metadata columns that are
                                          # in your final CCA model

### Creating CCA plots
cca.p <- plot(finalmodel,type = "none")

# Fitting of the environmental variables to the CCA plot
ef.cca<- envfit(cca.p,metadata[,c(7:10,12:19,21,23,24,31,33)]) # use the metadata columns that are
                                                               # in your final CCA model

# Creating R2 threshold for vectors
# Function: select.envfit - Setting r2 cutoff values to display in an
#                           ordination.
# function (select.envfit) filters the resulting list of function (envfit) based on their p values. This allows to display only significant values in the final plot.
# RUN FUNCTION FIRST!!
select.envfit<-function(fit, r.select){ #needs two sorts of input: fit= result of envfit, r.select= numeric, correlation minimum threshold
  for (i in 1:length(fit$vectors$r)) { #run for-loop through the entire length of the column r in object fit$vectors$r starting at i=1
    if (fit$vectors$r[i]<r.select) { #Check whether r<r.select, i.e. if the correlation is weaker than the threshold value. Change this Parameter for r-based selection
      fit$vectors$arrows[i,]=NA #If the above statement is TRUE, i.e. r is smaller than r.select, then the coordinates of the vectors are set to NA, so they cannot be displayed
      i=i+1 #increase the running parameter i from 1 to 2, i.e. check the next value in the column until every value has been checked
    }
  }
  return(fit)
}

# Running select function on actual data
ef.cca<- select.envfit(ef.cca, 0.3) # only includes significant variables
                                    # with an R2 of 0.3 or higher

# Setting up base plot
par(mar=c(5.1, 6.1, 3.1, 4.1))
plot(finalmodel,type = "none")
abline(h = 0, v = 0, col = "white", lwd = 2)
box()

# EDIT POINTS AS IT MATCHES YOUR ANALYSIS!
# Adding the points
points(cca.p,"sites", pch=19, col= "goldenrod3", select = metadata$Variable == "Level 1")
points(cca.p,"sites", pch=19, col= "mediumpurple2", select = metadata$Variable == "Level 2")
points(cca.p,"sites", pch=19, col= "springgreen4", select = metadata$Variable == "Level 3")
# Plotting envfit vectors
plot(ef.cca, col = "black", p.max=0.05)
# Add legend & Title
legend(locator(1),legend=c("Level 1","Level 2", "Level 3"),
       col=c("goldenrod3","mediumpurple2", "springgreen4"), pch=19, cex=1.2,
       title = "Variable")
# locator(1) allows you to choose a place to place the legend within your plot
title(main="Title")





