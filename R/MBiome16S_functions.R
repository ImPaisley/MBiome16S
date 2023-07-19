
#' @title Set Working Directory and Seed
#' @description
#' This function allows the user to set their working directory and seed using
#' only one function. This is a required step when performing any kind of
#' analyses in RStudio to allow for reproducible results.
#' @param file_path a string input, a file path to the folder where all generated files will be stored
#' @param seed_number a random number that is user-defined
#' @return This function does not return anything
#' @export
#' @examples
#' setwd_seed("C:/Users/Public/Documents",2023)
setwd_seed <- function(file_path,seed_number) {
  setwd(file_path)
  set.seed(seed_number)
}

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

batch_test <-function(abundance_data, metadata){
  library(vegan)
  bc.dist <- vegdist(abundance_data, method = "bray") # creates distance matrix; you can either use the dat.01per or dat.ra tables as input
  dis.Batch <- betadisper(bc.dist,metadata$Batch) # betadisper calculates dispersion (variances) within each group
  test <- permutest(dis.Batch, pairwise=TRUE, permutations=999) #determines if the variances differ by groups
  if (test$tab$`Pr(>F)`[1] <= 0.05){    #differences are SIGNIFICANT - use ANOSIM
    ano_sim <- anosim(bc.dist, metadata$Batch, permutations = 999)
    return(ano_sim)
  }
  else{            #differences are NOT SIGNIFICANT - use PERMANOVA (adonis))
    p_anova <- adonis2(bc.dist~metadata$Batch, permutations = 999)
    return(p_anova)
  }
}

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

calculate_abund.metadata_ADJ <- function(feature_csv,metadata_csv){
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

calculate_abund_ADJ <- function(feature_csv){
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

taxon_aggregate <- function(phyloseq_object,tax_level) {
  aggregation <- aggregate_taxa(phyloseq_object, tax_level) # aggregates the phyloseq object based on the taxonomic level provided
  counts <- as.data.frame(taxa_sums(aggregation)) # changes the aggregation to a table that can be exported and/or viewed
  colnames(counts)[1] <- paste(tax_level, "Count", sep="") # changes the name of first column in the table
  write.csv(counts, paste(tax_level, "Counts.csv", sep="_")) # exports the table as a csv file
  return(counts) # function will return the count table
}

top_taxa <- function(phyloseq_object, tax_level, n){
  top_taxanames <- sort(tapply(taxa_sums(phyloseq_object), tax_table(phyloseq_object)[, tax_level], sum), TRUE)[1:n]
  top_names <- as.data.frame(top_taxanames)
  colnames(top_names) <- paste(tax_level, "Abundance", sep=".")
  write.csv(top_names, paste("Top",as.character(n),tax_level, ".csv", sep="")) # exports the table as a csv file
  return(top_taxanames)
}

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

betadiv_stats <-function(abundance_data, metadata){
  library(vegan)
  bc.dist <- vegdist(abundance_data, method = "bray") # you can use either dat.ra or dat.01per to produce the distance matrix
  dis.Variable <- betadisper(bc.dist,metadata) # betadisper calculates dispersion (variances) within each group
  test <- permutest(dis.Variable, pairwise=TRUE, permutations=999) #determines if the variances differ by groups
  if (test$tab$`Pr(>F)`[1] <= 0.05){    #differences are SIGNIFICANT - use ANOSIM
    ano_sim <- anosim(bc.dist, metadata, permutations = 999)
    return(list(c(bc.dist,ano_sim)))
  }
  else{            #differences are NOT SIGNIFICANT - use PERMANOVA (adonis))
    p_anova <- adonis2(bc.dist~metadata, permutations = 999)
    return(list(c(bc.dist,p_anova)))
  }
}

select.envfit<-function(fit, r.select){ #needs two sorts of input: fit= result of envfit, r.select= numeric, correlation minimum threshold
  for (i in 1:length(fit$vectors$r)) { #run for-loop through the entire length of the column r in object fit$vectors$r starting at i=1
    if (fit$vectors$r[i]<r.select) { #Check whether r<r.select, i.e. if the correlation is weaker than the threshold value. Change this Parameter for r-based selection
      fit$vectors$arrows[i,]=NA #If the above statement is TRUE, i.e. r is smaller than r.select, then the coordinates of the vectors are set to NA, so they cannot be displayed
      i=i+1 #increase the running parameter i from 1 to 2, i.e. check the next value in the column until every value has been checked
    }
  }
  return(fit)
}
