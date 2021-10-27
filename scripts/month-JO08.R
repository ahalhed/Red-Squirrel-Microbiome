#---
#title: "PCNM for Squirrel Microbiome by Month (SHARCNET)"
#author: "Alicia Halhed"
#date: "10/04/2021"
#---

print("Set up (working directory, theme, and packages)")
# set working directory
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/Red-Squirrel-Microbiome")

# attach required packages
library(qiime2R)
library(phyloseq)
library(vegan)
library(zCompositions)
# devtools::install_github('ggloor/CoDaSeq/CoDaSeq')
library(CoDaSeq)
library(tidyverse)

# set theme for plots
theme_set(theme_bw())

print("Initiate functions for analysis")
# subset the metadata by grid and year
met_year <- function(metadata, grid, year) {
  met <- rownames_to_column(metadata, var = "SampleID")
  df1 <- subset(met, Grid == grid)
  df2 <- subset(df1, Year == year)
  df3 <- column_to_rownames(remove_rownames(df2), var = "SampleID")
  return(df3)
}
# subset the XY's by grid, year, and month
XY_month <- function(metadata, grid, year, month) {
  met <- rownames_to_column(metadata, var = "SampleID")
  df1 <- subset(met, Grid == grid, 
                select = c("SampleID", "Location.X", "Location.Y", "Year", "Month"))
  df2 <- subset(df1, Year == year, 
                select = c("SampleID", "Location.X", "Location.Y", "Month"))
  df3 <- subset(df2, Month == month, 
                select = c("SampleID", "Location.X", "Location.Y"))
  df4 <- column_to_rownames(remove_rownames(df3), var = "SampleID")
  return(df4)
}
# maximum distance
max_dist <- function(dm) {
  df1 <- as.data.frame(as.matrix(dm))
  # message: `summarise_each_()` is deprecated as of dplyr 0.7.0. (use across)
  summ <- summarise_each(df1, ~ max(df1, na.rm=TRUE))
  m <- apply(summ, 1, max)
  return(m)
}
# community object
comm_obj <- function(XY, c) {
  # subset the OTUs (c is OTU table being subset)
  comm <- c %>%
    subset(., rownames(.) %in% rownames(XY)) %>%
    .[ , colSums(.)>0 ]
  return(comm)
}
# subset full metadata by month
met_month <- function(XY, meta) {
  # subset for samples in the XY
  df1 <- subset(meta, rownames(meta) %in% rownames(XY))
  # remove NAs
  df2 <- df1 %>% select_if(~ !any(is.na(.)))
  # remove location information
  df3 <- select(df2, -c('Location.X', 'Location.Y', Location))
  # select columns with more than one level
  df4 <- df3[sapply(df3, function(x) length(unique(x))>1)]
  return(df4)
}

# get the data
print("Read in the Data")
print("Building phyloseq object")
ps <- qza_to_phyloseq(features = "../ASV-table-10-filtered.qza",
                      metadata = "../input/RS_meta.tsv") %>%
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.))
# based on the meta function from the microbiome package
# I don't want to load a whole package for one function
print("Read in the metadata")
rs_q2_metadata <- as(sample_data(ps), "data.frame")
rownames(rs_q2_metadata) <- sample_names(ps)

# example in https://github.com/ggloor/CoDaSeq/blob/master/Intro_tiger_ladybug.Rmd
print("Aitchison transformation")
# rows are OTUs, then transposed to OTUs as column
# impute the OTU table
OTUimp <- cmultRepl(otu_table(ps), label=0, method="CZM", output="p-counts") # all OTUs

# compute the aitchison values with absolute values
OTUclr <- codaSeq.clr(abs(OTUimp))

## Core and non-core divide
print("Extract Core")
# core OTUs
cOTU <- read.csv("data/core.csv") %>%
  # get the OTUs identified as core contributors to beta diversity
  .[which(.$fill == "core"),]

# make the new data frames
print("Subset the OTU table to find core and non-core OTUs")
OTU_core <- select(as.data.frame(OTUclr), one_of(cOTU$otu))
OTU_nc <- select(as.data.frame(OTUclr), -one_of(cOTU$otu))

## XY data
print("Accessing the XY data by month")
# start by getting the metadata for the grid/year
met <- met_year(rs_q2_metadata, "JO", 2008)
# loop to create individual month data frames
for (month in unique(met$Month)) {
  df <- XY_month(rs_q2_metadata, "JO", 2008, month)
  assign(paste('Month',month,sep = ' '),df)
  rm(df, month)
}
# make a list of the data frames generated from the loop
XY_list <- do.call("list",
                   # searching the global environment for the pattern
                   mget(grep("Month", names(.GlobalEnv), value=TRUE)))
print("Months in this Grid")
unique(met$Month)
print("Computing Euclidean Distances")
dist_list <- lapply(XY_list, dist)
print("Maximum Euclidean Distance by Month")
lapply(dist_list, max_dist)
## community objects
# subset the samples from the core microbiome
print("Build the community object (OTU table) for grid/year/month")
commFull <- lapply(XY_list, comm_obj, c=OTUclr)
commCore <- lapply(XY_list, comm_obj, c=OTU_core)
commNC <- lapply(XY_list, comm_obj, c=OTU_nc)
## metadata
# sample ID's are rownames
# get the metadata subset
print("Extract metadata for grid/year/month")
# extracting only non-location metadata columns
# with more than one unique level
met_list <- lapply(XY_list, met_month, meta=rs_q2_metadata)

## Clean up 1
# Remove objects we're done with
print("Removing pobjects that are no longer needed")
rm(rs_q2_metadata, cOTU, OTUclr, ps, OTU_core, OTU_nc)
# not removing the Months because they're small and not the same across all

## Analysis time!
# unweighted PCNM
print("Unweighted PCNM - for use with all OTU tables")
pcnm_list <- lapply(dist_list, pcnm)
# need to sort out how to generalize printing  month # with the vectors
for (month in pcnm_list) {
  print(month$vectors)
}
print("Acessing PCNM scores")
scores_list <- lapply(pcnm_list, scores)


# core OTUs
print("Analysis for Core OTUs")
print("Variance partitioning - Core OTUs")
vp_mod1_list <- mapply(varpart, commCore, scores_list, data=met_list, 
                       MoreArgs = list(~.),
                       SIMPLIFY = FALSE)
vp_mod1_list
# plot the partitioning
pdf(file = "./plots/core_JO2008_vp_mod1M.pdf")
# make plot
# plotted in numerical order by month
lapply(vp_mod1_list, plot)
dev.off()
#remove vp object, to repeat with new OTU table
rm(vp_mod1_list)

# test with RDA
print("Testing with RDA (full model) - core OTUS")
# create a tiny anonymous function to include formula syntax in call
abFrac <- mapply(function(x,data) rda(x~., data), 
                 commCore, met_list, SIMPLIFY=FALSE)
abFrac # Full model
lapply(abFrac, anova, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
lapply(abFrac, RsquareAdj)

# Test fraction [a] using partial RDA:
print("Testing with partial RDA (fraction [a]) - core OTUS")
# create a tiny anonymous function to include formula syntax in call
aFrac <- mapply(function(x,y,data) rda(x~.+Condition(scores(y)), data), 
                commCore, pcnm_list, met_list, SIMPLIFY=FALSE)
aFrac
lapply(aFrac, anova, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
lapply(aFrac, RsquareAdj)

# # forward selection for parsimonious model
# print("Forward selection for parsimonious model - core OTUs")
# # env variables
# print("Environmental variables - core OTUs")
# # create a tiny anonymous function to include formula syntax in call
# abFrac0 <- mapply(function(x,data) rda(x~1, data), 
#                   commCore, met_list, SIMPLIFY=FALSE) # Reduced model
# step.env <- mapply(function(x,y) ordiR2step(x, scope = formula(y)), 
#                    abFrac0, abFrac, SIMPLIFY=FALSE)
# step.env # an rda model, with the final model predictor variables
# 
# print("Summary of environmental selection process - core OTUs")
# lapply(step.env, function(x) x$anova)
# print("ANOVA on full environmental selection - core OTUs")
# lapply(step.env, anova)
# 
# # save plot
# pdf(file = "./plots/core_JO2008_step_envM.pdf")
# # make plot
# lapply(step.env, plot)
# dev.off()

# spatial variables
print("Spatial variables - core OTU")
pcnm_df <- lapply(pcnm_list, function(x) as.data.frame(scores(x)))
bcFrac <- mapply(function(x,data) rda(x~., data), 
                 commCore, pcnm_df, SIMPLIFY=FALSE) # Full model
bcFrac0 <- mapply(function(x,data) rda(x~1, data), 
                  commCore, pcnm_df, SIMPLIFY=FALSE) # Reduced model
step.space <- mapply(function(x,y) ordiR2step(x, scope = formula(y)), 
                     bcFrac0, bcFrac, SIMPLIFY=FALSE)
step.space

print("Summary of spatial selection process - core OTU")
lapply(step.space, function(x) x$anova)
print("ANOVA on full spatial selection - core OTU")
lapply(step.space, anova)

# save plot
pdf(file = "./plots/core_JO2008_step_spaceM.pdf")
# make plot
lapply(step.space, plot)
dev.off()

print("Partition Bray-Curtis dissimilarities - core OTUs")
vdist <- lapply(commCore, vegdist)
pbcd <- mapply(function(x,y,z) varpart(x, ~., y, data = z),
               vdist, scores_list, met_list, SIMPLIFY=FALSE)
pbcd

#cleanup
# remove objects to be replaced/no longer needed
rm(vdist,pbcd, commCore)
rm(abFrac, aFrac,abFrac0, pcnm_df, bcFrac, bcFrac0, step.space) #step.env, 


# non-core OTUs
print("Analysis for non-core OTUs")
print("Variance partitioning - non-core OTUs")
vp_mod1_list <- mapply(varpart, commNC, scores_list, data=met_list, 
                       MoreArgs = list(~.),
                       SIMPLIFY = FALSE)
vp_mod1_list
# plot the partitioning
pdf(file = "./plots/nc_JO2008_vp_mod1M.pdf")
# make plot
# plotted in numerical order by month
lapply(vp_mod1_list, plot)
dev.off()
#remove vp object, to repeat with new OTU table
rm(vp_mod1_list)

# test with RDA
print("Testing with RDA (full model) - non-core OTUS")
# create a tiny anonymous function to include formula syntax in call
abFrac <- mapply(function(x,data) rda(x~., data), 
                 commNC, met_list, SIMPLIFY=FALSE)
abFrac # Full model
# anova
lapply(abFrac, anova, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
lapply(abFrac, RsquareAdj)

# Test fraction [a] using partial RDA:
print("Testing with partial RDA (fraction [a]) - non-core OTUS")
# create a tiny anonymous function to include formula syntax in call
aFrac <- mapply(function(x,y,data) rda(x~.+Condition(scores(y)), data), 
                commNC, pcnm_list, met_list, SIMPLIFY=FALSE)
aFrac
# anova
lapply(aFrac, anova, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
lapply(aFrac, RsquareAdj)

# # forward selection for parsimonious model
# print("Forward selection for parsimonious model - non-core OTUs")
# # env variables
# print("Environmental variables - non-core OTUs")
# # create a tiny anonymous function to include formula syntax in call
# abFrac0 <- mapply(function(x,data) rda(x~1, data), 
#                   commNC, met_list, SIMPLIFY=FALSE) # Reduced model
# 
# step.env <- mapply(function(x,y) ordiR2step(x, scope = formula(y)), 
#                    abFrac0, abFrac, SIMPLIFY=FALSE)
# step.env # an rda model, with the final model predictor variables
# 
# print("Summary of environmental selection process - non-core OTUs")
# lapply(step.env, function(x) x$anova)
# print("ANOVA on full environmental selection - non-core OTUs")
# lapply(step.env, anova)
# 
# # save plot
# pdf(file = "./plots/nc_JO2008_step_envM.pdf")
# # make plot
# lapply(step.env, plot)
# dev.off()

# spatial variables
print("Spatial variables - non-core OTU")
pcnm_df <- lapply(pcnm_list, function(x) as.data.frame(scores(x)))
bcFrac <- mapply(function(x,data) rda(x~., data), 
                 commNC, pcnm_df, SIMPLIFY=FALSE) # Full model
bcFrac0 <- mapply(function(x,data) rda(x~1, data), 
                  commNC, pcnm_df, SIMPLIFY=FALSE) # Reduced model
step.space <- mapply(function(x,y) ordiR2step(x, scope = formula(y)), 
                     bcFrac0, bcFrac, SIMPLIFY=FALSE)
step.space

# summary of selection process
print("Summary of spatial selection process - non-core OTU")
lapply(step.space, function(x) x$anova)
print("ANOVA on full spatial selection - non-core OTU")
lapply(step.space, anova)

# save plot
pdf(file = "./plots/nc_JO2008_step_spaceM.pdf")
# make plot
lapply(step.space, plot)
dev.off()

print("Partition Bray-Curtis dissimilarities - non-core OTUs")
vdist <- lapply(commNC, vegdist)
pbcd <- mapply(function(x,y,z) varpart(x, ~., y, data = z),
               vdist, scores_list, met_list, SIMPLIFY=FALSE)
pbcd

#cleanup
# remove objects to be replaced
rm(vdist, pbcd, commNC)
rm(abFrac, aFrac,abFrac0, pcnm_df, bcFrac, bcFrac0, step.space) #step.env, 

# analysis for all OTUs
print("Analysis for All OTUs")
print("Variance partitioning - All OTUs")
vp_mod1_list <- mapply(varpart, commFull, scores_list, data=met_list, 
                       MoreArgs = list(~.),
                       SIMPLIFY = FALSE)
vp_mod1_list
# plot the partitioning
pdf(file = "./plots/JO2008_vp_mod1M.pdf")
# make plot
# plotted in numerical order by month
lapply(vp_mod1_list, plot)
dev.off()
#remove vp object, done with it now
rm(vp_mod1_list)

# test with RDA
print("Testing with RDA (full model) - all OTUS")
# create a tiny anonymous function to include formula syntax in call
abFrac <- mapply(function(x,data) rda(x~., data), 
                 commFull, met_list, SIMPLIFY=FALSE)

abFrac # Full model

lapply(abFrac, anova, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
lapply(abFrac, RsquareAdj)

# Test fraction [a] using partial RDA:
print("Testing with partial RDA (fraction [a]) - all OTUS")
# create a tiny anonymous function to include formula syntax in call
aFrac <- mapply(function(x,y,data) rda(x~.+Condition(scores(y)), data), 
                commFull, pcnm_list, met_list, SIMPLIFY=FALSE)
aFrac
# anova
lapply(aFrac, anova, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
lapply(aFrac, RsquareAdj)

# forward selection for parsimonious model
print("Forward selection for parsimonious model - all OTUs")

# # env variables
# print("Environmental variables - all OTUs")
# # create a tiny anonymous function to include formula syntax in call
# abFrac0 <- mapply(function(x,data) rda(x~1, data), 
#                   commFull, met_list, SIMPLIFY=FALSE) # Reduced model
# 
# step.env <- mapply(function(x,y) ordiR2step(x, scope = formula(y)), 
#                    abFrac0, abFrac, SIMPLIFY=FALSE)
# 
# step.env # an rda model, with the final model predictor variables
# 
# print("Summary of environmental selection process - all OTUs")
# lapply(step.env, function(x) x$anova)
# print("ANOVA on full environmental selection - all OTUs")
# lapply(step.env, anova)
# 
# # save plot
# pdf(file = "./plots/JO2008_step_envM.pdf")
# # make plot
# lapply(step.env, plot)
# dev.off()

# spatial variables
print("Spatial variables - all OTU")
pcnm_df <- lapply(pcnm_list, function(x) as.data.frame(scores(x)))
bcFrac <- mapply(function(x,data) rda(x~., data), 
                 commFull, pcnm_df, SIMPLIFY=FALSE) # Full model
bcFrac0 <- mapply(function(x,data) rda(x~1, data), 
                  commFull, pcnm_df, SIMPLIFY=FALSE) # Reduced model
step.space <- mapply(function(x,y) ordiR2step(x, scope = formula(y)), 
                     bcFrac0, bcFrac, SIMPLIFY=FALSE)
step.space

# summary of selection process
print("Summary of spatial selection process - all OTU")
lapply(step.space, function(x) x$anova)
print("ANOVA on full spatial selection - all OTU")
lapply(step.space, anova)

# save plot
pdf(file = "./plots/JO2008_step_spaceM.pdf")
# make plot
lapply(step.space, plot)
dev.off()

print("Partition Bray-Curtis dissimilarities - all OTUs")
vdist <- lapply(commFull, vegdist)
pbcd <- mapply(function(x,y,z) varpart(x, ~., y, data = z),
               vdist, scores_list, met_list, SIMPLIFY=FALSE)
pbcd

#cleanup
# remove objects to be replaced/no longer needed
rm(vdist,pbcd, commFull)
rm(abFrac, aFrac,abFrac0, pcnm_df, bcFrac, bcFrac0, step.space) #step.env,

# I have removed the variation decomposition with parsimonious variables, 
# since it was frequently failing and would likely cause issues.
