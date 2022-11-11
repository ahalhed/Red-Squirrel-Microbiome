#---
#title: "Squirrel Manuscript Figures"
#author: "Alicia Halhed"
#date: "10/04/2021"
#output: html_document
#---
# set working directory to output the plots into
# run on graham cluster interactively
# salloc --time=0-00:30:00 --mem=8G --account=def-cottenie
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
setwd("/home/ahalhed/MSc/Red-Squirrel-Microbiome/")
# attaching required packages for full analysis
# qiime2R to create phyloseq object
library(phyloseq)
library(zCompositions)
#devtools::install_github('ggloor/CoDaSeq/CoDaSeq')
library(CoDaSeq)
library(lubridate)
library(ggpubr)
library(qiime2R)
library(vegan)
library(tidyverse)
# set theme for ggplots
theme_set(theme_bw())

# calculate euclidean dissimilarity from an community matrix formatted ASV table
# relies on tidyverse, vegan
dis <- function(OTU, met) {
  # ASV is the community matrix containing the ASVs
  # met is a data frame containing the environmental metadata
  m <- rownames_to_column(met, var = "SampleID")
  df <- OTU %>% data.matrix(rownames.force = T) %>%
    vegdist(method = "euclidean") %>% 
    as.matrix %>% as.data.frame %>%
    rownames_to_column("sampleid1") %>%
    pivot_longer(-sampleid1, names_to = "sampleid2", values_to = "EucDis") %>%
    unique() %>% left_join(., m, 
                           by = c("sampleid1" = "SampleID")) %>%
    # .x for sample 1, .y for sample 2
    left_join(., m, by = c("sampleid2" = "SampleID"))
  return(df)
}


gr <- function(df){
  # df is a data frame containing metadata and Bray-Curtis distances
  # create a new grid column
  df$Grid <- ifelse(df$Grid.x == df$Grid.y, "Same Grid", "Different Grid")
  # add location information
  df$LocationIndividual <- ifelse(df$sampleid1 == df$sampleid2, "Same Sample",
                        ifelse(df$`Grid.x` != df$`Grid.y`, "Different Grid",
                               ifelse(df$`Squirrel.ID.x` == df$`Squirrel.ID.y` &
                                        df$`Location.X.x` == df$`Location.X.y` &
                                        df$`Location.Y.x` == df$`Location.Y.y`,
                                      "Same Squirrel, Same Location",
                                      ifelse(df$Squirrel.ID.x != df$Squirrel.ID.y &
                                               df$`Location.X.x` == df$`Location.X.y` &
                                               df$`Location.Y.x` == df$`Location.Y.y`,
                                             "Different Squirrel, Same Location",
                                             ifelse(df$Squirrel.ID.x != df$Squirrel.ID.y &
                                                      df$`Location.X.x` != df$`Location.X.y` |
                                                      df$`Location.Y.x` != df$`Location.Y.y`,
                                                    "Different Squirrel, Different Location",
                                                    ifelse(df$Squirrel.ID.x == df$Squirrel.ID.y &
                                                             df$`Location.X.x` != df$`Location.X.y` |
                                                             df$`Location.Y.x` != df$`Location.Y.y`,
                                                           "Same Squirrel, Different Location", NA))))))
  # add date information
  df$CollectionDate <- ifelse(df$CollectionDate.x == df$CollectionDate.y, "Same Date", "Different Date")
  return(df)
}

# add interval column
INT <- function(df) {
  # df is a dataframe containing jaccard distances and specific comparison information
  # first find the date range with interval
  df1 <- interval(df$`CollectionDate.x`, df$`CollectionDate.y`)
  # find length of the interval
  df2 <- as.period(df1, unit = "days")
  # extract number of days only
  df3 <- gsub("([0-9]+).*$", "\\1", df2)
  df4 <- df3 %>% as.numeric %>% abs
  return(df4)
}

XY_month <- function(metadata, grid, year, month) {
  m <- rownames_to_column(metadata, var = "SampleID")
  df1 <- subset(m, Grid == grid, 
                select = c("SampleID", "Location.X", "Location.Y", "Year", "Month"))
  df2 <- subset(df1, Year == year, 
                select = c("SampleID", "Location.X", "Location.Y", "Month"))
  df3 <- subset(df2, Month == month, 
                select = c("SampleID", "Location.X", "Location.Y"))
  df4 <- column_to_rownames(remove_rownames(df3), var = "SampleID")
  return(df4)
}

## get the data
print("Read in the Data")
print("Building phyloseq object")
ps <- qza_to_phyloseq(features = "../Red-Squirrel-MicrobiomeDATA/ASV-table-10-filtered.qza",
                      metadata = "../Red-Squirrel-MicrobiomeDATA/RS_meta.tsv") %>%
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), sample_data(.))

# based on the meta function from the microbiome package
print("Read in the metadata")
meta <- as(sample_data(ps), "data.frame")
rownames(meta) <- sample_names(ps)

print("Full ASV table")
print("Aitchison transformation")
# rows are ASVs, then transposed to ASVs as column
# impute the ASV table
OTUimp <- cmultRepl(otu_table(ps), label=0, method="CZM", output="p-counts") # all ASVs
# compute the aitchison values
OTU_full <- codaSeq.clr(abs(OTUimp)) %>% as.data.frame

## Core/non-core divide
print("Finding core microbiome")
print("Extract 75% Occupancy from BC Similarity Core")
# read in occupancy/abundance information
occ_abun <- read.csv("../Red-Squirrel-MicrobiomeDATA/core.csv")
# new column for just core and non-core
occ_abun$plot <- ifelse(occ_abun$Community == "Confirmed Core", "Core", "Non-core")
# get the ASVs identified as core contributors to beta diversity
# and greater than 75% occupancy (confirmed core)
cOTU <- occ_abun[which(occ_abun$Community == "Confirmed Core"),]
# make the new data frames
print("Subset the ASV table to find core and non-core ASVs")
OTU_core <- select(OTU_full, one_of(cOTU$otu))
OTU_nc <- select(OTU_full, -one_of(cOTU$otu))

## Figure 1 - Core Community Cutoff
# This plot shows the fraction of the ASVs included in the core microbiome
# create the framework for the plot
fig1 <- ggplot(occ_abun, aes(y = otu_occ, x = otu_rel, shape = plot, color = plot)) + #, color = plot
  geom_point() +
  # log transform the x axis
  scale_x_log10() +
  # add 75% threshold
  annotate("text", x = 0.00001, y = 0.78, label = ">75% Occupancy") +
  geom_hline(yintercept = 0.75, linetype = "dashed", size = 0.5) +
  scale_color_viridis_d() +
  # add axis labels
  labs(x = "Mean Relative Abundance of Each ASV (log10)", y = "Occupancy (Proportion of Samples)",
       color = "Community", shape = "Community")

# export plot 1 to a file
tiff("/home/ahalhed/MSc/Red-Squirrel-Microbiome/plots/figure1.tiff", width = 259, height = 169, units = 'mm', res = 400)
fig1 + theme(text = element_text(size = 20))
dev.off()

## Figure 2 - spatial pattern example
## KL 2008, June
# get XY data
XY_KL <- XY_month(meta, "KL", 2008, 6)
# euclidean distances
d_KL <- dist(XY_KL)
# PCNM
KL <- pcnm(d_KL)

# generate figure 2
tiff("/home/ahalhed/MSc/Red-Squirrel-Microbiome/plots/figure2.tiff", width = 240, height = 80, units = 'mm', res = 400)
par(mfrow=c(1,4))
# core (red)
ordisurf(XY_KL, scores(KL, choi=14), bubble = 4, col = "red", main = "PCNM 14")
# non-core (black)
ordisurf(XY_KL, scores(KL, choi=8), bubble = 4, col = "black", main = "PCNM 8")
ordisurf(XY_KL, scores(KL, choi=7), bubble = 4, col = "black", main = "PCNM 7")
ordisurf(XY_KL, scores(KL, choi=2), bubble = 4, col = "black", main = "PCNM 2")
dev.off()

## Figure 3 - Adjusted R2

# read in the data
# values taken from output files from monthly PCNM analyses
# pivot the data of interest into long format
adj <- read_csv("/home/ahalhed/MSc/Red-Squirrel-Microbiome/data/Figure3Data.csv")
adj <- adj[which(adj$Community != "Full"),]

# create plot for all adjusted R2 points
fig3 <- ggplot(adj, aes(Month, as.numeric(R2Adj), color = Community)) +
  geom_smooth(method = "lm") + #, aes(linetype = Community)
  geom_jitter(aes(shape = as.character(Year))) +
  scale_color_viridis_d() +
  #scale_color_manual(values=c("grey20", "black")) + 
  labs(y = expression(paste("Adjusted R"^"2")), shape = "Collection Year")

# exporting figure 3
tiff("/home/ahalhed/MSc/Red-Squirrel-Microbiome/plots/figure3.tiff", width = 240, height = 120, units = 'mm', res = 400)
fig3 + theme(text = element_text(size = 20))
dev.off()

print("is there a significant difference in the R2adj values based on the month and community of origin?")
print("All Adjusted R-squared Values - Spatial only")
adj %>% lm(R2Adj ~ Community*Month, data = .) %>% summary
print("All Adjusted R-squared Values - Core only")
adj[which(adj$Community=="Core"),] %>% lm(R2Adj ~ Month, data = .) %>% summary
print("All Adjusted R-squared Values - Non-core only")
adj[which(adj$Community=="Non-core"),] %>% lm(R2Adj ~ Month, data = .) %>% summary

## Figure 4 - LOESS regression
# calculate Aitchison dissimilarity (euclidean distance on CLR transformed ASV table)
core_dis <- dis(OTU_core, meta)
write.table(core_dis, file='./data/core-dis.tsv', quote=FALSE, sep='\t', row.names = F)
#core_dis <- read.table("./data/core-dis.tsv.gz", sep = "\t", header = T)
nc_dis <- dis(OTU_nc, meta)
write.table(nc_dis, file='./data/nonCore-dis.tsv', quote=FALSE, sep='\t', row.names = F)
#nc_dis <- read.table("./data/nonCore-dis.tsv.gz", sep = "\t", header = T)
full_dis <- dis(OTU_full, meta)
write.table(full_dis, file='./data/full-dis.tsv', quote=FALSE, sep='\t', row.names = F)
#full_dis <- read.table("./data/full-dis.tsv.gz", sep = "\t", header = T)

# mutate the distance df to include date columns
core_dis <- core_dis %>%
  mutate(CollectionDate.x = as_date(.$CollectionDate.x),
         CollectionDate.y = as_date(.$CollectionDate.y))
nc_dis <- nc_dis %>%
  mutate(CollectionDate.x = as_date(.$CollectionDate.x),
         CollectionDate.y = as_date(.$CollectionDate.y))
full_dis <- full_dis %>%
  mutate(CollectionDate.x = as_date(.$CollectionDate.x),
         CollectionDate.y = as_date(.$CollectionDate.y))
# add the interval column between collection dates and locations
# core
core_dis <- gr(core_dis)
core_dis$int <- INT(core_dis)
core_dis <- core_dis %>%
  mutate(Individual = word(core_dis$LocationIndividual, 1, sep = ","),
         Location = word(core_dis$LocationIndividual, 2, sep = ", "))
# non-core
nc_dis <- gr(nc_dis)
nc_dis$int <- INT(nc_dis)
nc_dis <- nc_dis %>%
  mutate(Individual = word(nc_dis$LocationIndividual, 1, sep = ","),
         Location = word(nc_dis$LocationIndividual, 2, sep = ", "))
# full
full_dis <- gr(full_dis)
full_dis$int <- INT(full_dis)
full_dis <- full_dis %>%
  mutate(Individual = word(full_dis$LocationIndividual, 1, sep = ","),
         Location = word(full_dis$LocationIndividual, 2, sep = ", "))

# break up different groups
# core
core_dsSL <- filter(core_dis, LocationIndividual %in% "Different Squirrel, Same Location")
core_dsDL <- filter(core_dis, LocationIndividual %in% "Different Squirrel, Different Location")
core_ssSL <- filter(core_dis, LocationIndividual %in% "Same Squirrel, Same Location")
core_ssDL <- filter(core_dis, LocationIndividual %in% "Same Squirrel, Different Location")
# non-core
nc_dsSL <- filter(nc_dis, LocationIndividual %in% "Different Squirrel, Same Location")
nc_dsDL <- filter(nc_dis, LocationIndividual %in% "Different Squirrel, Different Location")
nc_ssSL <- filter(nc_dis, LocationIndividual %in% "Same Squirrel, Same Location")
nc_ssDL <- filter(nc_dis, LocationIndividual %in% "Same Squirrel, Different Location")
# full
full_dsSL <- filter(full_dis, LocationIndividual %in% "Different Squirrel, Same Location")
full_dsDL <- filter(full_dis, LocationIndividual %in% "Different Squirrel, Different Location")
full_ssSL <- filter(full_dis, LocationIndividual %in% "Same Squirrel, Same Location")
full_ssDL <- filter(full_dis, LocationIndividual %in% "Same Squirrel, Different Location")

# bind the four together (10% of different location, different location)
linesC <- rbind(core_ssDL, core_dsSL) %>%
  rbind(., core_ssSL) %>% rbind(., core_dsDL %>% sample_frac(.1))
linesR <- rbind(nc_ssDL, nc_dsSL) %>%
  rbind(., nc_ssSL) %>% rbind(., nc_dsDL %>% sample_frac(.1))
linesF <- rbind(full_ssDL, full_dsSL) %>%
  rbind(., full_ssSL) %>% rbind(., full_dsDL %>% sample_frac(.1))
# within each of the years
linesYR <- linesR[which(linesR$Year.x == linesR$Year.y),]
linesYC <- linesC[which(linesC$Year.x == linesC$Year.y),]
linesYF <- linesF[which(linesF$Year.x == linesF$Year.y),]
# remove objects no longer needed
rm(linesC, linesR, linesF,
   core_dsSL, core_dsDL, core_ssSL, core_ssDL,
   nc_ssDL, nc_ssSL, nc_dsDL, nc_dsSL,
   full_ssDL, full_ssSL, full_dsDL, full_dsSL)

# non-core by year
# changing the order of facets/lines
linesYR$Individual_f = factor(linesYR$Individual, levels=c('Same Squirrel','Different Squirrel'))
linesYR$Location_f = factor(linesYR$Location, levels=c('Same Location','Different Location'))
# generate plot
ncYP <- ggplot(linesYR, aes(x = int, y = EucDis, linetype = Location_f)) +
  geom_smooth(method='loess', formula= y~x, color="black") + facet_grid(~ Individual_f) +
  labs(x = "Days between Sample Collection", y = "Aitchison Distance", 
       linetype = "Sampling Location") + theme(text = element_text(size = 20))

# core by year
# changing the order of facets/lines
linesYC$Individual_f = factor(linesYC$Individual, levels=c('Same Squirrel','Different Squirrel'))
linesYC$Location_f = factor(linesYC$Location, levels=c('Same Location','Different Location'))
# generate plot
coreYP <- ggplot(linesYC, aes(x = int, y = EucDis, linetype = Location_f)) +
  geom_smooth(method='loess', formula= y~x, color="black") + facet_grid(~ Individual_f) +
  labs(x = "Days between Sample Collection", y = "Aitchison Distance", 
       linetype = "Sampling Location") + theme(text = element_text(size = 20))

# full by year
# changing the order of facets/lines
linesYF$Individual_f = factor(linesYF$Individual, levels=c('Same Squirrel','Different Squirrel'))
linesYF$Location_f = factor(linesYF$Location, levels=c('Same Location','Different Location'))
# generate plot
fullYP <- ggplot(linesYF, aes(x = int, y = EucDis, linetype = Location_f)) +
  geom_smooth(method='loess', formula= y~x, color="black") + facet_grid(~ Individual_f) +
  labs(x = "Days between Sample Collection", y = "Aitchison Distance", 
       linetype = "Sampling Location") + 
  ggtitle("Full Microbial Community") + theme(text = element_text(size = 20))

# export figure 4
tiff("/home/ahalhed/MSc/Red-Squirrel-Microbiome/plots/figure4.tiff", width = 240, height = 240, units = 'mm', res = 400)
ggarrange(coreYP, ncYP, labels = c("A", "B"),
          nrow=2, common.legend = T)
dev.off()

# putting full in a supplemental figure
tiff("/home/ahalhed/MSc/Red-Squirrel-Microbiome/plots/supp4full.tiff", width = 110, height = 80, units = 'mm', res = 400)
fullYP
dev.off()