setwd("/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/Red-Squirrel-Microbiome/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(vegan)
library(tidyverse)

theme_set(theme_bw())

nReads <- 4000                                                            # input dataset needs to be rarified and the rarifaction depth included 
otu <- read_qza("../filtered-table-10.qza")$data
map <- read_q2metadata("../input/RS_meta.tsv") # this is metadata

otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # mean relative abundance
occ_abun <- rownames_to_column(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance data frame

# Ranking OTUs based on their occupancy
# For caluclating raking index we included following conditions:
#   - time-specific occupancy (sumF) = frequency of detection within time point (genotype or site)
#   - replication consistency (sumG) = has occupancy of 1 in at least one time point (genotype or site) (1 if occupancy 1, else 0)

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
  gather(SampleID, abun, -otu) %>%
  left_join(map, by = 'SampleID') %>%
  group_by(otu, CollectionDate) %>%
  summarise(time_freq=sum(abun>0)/length(abun),            # frequency of detection between time points
            coreTime=ifelse(time_freq == 1, 1, 0)) %>%     # 1 only if occupancy 1 with specific time, 0 if not
  group_by(otu) %>%
  summarise(sumF=sum(time_freq),
            sumG=sum(coreTime),
            nS=length(CollectionDate)*2,           
            Index=(sumF+sumG)/nS)                 # calculating weighting Index based on number of time points detected and 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by='otu') %>%
  transmute(otu=otu,
            rank=Index) %>%
  arrange(desc(rank))

# Calculating the contribution of ranked OTUs to the BC similarity
BCaddition <- NULL

# calculating BC dissimilarity based on the 1st ranked OTU
otu_start=otu_ranked$otu[1]                   
start_matrix <- as.matrix(otu[otu_start,])
start_matrix <- t(start_matrix)
x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 
BCaddition <- rbind(BCaddition,df_s)
# calculating BC dissimilarity based on additon of ranked OTUs from 2nd to 500th. Can be set to the entire length of OTUs in the dataset, however it might take some time if more than 5000 OTUs are included.
for(i in 2:500){                              
  otu_add=otu_ranked$otu[i]                       
  add_matrix <- as.matrix(otu[otu_add,])
  add_matrix <- t(add_matrix)
  start_matrix <- rbind(start_matrix, add_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
}
# calculating the BC dissimilarity of the whole dataset (not needed if the second loop is already including all OTUs) 
x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))   
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names,x)
names(df_full)[2] <- length(rownames(otu))
BCfull <- left_join(BCaddition,df_full, by='x_names')

rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%            # mean Bray-Curtis dissimilarity
  arrange(desc(-MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))   # proportion of the dissimilarity explained by the n number of ranked OTUs
Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

#Creating thresholds for core inclusion 

#Method: 
#A) Elbow method (first order difference) (script modified from https://pommevilla.github.io/random/elbows.html)
fo_difference <- function(pos){
  left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
  right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
  return(left - right)
}
BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

#B) Final increase in BC similarity of equal or greater then 2% 
lastCall <- last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))

#Creating occupancy abundance plot
occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))]] <- 'core'
# add 95% occupancy threshold for core
occ_abun$Community <- ifelse(occ_abun$otu_occ >= 0.95 & occ_abun$fill == "core", "Confirmed Core",
                             ifelse(occ_abun$otu_occ < 0.95 & occ_abun$fill == "core", "Core Candidate",
                                    "Confirmed Non-core"))
# add a taxonomy column
tax <- read_qza("../taxonomy/SILVA-taxonomy-10.qza")$data %>%
  rename("otu" = "Feature.ID")
# clean up/separatee taxonomy labels
tax$Taxon <- tax$Taxon %>%
  str_replace_all("D_0__", "") %>%
  str_replace_all("D_1__", "") %>%
  str_replace_all("D_2__", "") %>%
  str_replace_all("D_3__", "") %>%
  str_replace_all("D_4__", "") %>%
  str_replace_all("D_5__", "") %>%
  str_replace_all("D_6__", "")
occ_abunT <- tax %>% 
  mutate("Kingdom" = word(.$Taxon, 1, sep = ";"), #k__
         "Phylum" = word(.$Taxon, 2, sep = ";"), #p__
         "Class" = word(.$Taxon, 3, sep = ";"), #c__
         "Order" = word(.$Taxon, 4, sep = ";"), #o__
         "Family" = word(.$Taxon, 5, sep = ";"), #f__
         "Genus" = word(.$Taxon, 6, sep = ";"), #g__
         "Species" = word(.$Taxon, 7, sep = ";")) %>% #s__
  # join to core labels
  right_join(occ_abun)

# exporting the data frame with which are core
# to load into manuscript figure file
write.table(occ_abunT, file = "./data/core.csv", sep = ",", quote = F, row.names = F)
