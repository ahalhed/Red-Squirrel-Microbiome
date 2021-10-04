# Author: Alicia Halhed
# Species: North American Red Squirrel
# Sample: Fecal samples
# 16S rRNA
# working with qiime2-2019.10

# script starts here
# ________________________________________

# Working in the below one drive folder
cd /Users/aliciahalhed/OneDrive\ -\ University\ of\ Guelph/Alicia\'s\ Thesis/red-squirrel-data/rs-QIIME2-AH 
# /home/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome for SHARCNET directory
# Data in /Guelph/red-squirrel-data/original
# Add '&' to the end to put the process in the background
# To move a job currently running in the forground to the background, pause it with CTRL-Z and then run bg (moves last paused job to the background)
# Use disown to allow terminal to close 

```{r}
library(Biostrings)
library(dplyr)
# read the fasta file into R (this is the unzipped version of the zip archive provided by the authors)
rs_fasta <- readDNAStringSet("/Users/aliciahalhed/OneDrive\ -\ University\ of\ Guelph/Alicia\'s\ Thesis/red-squirrel-data/rs-QIIME2-AH/RS_seqs.fasta")
seq_name <- names(rs_fasta)
sequence <- paste(rs_fasta)
# put the fasta file into a data frame
rs_fasta_df <- data.frame(seq_name, sequence)
# write the full data frame out to a file for QIIME2
# 'SampleData[Sequences]' format in Q2 doesn't seem to like the sequences being "chunked" but would take the continuous (what I did locally)
rs_fasta_df %>%
  mutate(seq_name = paste(">", rs_fasta_df$seq_name, sep="")) %>% 
  write.table("RS_seqsR.fasta", sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

# Activate QIIME2
conda activate qiime2-2019.10

# Convert FASTA file to a QIIME2 artifact (need to use an unzip version of the fasta file)
#import sequences for OTU picking
qiime tools import \
  --input-path ./input/RS_seqsR.fasta \
  --output-path seqs.qza \
  --type 'SampleData[Sequences]'

# dereplicate the sequences
qiime vsearch dereplicate-sequences \
  --i-sequences seqs.qza \
  --o-dereplicated-table table.qza \
  --o-dereplicated-sequences rep-seqs.qza 

qiime metadata tabulate \
  --m-input-file table.qza \
  --o-visualization table.qzv

# tabulate the metadata
qiime metadata tabulate \
  --m-input-file ./input/RS_meta.tsv \
  --o-visualization tabulated-metadata.qzv
# to look at a visualization
qiime tools view tabulated-metadata.qzv

# de novo OTU clustering at 99% identity
qiime vsearch cluster-features-de-novo \
  --i-table table.qza \
  --i-sequences rep-seqs.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table OTU-table-dn-99.qza \
  --o-clustered-sequences OTU-rep-seqs-dn-99.qza

# Generation of Tree for Phylogenetic diversity analysis

# de novo tree
# threads to 0 to use all available cores
qiime alignment mafft \
  --p-parttree \
  --p-n-threads 0 \
  --i-sequences OTU-rep-seqs-dn-99.qza \
  --o-alignment aligned_sequences.qza 
# https://docs.qiime2.org/2019.10/tutorials/phylogeny/#reducing-alignment-ambiguity-masking-and-reference-alignments
# below took ~29.5 hours to run (8 cores, 16G RAM)
qiime alignment mask \
  --i-alignment aligned_sequences.qza \
  --o-masked-alignment masked_sequences.qza

qiime phylogeny fasttree \
  --i-alignment aligned_sequences.qza \
  --o-tree ./trees/unrooted_tree.qza 

qiime phylogeny midpoint-root \
  --i-tree ./trees/unrooted_tree.qza \
  --o-rooted-tree ./trees/rooted_tree.qza 


# so apparently not all the samples have metadata... see ./output/server5-jan16-26377006.out
# I will therefore filter the table and then run the rest of the script
qiime feature-table filter-samples \
  --i-table OTU-table-dn-99.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-filtered-table filtered-table.qza
# 541922 rows in this table

# June 26, 2020 - dropping features in less than 10 samples
qiime feature-table filter-features \
  --i-table filtered-table.qza \
  --p-min-samples 10 \
  --o-filtered-table filtered-table-10.qza
# 23425 rows in this table
# labelled as FilteredOTUtable.qza on GitHub repo

# filter representative sequences for this group
qiime feature-table filter-seqs \
  --i-data rep-seqs.qza \
  --i-table filtered-table-10.qza \
  --o-filtered-data rep-seqs-10.qza

# ASSIGN TAXONOMY
# Obtaining SILVA reference database (much larger database, will likely do a better job at classifying)
wget -O "silva-132-99-nb-classifier.qza" "https://data.qiime2.org/2019.10/common/silva-132-99-nb-classifier.qza"
# just doing this for OTUs that occur in ten or more samples (was taking forever to run otherwise)
# Classifying taxonomies
qiime feature-classifier classify-sklearn \
  --i-classifier ./references/silva-132-99-nb-classifier.qza \
  --i-reads rep-seqs-10.qza \
  --o-classification ./taxonomy/SILVA-taxonomy-10.qza
# Generating taxonomy visualization
qiime taxa barplot \
  --i-table filtered-table-10.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy-10.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./taxonomy/SILVA-dn-taxa-bar-plots-10.qzv
# Extracting Taxonomic Clasification
# Phylum
qiime taxa collapse \
  --i-table filtered-table-10.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy-10.qza \
  --p-level 2 \
  --o-collapsed-table ./taxonomy/SILVA-table-10-l2.qza
# Class
qiime taxa collapse \
  --i-table filtered-table-10.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy-10.qza \
  --p-level 3 \
  --o-collapsed-table ./taxonomy/SILVA-table-10-l3.qza
# Order
qiime taxa collapse \
  --i-table filtered-table-10.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy-10.qza \
  --p-level 4 \
  --o-collapsed-table ./taxonomy/SILVA-table-10-l4.qza
# Family
qiime taxa collapse \
  --i-table filtered-table-10.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy-10.qza \
  --p-level 5 \
  --o-collapsed-table ./taxonomy/SILVA-table-10-l5.qza
# Genus
qiime taxa collapse \
  --i-table filtered-table-10.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy-10.qza \
  --p-level 6 \
  --o-collapsed-table ./taxonomy/SILVA-table-10-l6.qza
# Species
qiime taxa collapse \
  --i-table filtered-table-10.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy-10.qza \
  --p-level 7 \
  --o-collapsed-table ./taxonomy/SILVA-table-10-l7.qza

# Close QIIME2
conda deactivate
