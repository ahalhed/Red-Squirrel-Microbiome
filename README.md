# Red-Squirrel-Microbiome
Repository for scripts associated with the publication for the North American red squirrel microbiome project.

## What can you find in this directory?
The four character alphanumeric identifier in each filename (AG08, for example) corresponds to the sampling grid and collection year for the samples of interest in that script. There are 8 possible combinations (3 years, 6 grids). Files labelled with "core" were run with a community matrix containing only the core OTUs (present in >95% of collected samples and identified as significant contributors to microbial diversity using an occupancy-abundance model). Those labelled "nc" were run with a community matrix containing those OTUs NOT present in the core microbiome (aka non-core microbiome).

### Script files
There is a single QIIME2 script in this directory, labelled *rs-q2.sh*. This particular script was run in chunks as individual batch jobs, as many steps were quite long running. These individual jobs are stored in a separate private GitHub repository. The resulting QIIME2 visualizations and artifacts are stored in a directory on the GRAHAM cluster of SHARCNET (ComputeCanada).

The remaining files in the scripts folder are a combination of shell and R scripts; corresponding scripts have the same naming but different extensions. All of these files, with the exception *rs-q2.sh*, are post-QIIME2 analysis in R. 

### Plots
In addition to STDOUT, discussed below, most scripts output a variety of figures corresponding to different steps of the analysis. The labelling of these files is consistent with the type of plot, microbial community (full (no added label), core, or non-core (nc)), and grid/year combination.

### Output files
Files labelled with the ".out" extension are text files containing the STDOUT from all SHARCNET batch jobs. These files are named according to the grid/year combination whose spatial analysis STDOUT is contained within the file and for the SHARCNET job number.

### Data files
These are small data files generated as part of some analyses.
