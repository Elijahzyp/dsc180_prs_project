## Build and validate a PRS model


This part of the repo stored all of the file that are needed for the following steps. Please refer to [this notebook](PRS_Single_Gene_Pipeline/prs_final.ipynb) for the details of each step. 

1. **Prepare Data**: Download the 1000 Genomes Project PLINK files, gene expression data, gene annotation data, and GWAS summary statistics, ensuring all files are organized in the same directory.

2. **Format Individual Genotype Data**: Convert 23andMe genotype data into PLINK format for compatibility with downstream analysis.

3. **Define cis-regions**: Use the gene annotation file to define cis-regions as ±500 kb around the transcription start and end sites of each gene.

4. **Perform cis-eQTL Analysis**: Run linear regression models in R or Python to identify SNP-gene associations, adjusting for population structure and other covariates.

5. **Visualize Results**: Generate LocusZoom plots to visualize significant SNP-gene associations and linkage disequilibrium patterns.

6. **Prune SNPs for PRS**: Apply pruning and thresholding (P+T) using PLINK to select SNPs based on linkage disequilibrium and p-value thresholds.

7. **Calculate PRS**: Use PLINK to compute PRS by summing weighted SNP effect sizes for each individual.

8. **Validate PRS Models**: Split the dataset into training and testing sets, optimize p-value thresholds using cross-validation, and evaluate performance on the testing set.
