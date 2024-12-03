# dsc180_prs_project
A repo used for dsc180 capstone project.

## Introduction
This capstone project, led by Professor Tiffany Amariuta, investigates the genetic regulation of complex traits by integrating genotype and expression data. The project applies genome-wide association studies (GWAS), cis-expression quantitative trait loci (cis-eQTL) analyses, and polygenic risk score (PRS) modeling to understand how genetic variation influences regulatory mechanisms and trait variability.

For more details about the project, please refer to the [website](https://tiffanyamariuta.github.io/capstone-genetic-risk-prediction/).



## Environment
All project dependencies and environment configurations are specified in this [yml file](https://github.com/Elijahzyp/dsc180_prs_project/blob/main/environment.yml). Please create the environment using conda to ensure compatibility with the tools and libraries used in the analysis.

## Data
The datasets utilized in this project include genotypic, phenotypic, and regulatory information from multiple sources. Due to their large sizes, these datasets are hosted on Google Drive. Below is a detailed description of each dataset:

**1000 Genomes Project (PLINK Format)**

- Provides genotype data for population-level analysis.
- Files contain .bed, .bim, and .fam files, containing SNP information and sample metadata.

https://drive.google.com/file/d/1ncj0wm8ll0eOEunqNYuZxOosg8kr_fXx/view?usp=sharing

**Gene Annotation Data (gene_annot.txt)**

-  Defines genomic coordinates for protein-coding genes.

- Includes columns for gene ID, chromosome, start/end positions, and gene type.

https://drive.google.com/file/d/1JtAdensaa0iikOpRbUW0DPgXFyEsbNAb/view?usp=drive_link

**Gene Annotation Data (gene_annot.txt)**

- Defines genomic coordinates for protein-coding genes.

- Includes columns for gene ID, chromosome, start/end positions, and gene type.

https://drive.google.com/file/d/1qbDf-ZrSF3aMGM8djNhVcZuNn0gBUTvB/view?usp=drive_link

**Individual Genotype Data (23andMe)**

- Integrates personal-level genetic information for PRS modeling.

- Converted to PLINK format for compatibility with GWAS and PRS pipelines.

https://drive.google.com/file/d/1L2Zblx6j5SyLYoqYKt1Tb1WuP1Bkq4gc/view?usp=drive_link

**Gene Assembly**

- Provides assembled gene sequences and metadata for genomic analysis.

- Includes annotations necessary for aligning and interpreting gene regions.

https://drive.google.com/file/d/1KZ2TMokSjfCJ33Cgm7uYZYEDQkVMsS_H/view?usp=drive_link

**GWAS Summary Statistics**

- **Purpose**: Used to calculate PRS for height, asthma, and ankylosing spondylitis.

​	Height GWAS Summary Statistics

​	https://drive.google.com/file/d/16eWxZxwjfbeNY6qfmPCl-tNMYF3blYnr/view?usp=drive_link

​	Ankylosing Spondylitis GWAS Summary Statistics

​	https://drive.google.com/file/d/1mreVcKG-R8NIPr6BvSDZrv7BXvLEMV0k/view?usp=sharing

​	Child_Onset_Asthma GWAS Summary Statistics

​	https://drive.google.com/file/d/1KhrCmHg08br9UiIzaGB_CZlQJDtNoG1O/view?usp=sharing

**cis-eQTL Results (cis_eQTL_results_all_chromosomes.txt)**

- **Purpose**: Genome-wide cis-eQTL summary statistics for downstream analysis.

https://drive.google.com/file/d/1YXI-wSB-5ucDwJWUaqtyU24y9EsQH0at/view?usp=sharing

**Functional Annotation Data (ATAC-seq for B Cells)**

- **Purpose**: Identifies open chromatin regions in B cells for functional overlap analysis.

https://drive.google.com/file/d/1Qz8bX54PBC6dfY91MNkuipj4l0ZpLQ07/view?usp=sharing



## **Methods**

1. **Prepare Data**: Download the 1000 Genomes Project PLINK files, gene expression data, gene annotation data, and GWAS summary statistics, ensuring all files are organized in the same directory.

2. **Format Individual Genotype Data**: Convert 23andMe genotype data into PLINK format for compatibility with downstream analysis.

3. **Define cis-regions**: Use the gene annotation file to define cis-regions as ±500 kb around the transcription start and end sites of each gene.

4. **Perform cis-eQTL Analysis**: Run linear regression models in R or Python to identify SNP-gene associations, adjusting for population structure and other covariates.

5. **Visualize Results**: Generate LocusZoom plots to visualize significant SNP-gene associations and linkage disequilibrium patterns.

6. **Prune SNPs for PRS**: Apply pruning and thresholding (P+T) using PLINK to select SNPs based on linkage disequilibrium and p-value thresholds.

7. **Calculate PRS**: Use PLINK to compute PRS by summing weighted SNP effect sizes for each individual.

8. **Validate PRS Models**: Split the dataset into training and testing sets, optimize p-value thresholds using cross-validation, and evaluate performance on the testing set.

9. **Integrate 23andMe Data**: Incorporate individual genotype data into the PRS pipeline and compare results with population-level PRS distributions.
10. **Annotate Functional Variants**: Cross-reference cis-eQTL results with functional annotation datasets (e.g., ATAC-seq) to identify biologically relevant variants.

11. **Interpret Results**: Summarize findings through statistical analysis, visualizations, and population-level comparisons for traits like height, asthma, and ankylosing spondylitis.



## Commands to reproduce

Please refer to data_combine.sh to get to know how to combine individual genomic data with a group

Please refer to prep_prs to get to know how to convert a gwas data into correct format.

Please refer to PRS.sh to get to know how to perform a final PRS calcualtion.

Please refer to prs_graph.py to get to know how to make a graph based on your previous output

Readme. On Construction...
