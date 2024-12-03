import pandas as pd
import pybedtools

# Load cis-eQTL data
cis_eqtl_df = pd.read_csv('cis_eQTL_results_all_chromosomes.txt', sep='\t')

# Convert cis-eQTL data to BED format (chromosome, start, end, SNP ID)
# Use SNP position as genomic location, with end = SNP_Pos + 1
cis_eqtl_bed_df = cis_eqtl_df[['Chr', 'SNP_Pos', 'SNP_Pos', 'SNP_ID']].copy()
cis_eqtl_bed_df.columns = ['chrom', 'start', 'end', 'name']  # Rename columns to match BED format
cis_eqtl_bed_df['end'] = cis_eqtl_bed_df['end'] + 1  # Set the end position as SNP position + 1

# If the ATAC-seq data uses a "chr" prefix, add the "chr" prefix to cis-eQTL data
cis_eqtl_bed_df['chrom'] = 'chr' + cis_eqtl_bed_df['chrom'].astype(str)

# Convert DataFrame to BedTool object
cis_eqtl_bed = pybedtools.BedTool.from_dataframe(cis_eqtl_bed_df)

# Load ATAC-seq data
atac_bed = pybedtools.BedTool('ENCFF421XIL.bed')

# Perform overlap analysis using bedtools
# `intersect` will identify all SNPs that overlap with open chromatin regions
overlap = cis_eqtl_bed.intersect(atac_bed, wa=True, wb=True)

# Check if the overlap file is empty
if overlap.count() == 0:
    print("No overlapping regions found. The output file is empty.")
else:
    # Load the overlap result into a DataFrame
    overlap_df = pd.read_csv(overlap.fn, sep='\t', header=None)
    print("The number of columns in the overlap file:", overlap_df.shape[1])
    print(overlap_df.head())  # Display the first few rows to verify the format

    # Set column names based on the number of columns
    if overlap_df.shape[1] == 13:  # If the file has 13 columns
        overlap_df.columns = ['Chr_eqtl', 'Start_eqtl', 'End_eqtl', 'SNP_ID', 
                              'Chr_atac', 'Start_atac', 'End_atac', 'Score', 
                              'Strand', 'SignalValue', 'pValue', 'qValue', 'Peak']
    elif overlap_df.shape[1] == 14:  # If the file has 14 columns
        overlap_df.columns = ['Chr_eqtl', 'Start_eqtl', 'End_eqtl', 'SNP_ID', 
                              'Chr_atac', 'Start_atac', 'End_atac', 'Score', 
                              'Strand', 'SignalValue', 'pValue', 'qValue', 'Peak', 'Extra']

    # Display results
    print("Number of overlapping SNPs:", overlap_df.shape[0])
    print(overlap_df.head())

    # Save the results to a file
    overlap_df.to_csv('overlapping_eqtl_atac_peaks.txt', sep='\t', index=False)
    print("The overlap analysis results have been saved to 'overlapping_eqtl_atac_peaks.txt'.")