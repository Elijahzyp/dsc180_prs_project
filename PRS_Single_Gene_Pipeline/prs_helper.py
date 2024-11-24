

import pandas as pd
import numpy as np
from pandas_plink import read_plink1_bin
import matplotlib.pyplot as plt
from scipy import stats

hm3snps = pd.read_csv('data/w_hm3.snplist',delimiter='\t')


def cis_eQTL(window_size, gene_name, G_train, expression, gene_annotation):
    """
    Find cis-eQTLs for a specific gene within a given window size.
    
    Args:
        window_size (int): Window size in base pairs.
        gene_name (str): Name of the gene to analyze.
        G_train (pandas.DataFrame): Genotype data from PLINK files.
        expression (pandas.DataFrame): Expression data with samples as rows and genes as columns.
        gene_annotation (pandas.DataFrame): Annotation data with gene information (gene_name, start, end, chromosome).
    
    Returns:
        pandas.DataFrame: cis-eQTL results for the specified gene.
    """
    # Find the gene's start, end, and chromosome
    gene_info = gene_annotation[gene_annotation['SYM'] == gene_name]
    if gene_info.empty:
        raise ValueError(f"Gene {gene_name} not found in annotation data.")
    
    gene_start = gene_info.iloc[0]['START']
    gene_end = gene_info.iloc[0]['STOP']
    gene_chromosome = str(gene_info.iloc[0]['CHR'])
    
    
    # Define the window around the gene
    snp_start = gene_start - window_size
    snp_end = gene_end + window_size

    # Filter genotype data for SNPs within the window on the same chromosome
    snps_in_window = G_train.snp.loc[
        (G_train.snp['chrom'] == gene_chromosome) & 
        (G_train.snp['pos'] >= snp_start) & 
        (G_train.snp['pos'] <= snp_end)
    ]
    #return snps_in_window
    #if snps_in_window.empty:
    #    raise ValueError(f"No SNPs found within the window for gene {gene_name}.")
    
    snp = snps_in_window
    selected_expression = expression[expression['gene'] == 'ENSG00000177000'][np.array(G_train.fid)].T
    filtered_snps = G_train.where(G_train['snp'].isin(np.array(snp)), drop=True)
    #filtered_snps_bfile = filtered_snps.compute()
    filtered_snps_bfile = filtered_snps
    _columns = np.array(filtered_snps_bfile.snp)
    _index = np.array(filtered_snps_bfile.sample)
    snp_df = pd.DataFrame(np.array(filtered_snps_bfile), index = _index, columns = _columns)
    entire_df = pd.merge(selected_expression, snp_df, left_index=True, right_index=True, how='inner')
    
    
    snp_position = pd.DataFrame({
        'SNP': G_train['snp'],  # The SNP column
        'position': G_train['pos']  # The SNP positions (base pair positions)
    })

    results = []

    for i in range(len(snp)):

        y = entire_df.iloc[:,0]

        x = entire_df[np.array(snp)[i]] 
        if len(np.unique(x)) == 1:
            print(f"Skipping SNP at index {i}: All x values are identical.")
            continue


        each_states = stats.linregress(x, y)
        results.append([np.array(snp)[i], each_states.slope, each_states.pvalue])
    
    results_df = pd.DataFrame(results, columns=["SNP", "Beta", "P-value"])
    results_df = pd.merge(results_df, snp_position, on='SNP')
    
    return results_df





def plot_locuszoom(result_position_df, pval_col='P-value', position_col='position', figsize=(10, 6)):
    """
    Generate a LocusZoom-style plot for eQTL analysis.
    
    Args:
        result_position_df (pd.DataFrame): DataFrame containing genomic positions and P-values.
        pval_col (str): Column name for P-values in the DataFrame. Default is 'P-value'.
        position_col (str): Column name for genomic positions in the DataFrame. Default is 'position'.
        figsize (tuple): Figure size for the plot. Default is (10, 6).
    """
    # Calculate -log10(P-value)
    result_position_df['-log10(P-value)'] = -np.log10(result_position_df[pval_col])
    
    # Create the plot
    plt.figure(figsize=figsize)
    plt.scatter(result_position_df[position_col], result_position_df['-log10(P-value)'], 
                c='blue', alpha=0.6)
    
    # Add labels, title, and grid
    plt.xlabel('Genomic Coordinate (Position)', fontsize=12)
    plt.ylabel('-log10(P-value)', fontsize=12)
    plt.title('LocusZoom Plot for eQTL Analysis', fontsize=14)
    plt.grid(True)
    
    # Display the plot
    plt.show()
    
    
    