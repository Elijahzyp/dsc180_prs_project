import numpy as np
import pandas as pd
from pandas_plink import read_plink1_bin
from sklearn.linear_model import LinearRegression
from scipy import stats



def eQtl(gene_index, expression1, gene1, chr = 1):
    each_gene_expression = expression1[expression1['ID'] == gene1.iloc[gene_index]['SYM']]
    each_gene_expression_data = each_gene_expression.T.iloc[4:-1,]
    
    if each_gene_expression.shape[0] != 0:
    
        each_gene_expression_data.columns = [gene1.iloc[gene_index]['SYM']]

        each_gene_name = gene1.iloc[gene_index]['ID']
        each_G = read_plink1_bin(f"cis_snp/1000G.EUR.1.{each_gene_name}.bed", f"cis_snp/1000G.EUR.1.{each_gene_name}.bim", f"cis_snp/1000G.EUR.1.{each_gene_name}.fam", verbose=False)
        each_genotype_matrix = each_G.compute()
        _columns = np.array(each_G.snp)
        _index = np.array(each_G.sample)
        each_snp_df = pd.DataFrame(np.array(each_genotype_matrix), index = _index, columns = _columns)
        each_entire_df = pd.merge(each_gene_expression_data, each_snp_df, left_index=True, right_index=True, how='inner')
        each_gene_expression_data = each_entire_df.iloc[:,0]
        each_snp_df = each_entire_df.iloc[:,1:]
        results = []
        for snp_ in range(len(each_snp_df.columns)):
            
            the_expression = list(np.array(each_gene_expression_data.values))
            snp_genotype = list(np.array(each_snp_df.iloc[:,snp_].values))
            snp_genotype_reshaped = np.array(snp_genotype).reshape(-1, 1)
            model = LinearRegression()
            model.fit(snp_genotype_reshaped, the_expression)
            stats.linregress(snp_genotype, the_expression)
            each_states = stats.linregress(snp_genotype, the_expression)
            results.append([each_snp_df.columns[snp_], model.coef_[0], each_states.pvalue])


        results_df = pd.DataFrame(results, columns=["SNP", "Beta", "P-value"])
        snp_positions = pd.DataFrame({
        'snp': each_G['snp'],  # The SNP column
        'position': each_G['pos']  # The SNP positions (base pair positions)
        })
        result_position_df = pd.merge(snp_positions, results_df, left_on="snp", right_on='SNP')
        result_position_df['-log10(P-value)'] = -np.log10(result_position_df['P-value'])
        result_position_df['cis_gene'] = gene1.iloc[gene_index]['ID']
        return result_position_df
    else:
        return None