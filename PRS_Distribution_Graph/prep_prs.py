import os
import pandas as pd
main_dir = os.getcwd()[:-3]

## This function could help you to convert the gwas summary stats file into correct format. 
def summary_stats_prep(input_file_name, input_columns, output_file_name, separator):

  gwas_summary = pd.read_csv(input_file_name, sep=separator)
  prs_data = gwas_summary[input_columns].copy()
  prs_data.columns = ['SNP', 'ALLELE', 'BETA']
  prs_data['ALLELE'] = prs_data['ALLELE'].str.upper()
  prs_data = prs_data.drop_duplicates(subset='SNP', keep='first')
  prs_data = prs_data.dropna()
  prs_data.to_csv(output_file_name, sep=' ', index=False)
  
