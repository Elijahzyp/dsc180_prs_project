import os
import pandas as pd
main_dir = os.getcwd()[:-3]

## This function could help you to convert the gwas summary stats file into correct format. 
def prepare_summary_stats_for_prs(input_file_name, input_columns, output_file_name, separator):

  gwas_summary = pd.read_csv(main_dir + "PRS/data/" + input_file_name, sep=separator)
  prs_data = gwas_summary[input_columns].copy()
  prs_data.columns = ['SNP', 'ALLELE', 'BETA']
  prs_data['ALLELE'] = prs_data['ALLELE'].str.upper()
  prs_data.to_csv(main_dir + "PRS/" + output_file_name, sep=' ', index=False)
  
