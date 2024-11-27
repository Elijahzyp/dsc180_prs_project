## After you finished date combining and gwas data preparation(prep_prs.py), you can do this process to run PRS



plink --bfile final_combined --score hight_data.gwas.txt 1 2 3 --out prs_for_hight --allow-extra-chr


## change the file after score flag to be your own gwas data, and bfile to be your own genomic data. 

