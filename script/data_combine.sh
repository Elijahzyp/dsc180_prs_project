## Command to combine one people's genome data into the data set that includes a group of people. 

## Please change the Tiffany into your own name

plink --bfile new_merge --bmerge Tiffany_sample.bed Tiffany_sample.bim Tiffany_sample.fam --make-bed --out final_combined --allow-extra-chr

## after running above command, you might see an error message, no worry , in the same directory, you can also see a 
## file called combined_with_Tiffany-merge.missnp
## run the following command

./plink --bfile 1000G_combined --exclude combined_with_Tiffany-merge.missnp --make-bed --out new_merge

## Then rerun the very first command
plink --bfile new_merge --bmerge Tiffany_sample.bed Tiffany_sample.bim Tiffany_sample.fam --make-bed --out final_combined --allow-extra-chr

## you should be good!