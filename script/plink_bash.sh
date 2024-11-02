p_val=0.0001
r2=0.1
plink --bfile 1000G_eur_train_${gene} --clump-p1 $p_val --clump-r2 $r2 --clump-kb 250 --clump ${eqtl_sumstats_file} --clump-snp-field SNP --clump-field P --out 1000G_eur_train_${gene}
#extract SNPs from clumped output
awk 'NR!=1{print $3}' 1000G_eur_${gene}.clumped > PRS.SNPs.${gene}
plink --bfile 1000G_eur_train_${gene} --extract PRS.SNPs.${gene} --make-bed --out 1000G_eur_PRS_${gene}
#make score file in R
  for (gene in genes){
  y <- fread(paste0(“PRS.SNPs.“,gene), header = F)$V1
  w <- which(nchar(y) > 0)
  y <- y[w]
  m <- match(y, sumstats$SNP)
  score <- cbind(sumstats$SNP[m], sumstats$A1[m], sumstats$BETA[m]) #want to make sure BETA corresponds to A2 (not A1)
  write.table(score, file = paste0(gene, “_score_file.txt”), row.names = F, col.names = F, sep = “\t”, quote = F)
  }
#make plink files for test set + score people in the test set
plink --bfile 1000G_eur_PRS_test_${gene} --out PRS_test_${gene} --score ${gene}_score_file.txt 1 2 3