#Convert to plink bed format
plink --vcf pumilio_snpCleanv4.vcf --make-bed --out pumilio.v4 --double-id --allow-extra-chr&


#Calculate relatedness matrix 
gemma -bfile pumilio.v4 -gk 2 -o pum.v4
 
#Use only the relatedness matrix, not PCs

#pumilio.v4.fam has sex as phenotype, and the plink column for sex is all set to 0
#Sex gwas remove covariate of sex!
nohup gemma -bfile pumilio.v4 -lmm 4 -o sex_nopcs_rescaffold -k pum.v4.sXX.txt 1> sex.out 2> out.sexgemma &  

#pumilio.v4.fam has sex all the phenotypes, and the sex column is filled out
#Columns
#["pop","bam","mom","dad","sex","class1","B1","S1U","S1V","S1B","S1G","S1Y","S1R","blackproportion","vgg16_k6","PC1","PC2","sexgwas"]]

for n in $(seq 1 13)
do
echo $n
nohup gemma -bfile pumilio.v4 -lmm 4 -o n${n}_rescaffold -k pum.v4.sXX.txt -n $n 1> n$n.out 2> out.gemma.n$n &  
done

