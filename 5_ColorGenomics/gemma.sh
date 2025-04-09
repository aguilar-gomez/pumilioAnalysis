#relatedness file
nano relatedIndv

#New version ran with plink2
plink2 --bcf pumilio_version11.bcf  --make-bed --out pumilio.v11 --double-id --allow-extra-chr&

#Add phenotypes
cat binaryglarek6_PCs.tsv |tail -n +2>pumilio.v11.fam

#Extract related individuals
while read line
do
grep "$line" pumilio.v11.fam
done < relatedIndv > related2Exclude.fam

#Remove samples
plink2 --bfile pumilio.v11 \
       --remove related2Exclude.fam \
       --make-bed \
       --out unrelated_pumilio --allow-extra-chr

#Removing individuals gets rid of the phenotypes
grep -v -f relatedIndv pumilio.v11.fam > unrelated_pumilioPhenotypesComplete.fam
#The following outputs nothing so I can replace the file
diff <(cut -f1,2 unrelated_pumilio.fam) <(cut -f1,2 unrelated_pumilioPhenotypesComplete.fam)
#Remove CLASS phenotype because it has missing values 
cut -f1-5,7- unrelated_pumilioPhenotypesComplete.fam >unrelated_pumilio.fam

#Calculate relatedness matrix
gemma -bfile unrelated_pumilio -gk 2 -o pum_unrelated
#GEMMA 0.98.5 (2021-08-25) by Xiang Zhou, Pjotr Prins and team (C) 2012-2021
#Reading Files ... 
## number of total individuals = 321
## number of analyzed individuals = 321
## number of covariates = 1
## number of phenotypes = 1
## number of total SNPs/var        = 10231695
## number of analyzed SNPs         =  1887184
#Calculating Relatedness Matrix ... 
#================================================== 100%
#**** INFO: Done.

#pumilio.v4.fam has sex all the phenotypes, and the sex column is filled out
#Columns
mv output/pum_unrelated.sXX.txt .
#["pop","bam","mom","dad","sex","B1","S1U","S1V","S1B","S1G","S1Y","S1R","blackproportion","vgg16_k6","PC1","PC2","sexgwas"]]
for n in $(seq 1 12)
do
echo $n
nohup gemma -bfile unrelated_pumilio -lmm 4 -o n${n}_unrelated -k pum_unrelated.sXX.txt -n $n 1> n$n.out 2> out.gemma.n$n &  
done

nohup sh runGemma.sh > outGemma &

#Change gemma default setting to lose less SNPs for asso
#-miss .1 (allow missing up to 10%, default 5%)

for n in $(seq 4 7)
do
echo $n
nohup gemma -bfile unrelated_pumilio -lmm 4 -o n${n}_filters -k pum_unrelated.sXX.txt -n $n -miss .1 -r2 0.999999 -maf 0.005 1> n$n.out 2> out.gemmaFilter.n$n &  
done

nohup sh runGemmaFilter.sh > outGemmaFilter &
