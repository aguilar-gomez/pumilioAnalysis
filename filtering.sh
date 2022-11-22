#Generate vcf, doing calls was too time consuming and I used genotype likelihoods anyway
bcftools mpileup --threads 16 -f $REF -b pumilio.r.bamlist -I -q 10 -Q 20 --rf 2 –ff 1804 -Oz -o pumilio_likelihoods.vcf.gz &


#Angsd filtering
angsd -bam $bamlist -out pumilio_angsd -minInd 250 -setMinDepthInd 1 -minMapQ 25 -minQ 25 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -GL 1 -doMaf 1  -doMajorMinor 4 -doGlf 2 -SNP_pval 1e-6 -nThreads 10 -ref $REF -setMaxDepthInd 50 -doCounts 1 -skipTriallelic 1

#Getting sites
perl maf2pos.pl pumilio_angsd.mafs.gz> pumilio_angsd_rescaffold &
angsd sites index pumilio_angsd_rescaffold
cut -f1 pumilio_angsd_rescaffold |uniq >pumilio_angsd_rescaffold.rf

#snpCleaner time
bcftools index pumilio_likelihoods.vcf.gz
#Call and filter sites
bcftools view pumilio_likelihoods.vcf.gz --threads 16 -Oz -R pumilio_angsd_rescaffold| bcftools call -f GQ -O v -c >pumilio_called_nomaf_secondattempt.vcf &

#Run snpCleaner
# -u 1 minimum depth for individual to be covered, -k min number of covered individuals
# -H excess of heterozygous test
#-e end distance bias
# -b min pvalue for base quality bias, -S min pvalue strand bias, -f min pvalue map quality bias
#-h HWE
cat pumilio_called_nomaf_secondattempt.vcf | snpCleaner.pl -u 1 -k 200 -H 1e-6 -b 1e-100 -S 1e-4 -f 1e-4 -e 1e-4 -h 0 -B pumilio.snpCleaner.format -p pumilio_failedsites.txt.bz > pumilio_snpCleanv4.vcf
 
