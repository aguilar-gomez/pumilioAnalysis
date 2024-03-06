
#Generate vcf with outgroup
cat pumilio.r.bamlist sylvatica.r.bamlist > pumiSyl.r.bamlist
bamlist=pumiSyl.r.bamlist
REF=/space/s2/diana/pumilio/rescaffolded_genome/Opum.rescaffold.fasta
angsd -bam $bamlist -out pumilio_outSyl -minInd 340 \
      -setMinDepthInd 1 -minMapQ 25 -minQ 25 -remove_bads 1 \
      -uniqueOnly 1 -only_proper_pairs 1  -GL 1 -doMaf 1 \
      -doMajorMinor 4 -doGlf 2 -SNP_pval 1e-4 -nThreads 20 \
      -ref $REF -setMaxDepthInd 50 -doCounts 1 -skipTriallelic 1 \
      -dosnpstat 1 -doHWE 1 -sb_pval 1e-4 -doGeno 3 \
      -hetbias_pval 1e-6 -sites pumilio_sites_all.pos \
      -rf pumilio_sites_all.rf -doBcf 1 -doPost 1 --ignore-RG 0 \
      -edge_pval 1e-4 -mapQ_pval 1e-4

#Convert to plink format
plink --vcf pumilio_outSyl.vcf --make-bed --out pumSyl_2024 \
      --geno 0.0 --allow-no-sex --allow-extra-chr --double-id \
      --set-missing-var-ids @:# 

#Remove find, to keep outgroup with one individual
plink --bfile pumSyl_2024 --geno 0.0 --allow-no-sex --allow-extra-chr\
      --double-id --maf 0.0057 --remove Osyl91 --out OsylOpum \
      --make-bed --set-missing-var-ids @:# 
      
#Generate frq file
plink --bfile OsylOpum --geno 0.0 --freq --missing --within pumiSyl.clust \
      --out OsylOpum2024 --double-id --allow-extra-chr

#Convert to treemix format
gzip OsylOpum2024.frq.strat
plink2treemix.py OsylOpum2024.frq.strat.gz OsylOpum2024.frq.gz


#Run treemix
FILE=OsylOpum2024
for i in {0..5}
do
 treemix -i $FILE.frq.gz -m $i -o $FILE.$i -bootstrap -k 500  -root Osylvatica> treemix_${i}_log &
done 
