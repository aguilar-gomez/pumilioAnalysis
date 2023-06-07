#No SNP filtering, keep invariable
bamlist=pumilio.r.bamlist
REF=/space/s2/diana/pumilio/rescaffolded_genome/Opum.rescaffold.fasta
/home/diana/bin/angsd -bam $bamlist -out pumilio_sites_all -minInd 250 -setMinDepthInd 1 -minMapQ 25 -minQ 25 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -GL 1 -doMaf 1 -doMajorMinor 4 -doGlf 2 -nThreads 10 -ref $REF -setMaxDepthInd 50 -doCounts 1 -skipTriallelic 1

#SNP filtering
#Do all filters and generate beagle simultaneously
/home/diana/bin/angsd -bam $bamlist -out pumilio_version11 -minInd 250 -setMinDepthInd 1 -minMapQ 25 -minQ 25 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -GL 1 -doMaf 1 -doMajorMinor 4 -doGlf 2 -SNP_pval 1e-6 -nThreads 20 -ref $REF -setMaxDepthInd 50 -doCounts 1 -skipTriallelic 1 -dosnpstat 1 -doHWE 1 -sb_pval 1e-4 -hetbias_pval 1e-6 -doGeno 3 -sites pumilio_sites_all.pos -rf pumilio_sites_all.rf -doBcf 1 -doPost 1 --ignore-RG 0 -edge_pval 1e-4 -mapQ_pval 1e-4

#Do filtering with outgroup sylvatica
cat pumilio.r.bamlist sylvatica.r.bamlist > pumiSyl.r.bamlist
bamlist=pumiSyl.r.bamlist
/home/diana/bin/angsd -bam $bamlist -out pumilio_outSyl -minInd 340 -setMinDepthInd 1 -minMapQ 25 -minQ 25 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1  -GL 1 -doMaf 1 -doMajorMinor 4 -doGlf 2 -SNP_pval 1e-4 -nThreads 20 -ref $REF -setMaxDepthInd 50 -doCounts 1 -skipTriallelic 1 -dosnpstat 1 -doHWE 1 -sb_pval 1e-4 -doGeno 3 -hetbias_pval 1e-6 -sites pumilio_sites_all.pos -rf pumilio_sites_all.rf -doBcf 1 -doPost 1 --ignore-RG 0 -edge_pval 1e-4 -mapQ_pval 1e-4
