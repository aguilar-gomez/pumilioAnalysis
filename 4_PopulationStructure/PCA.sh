#Generate beagle
angsd -bam $bamlist -out pumilio_v4 -setMinDepthInd 1 -minMapQ 25 -minQ 25 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -GL 1 -doMaf 1 -doMajorMinor 4 -doGlf 2 -nThreads 10 -ref $REF -doCounts 1 -sites pumilio_oksites_angsdv4 -rf pumilio_oksites_angsdv4.rf

#Run PCAngsd with 10 million sites
python ~/bin/pcangsd.py -beagle pumilio_v4.beagle.gz -o pca_pumilio_v4 -threads 16 &


