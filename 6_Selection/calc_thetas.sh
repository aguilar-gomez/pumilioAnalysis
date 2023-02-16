#Calculate site frequency spectrum
outdir=$1
sites_file=$2
regions=$3
pop=$4
REF=/space/s2/diana/pumilio/rescaffolded_genome/Opum.rescaffold.fasta
/home/diana/bin/angsd -nThreads 3 -bam $pop.r.bamlist -ref $REF -anc $REF -out $outdir/$pop.shared. \
                    	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -baq 1 \
                    	-minMapQ 25 -minQ 25 \ 
                    	-GL 1 -doSaf 1 -sites $sites_file -rf $regions\


#Tajimas D and thetas 
#The values of the different thetas are negative because they are logscale!

#Calculate thetas
realSFS saf2theta $pop.saf.idx -sfs $pop.folded.sfs -outname $pop
thetaStat print $pop.thetas.idx > $pop.thetas.persite.txt

#For Dxy
nohup realSFS -P 4 $POP.r.saf.idx -sites $pop1.$pop2.shared.pos -fold 1 1> Results/$POP.sfs 2>Results/$POP.outsfs &

#Calculate Dxy
$NGSTOOLS/ngsStat -npop 2 -postfiles $pop1.shared.saf $pop2.shared.saf -nsites $NSITES -nind $pop1_nind $pop2_nind -outfile $pop1.$pop2.stats.txt


#Do windows of SNPS
Tajimas_SNPmidpoint.py $pop.thetas.persite.shared.sorted $pop.Tajimas.midpoint.SNP$sites $sites $popsize &
SNPwindows_midpoint.py $pop1.$pop2.persite.Dxy $pop1.$pop2.SNP200.Dxy.midpoint 200 Dxy &
