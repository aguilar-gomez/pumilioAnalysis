#Fst index and scan use unfolded saf
pop1=$1
pop2=$2
#joint SFS
realSFS $pop1.r.saf.idx $pop2.r.saf.idx -P 5 >$pop1.$pop2.ml
#Fst index
realSFS fst index $pop1.r.saf.idx $pop2.r.saf.idx -sfs $pop1.$pop2.ml  -fstout $pop1.$pop2
realSFS fst stats $pop1.$pop2.fst.idx > fstindex$pop1$pop2
#Fst scan
realSFS fst stats2 $pop1.$pop2.fst.idx -win 100000 -step 20000 -type 2 > fst_w100kb$pop1$pop2


#PBS code and Windows of SNPs
#PBS and fst all together
pop1=$1
pop2=$2
pop3=$3
SNP=$4
realSFS fst index $pop1.r.saf.idx $pop2.r.saf.idx $pop3.r.saf.idx -sfs $pop1.$pop2.ml -sfs $pop1.$pop3.ml -sfs $pop2.$pop3.ml -fstout $pop1.$pop2.$pop3
realSFS fst stats  $pop1.$pop2.$pop3.fst.idx>pbs.$pop1.$pop2.$pop3
realSFS fst stats2 $pop1.$pop2.$pop3.fst.idx -win 100000 -step 20000 -type 2 >slidingwindow100Ks20k.$pop1.$pop2.$pop3
#Make SNP windows
realSFS fst print $pop1.$pop2.$pop3.fst.idx  > $pop1.$pop2.$pop3.printed.txt
PBS_SNP_windows.py $pop1.$pop2.$pop3.printed.txt $pop1.$pop2.$pop3.wSNP$SNP $SNP
