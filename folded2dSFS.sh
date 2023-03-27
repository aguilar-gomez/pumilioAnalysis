Calculate 2d folded SFS
for pop1 in AG CM DB PO CL HP SC TB
do 
for pop2  in AG CM DB PO CL HP SC TB
do 
if   [ $pop1 \< $pop2 ]
then
echo $pop1 $pop2
#2d sfs for pop1 and pop2 doing proper folding
nohup realSFS $pop1.r.saf.idx $pop2.r.saf.idx -P 1 -fold 1 1>$pop1$pop2.2dsfs.sfs 2>$pop1$pop2.2dsfs.sfs.out &
fi
done 
done
