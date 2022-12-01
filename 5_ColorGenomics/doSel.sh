#Selscan
#selscan genotypematrix.dgm fmatrix cmatrix
#Number of markers m
#Number of components k
#Fmatrix has the allele frequencies within each component (k x m)
#Cmatrix has the covariance matrix between components ([k-1] x [k-1])

#General selection, to do it specific for a population, a constant needs to be added:
#https://github.com/jade-cheng/ohana/wiki/Population-or-ancestry-specific-selection-scan

for k in $(seq 0 6);
do
cmatrixK.py c.matrix_k7 $k 100
done


nohup $OHANA/qpas pumiliov4maf5.lgm -k 7 -qi pumilio1_k7_e.08_q.matrix -fo pumilio_fullv4.f.matrix -e 0.08 -mi 450 >out.qpas.full &

OHANA=/home/diana/programs/ohana/bin
for k in $(seq 0 6);
do 
echo "running selscan with selection in component k $k"
$OHANA/selscan pumiliov4maf5.lgm pumilio_fullv4.f.matrix c.matrix_k7 -cs c.matrix_k7selK${k}_h10 > scan.pumiliok$k.txt &
done

