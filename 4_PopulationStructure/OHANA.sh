
#OHANA
OHANA=/home/diana/programs/ohana/bin
nohup $OHANA/convert bgl2lgm ./pumilio_v4_maf5.beagle ./pumiliov4maf5.lgm &

#1%
/home/diana/programs/ohana/tools/sample-sites.py  ./pumiliov4maf5.lgm  1 ./pumiliov4maf5onepercent.lgm &
########################################
Running qpas
Kmin=4 # min value of k
Kmax=9 # maximum value of k

## run ohana structure (qpas)
OHANA=/home/diana/programs/ohana/bin
for k in $(seq $Kmin $Kmax);
    do echo "running qpas with k $k"
    $OHANA/qpas pumiliov4maf5onepercent.lgm -k $k -qo pumilio1_k${k}_e.08_q.matrix -fo pumilio1_k${k}_e.08_f.matrix -e 0.08 -mi 450 >out1.k${k} &
done
########################################
nohup sh qpas.sh > out &

########################################
Running nemeco and tree
Kmin=4 # min value of k
Kmax=7 # maximum value of k
OHANA=/home/diana/programs/ohana/bin
for k in $(seq $Kmin $Kmax);
do
echo "running nemeco with k $k"
$OHANA/nemeco pumiliov4maf5onepercent.lgm pumilio1_k${k}_e.08_f.matrix -co c.matrix_k${k} -mi 5 > out1.nemk${k} &
$OHANA/convert cov2nwk c.matrix_k${k} pum1_k${k}.nwk
tail -n +2 pumilio1_k${k}_e.08_q.matrix  > pum_k${k}.Q
done
########################################
nohup sh donemeco.sh >out1.tree



########################################
pong -m filemap -n poporder2 -i popindividuals_3code_fix -l colores3
pong -m filemap -n popsOrder -i popindividuals_3code_fix -l colores
pong -m filemap89 -n popsOrderv3 -i popindividuals_3code_v2 -l colores3

########################################
On my computer:
ssh -N -L 4000:localhost:4000 diana@tilden.biol.berkeley.edu
Go to:
http://localhost:4000/

########################################
#Selscan
#selscan genotypematrix.dgm fmatrix cmatrix
#Number of markers m
#Number of components k
#Fmatrix has the allele frequencies within each component (k x m)
#Cmatrix has the covariance matrix between components ([k-1] x [k-1])

That is just general selection, to do it specific for a population, a constant needs to be added:
https://github.com/jade-cheng/ohana/wiki/Population-or-ancestry-specific-selection-scan

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
