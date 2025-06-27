pop=$1
thetaStat print $pop.thetas.idx  >$pop.thetas.persite.txt


for pop in DBB DBR CMY CMR CL SK
do
nohup sh thetaPrint.sh $pop > out.thetaPrint$pop &
done

for pop in CL SK
do
nohup sh thetaPrint.sh $pop > out.thetaPrint$pop &
done




DBBpopsize=32
DBRpopsize=31
CMYpopsize=30
CMRpopsize=33
CLpopsize=15
SKpopsize=30

winsize=1000
for pop in DBB DBR CMY CMR CL SK
do
size="${pop}popsize"
value="${!size}"
echo $size $value
nohup ~/github_repo/angsdSupplement/TajimasD_SNPwindows.py $pop.thetas.persite.txt $pop.Tajimas.SNP$winsize $winsize $value &
done
