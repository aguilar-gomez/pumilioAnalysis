NUMBEROFSIMS=20   #PLACE 7

rm -r TemporaryFiles
mkdir TemporaryFiles
for (( c=1; c<=$NUMBEROFSIMS; c++ ))
do
python3.9 run_simulation.py
Rscript ConvertData.R
mv TemporaryFiles/Data.vcf TemporaryFiles/Data${c}.vcf
mv TemporaryFiles/TrueTree.txt TemporaryFiles/TrueTree${c}.txt
mv TemporaryFiles/adbayesinput.txt TemporaryFiles/adbayesinput${c}.txt
done

for (( c=1; c<=$NUMBEROFSIMS; c++ ))
do

python3.10 ~/desktop/admixturebayes/runMCMC.py --input_file TemporaryFiles/adbayesinput${c}.txt --outgroup outgroup  --n 20000    --max_admixes  0

python3.10 GetAdBayesResults.py

cp MAPadd.txt TemporaryFiles/MAPadd${c}.txt 
cp MAPtree.txt TemporaryFiles/MAPtree${c}.txt

for (( d=1; d<=100; d++ ))
do
cp TemporaryFiles/Tree${d}.txt TemporaryFiles/Tree${d}_${c}.txt 

done

done

for (( c=1; c<=$NUMBEROFSIMS; c++ ))
do

cp TemporaryFiles/TrueTree${c}.txt TemporaryFiles/TrueTree.txt
cp TemporaryFiles/MAPTree${c}.txt TemporaryFiles/MAPTree.txt

for (( d=1; d<=100; d++ ))
do
cp TemporaryFiles/Tree${d}_${c}.txt TemporaryFiles/Tree${d}.txt
done

Rscript CheckTopologyEqualityAd.R
cp TemporaryFiles/TopEq.txt TemporaryFiles/TopEq${c}.txt
done

rm -f AdBayesTop.txt
touch AdBayesTop.txt

for (( c=1; c<=$NUMBEROFSIMS; c++ ))
do

sed '2,2!d' TemporaryFiles/TopEq${c}.txt > temp.txt
cat AdBayesTop.txt temp.txt > tempp.txt 
mv tempp.txt AdBayesTop.txt 

done

rm temp.txt
rm mcmc_samples.csv

rm MAPadd.txt
rm MAPtree.txt

rm -r -f TemporaryFiles
