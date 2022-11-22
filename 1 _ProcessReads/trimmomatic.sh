pop=$1
for sample in $pop*_R1_*fastq.gz
do 
name=${sample%_*R1*}
echo $name 
java -jar ~/bin/trimmomatic-0.39.jar PE -threads 5 $sample ${name}_R2_001.fastq.gz $name.trim.R1.fq.gz output_forward_unpaired.fq.gz $name.trim.R2.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:75
rm *unpaired*
done
