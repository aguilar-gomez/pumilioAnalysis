pop=$1
for sample in $pop*trim.R1.fq.gz;
do  
name=${sample%%.*}
echo $name 
gunzip $name.trim.R1.fq.gz $name.trim.R2.fq.gz
     prinseq-lite.pl -fastq $name.trim.R1.fq -fastq2 $name.trim.R2.fq -lc_method dust -lc_threshold 7 -custom_params "G 50" -out_good $name.filter
     fastqc -o fastqc_prinseq -f fastq -t 5 $name.filter_1.fastq
     fastqc -o fastqc_prinseq -f fastq -t 5 $name.filter_2.fastq
rm *bad* *single*
gzip $name.trim.R1.fq $name.trim.R2.fq
 
 
done
