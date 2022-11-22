pop=$1
ref=$2 #fasta file
for sample in $pop**.filter_1*
do
name=${sample%.filter*}
name_o=${sample%_*}
echo $name
bwa mem -t 5 $ref $name_o"_1"* $name_o"_2"* >$name.sam
samtools view -Sb -F 1804 $name.sam|samtools sort -@ 5 - >$name.rescaffold.bam
samtools index $name.rescaffold.bam
    rm $name.sam
done
