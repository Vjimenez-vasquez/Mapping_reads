#!/usr/bin/bash

region=$1

if [ "$region" == "aligned" ]; then

for b1 in *.bam
do
prefix=$(echo $b1 | sed -e "s/.GC.*//g")
echo "sample: $prefix"
echo "extrayendo secuencias 'alineadas'"
echo " ..."
samtools view -@ 15 -b -f 1 -F 4 -F 8 $b1 > ${prefix}.${region}.bam ;
samtools sort -@ 15 -n ${prefix}.${region}.bam -o ${prefix}.${region}.sorted.bam ;
samtools fastq -@ 15 -1 ${prefix}.${region}.f.fastq -2 ${prefix}.${region}.r.fastq -0 /dev/null -s /dev/null -n ${prefix}.${region}.sorted.bam ;
rm ${prefix}.${region}.bam ${prefix}.${region}.bam.bai ${prefix}.${region}.sorted.bam ${prefix}.${region}.sorted.bam.bai ;
gzip ${prefix}.${region}.f.fastq ${prefix}.${region}.r.fastq ;
done ;

else


for b1 in *.bam
do
prefix=$(echo $b1 | sed -e "s/.GC.*//g")
echo "sample: $prefix"
echo "extrayendo secuencias 'no-alineadas'"
echo " ..."
samtools view -@ 15 -b -f 4 $b1 > ${prefix}.${region}.bam ;
samtools sort -@ 15 -n ${prefix}.${region}.bam -o ${prefix}.${region}.sorted.bam ;
samtools fastq -@ 15 -1 ${prefix}.${region}.f.fastq -2 ${prefix}.${region}.r.fastq -0 /dev/null -s /dev/null -n ${prefix}.${region}.sorted.bam ;
rm ${prefix}.${region}.bam ${prefix}.${region}.bam.bai ${prefix}.${region}.sorted.bam ${prefix}.${region}.sorted.bam.bai ;
gzip ${prefix}.${region}.f.fastq ${prefix}.${region}.r.fastq ;
done ;

fi

## flags "-f 1 -F 4 -F 8" para obtener reads mapeados ##
## flags "-f 4 " para obtener reads no-mapeados ##
## ./command3.sh aligned ##
## ./command3.sh no_aligned ##