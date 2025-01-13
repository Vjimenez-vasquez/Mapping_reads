#!/usr/bin/bash

reference=$1
flag=$2
region=$(basename $reference .fasta)

#0# indexar el genoma de referencia (en caso del genoma humano, contar con todo el genoma indexado en una carpeta fija)#
bwa index $reference ;

#1# preparar las instrucciones generales#
for r1 in *fastq.gz
do
prefix=$(basename $r1 .R1.fastq.gz)
r2=${prefix}.R2.fastq.gz

#2# instrucciones para generar el archivo .bam#
bwa mem -t 4 $reference $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 4 -bS -T $reference ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 4 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -@ 4 -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort -@ 4 ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup -@ 4 ${prefix}_cuatroa.bam ${prefix}.${region}.bam ;
samtools index -@ 4 ${prefix}.${region}.bam ;

#3# remover los archivos intermediarios#
rm ${prefix}_uno.bam ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam ;

#4# obtencion de: 1)archivos bam conformado por solo reads "mapeados" o "no-mapeados" (dependiendo del flag elegido) y 2) fastq files "f" y "r" #
samtools view -@ 4 -b -f $flag ${prefix}.${region}.bam > ${prefix}.${region}.mapped.bam ;
samtools index -@ 4 ${prefix}.${region}.mapped.bam ;
samtools fastq -1 ${prefix}.${region}.f.fastq -2 ${prefix}.${region}.r.fastq -0 /dev/null -s /dev/null -n ${prefix}.${region}.mapped.bam ;
rm ${prefix}.R2.fastq.gz.${region}.r.fastq ${prefix}.R2.fastq.gz.${region}.f.fastq ${prefix}.R2.fastq.gz.${region}.mapped.bam ;

done ;
