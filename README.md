# Mapping_reads
Mapping and estimating genome coverage, mean depths and other basic measures. 

A code designed by Victor Jimenez-Vasquez - vr.jimenez.vs@gmail.com

![Captura de pantalla de 2023-03-28 10-22-20](https://user-images.githubusercontent.com/89874227/228287232-802e9517-6f03-45d7-9f9c-3b8ed53e075b.png)


# The code 
```r
#!/bin/bash

#1# indexar el genoma de referencia#
bwa index reference.fasta ;

#2# preparar las instrucciones generales#
for r1 in *fastq.gz
do
prefix=$(basename $r1 _L004_R1_001.fastq.gz)
r2=${prefix}_L004_R2_001.fastq.gz

#3# instrucciones para generar el archivo .bam#
bwa mem -t 15 reference.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 15 -bS -T reference.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 15 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -@ 15 -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort -@ 15 ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup -@ 15 ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index -@ 15 ${prefix}.bam ;

#4# remover los archivos intermediarios#
rm ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam ;
done ;
ls -lh ; 

#5# estimate genome coverage #
for r1 in *.bam
do
prefix=$(basename $r1 .bam)
samtools coverage $r1 -o ${prefix}.tsv
done ;
mkdir tsv ;
mv *.tsv tsv/ ;
cd tsv/ ;
grep -v "#" *tsv | tr "\t" "," | tr ":" "," | sed -e 's/.tsv//g' | sed '1i sample reference startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq' | tr " " "," > coverage.csv
cat coverage.tsv ;
ls ;
```

#6# R summary #
```r
Rscript coverage_summary.R ; 
ls -lh ;
exit 
```

# Output 
```r
file : coverages.tsv
```
