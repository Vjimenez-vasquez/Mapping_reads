# Mapping_reads
Mapping and estimating genome coverage, mean depths and other basic measures. 

A code designed by Victor Jimenez-Vasquez - vr.jimenez.vs@gmail.com

![Captura de pantalla de 2023-03-28 10-22-20](https://user-images.githubusercontent.com/89874227/228287232-802e9517-6f03-45d7-9f9c-3b8ed53e075b.png)


# 1: The code 
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

# 2: R summary #
```r
Rscript coverage_summary.R ; 
ls -lh ;
exit 
```

# 3: R Output 
```r
file : coverages.tsv
```

# 4: Extract aligned or unaligned reads (fastq format) "mapping_code.sh"
```r
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
```

# 5: "mapping_code.sh" usage
```r

./mapping_code.sh reference flag

example:

mapping_code.sh genome.fasta 3
genome = genome in fasta format
flag = number according to https://broadinstitute.github.io/picard/explain-flags.html
```

# 6: ACTUALIZACION DEL 21/07/26: "command3.sh"
```r
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
```
