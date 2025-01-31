#!/bin/bash
build=SGD
index1="/home/Database/bowtie-indexes/C_glabrata"
index2="/home/Database/bowtie-indexes/S_cerevisiae"
run="230607"
post="-n2-k1"
param="--chunkmbs 2048 -p8 -S"

func(){
    bowtie $index1 -1 $fq1 -2 $fq2 $param --un Cgla-un_fastq/Cgla_un_${run}_${prefix}.fastq
    bowtie $index2 -1 Cgla-un_fastq/Cgla_un_${run}_${prefix}_1.fastq -2 Cgla-un_fastq/Cgla_un_${run}_${prefix}_2.fastq $param > sam/Scer/${prefix}_Scer.sam
}

for prefix in 54I1 54I2 54W1 54W2
do
    fq1=${fastq}/${run}_${prefix}_1.fastq
    fq2=${fastq}/${run}_${prefix}_2.fastq
    func $prefix &>> logs/bowtie_Scer_$prefix.log
done

rm *~
