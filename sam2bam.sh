#!/bin/bash

mkdir -p bam

sam(){
    samtools view -S -b sam/Scer/${prefix}_Scer.sam -o bam/${prefix}_Scer.bam
    samtools sort bam/${prefix}_Scer.bam -o bam/${prefix}_Scer.sorted.bam
    samtools index bam/${prefix}_Scer.sorted.bam
}

for prefix in 54I1 54I2 54W1 54W2
do
    sam $prefix &>> logs/sam2bam_${prefix}.log
done

rm *~

    
