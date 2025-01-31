database="/home/Database/others/S_cerevisiae"
gt="/home/Database/others/S_cerevisiae/genome_table"
index="$database/bowtie-indexes/S_cerevisiae"
run="230607"
name="Sao"
post="-n2-k1"
bin=100
genes="$database/SGD_features.tab"
ARS="$database/ARS-oriDB.txt"
TER="$database/TERs.txt"

for prefix in 54W1 54W2 54I1 54I2

do
    for command in "parse2wig+ -i bam/${prefix}_Scer.sorted.bam -o $prefix --gt $gt -n GR --nrpm 1500000 --ncmp 1500000 --pair %>> parse2wig+.log"
    do
        echo $command
        eval $command
    done
done


rm *~



