sing="singularity exec --bind /work,/work2 /work/SingularityImages/rnakato_ssp_drompa.2022.08.sif"
database="/home/Database/others/S_cerevisiae"
gt="/home/Database/others/S_cerevisiae/genome_table"
index="$database/bowtie-indexes/S_cerevisiae"
run="230721"
name="Sao"
post="-n2-k1"
bin=100
genes="$database/SGD_features.tab"
ARS="$database/ARS-oriDB.txt"
TER="$database/TERs.txt"


C1=74PK_W
C2=75PK_W

#IP                                                                                                                                             
C3=74PK_IP
C4=75PK_IP


#WCE
W1="bin.norm/$C1"
W2="bin.norm/$C2"

#IP
IP1="bin.norm/$C3" #74_dd
IP2="bin.norm/$C4" #75_wt


#input
s1="-i $IP1,$W1,$C3"
s2="-i $IP2,$W2,$C4"

odir=bg.norm
mkdir $odir
drompa_peakcall PC_ENRICH $s1 -p $run-PCE-sm1000_${C3}.norm -gt $gt -binsize $bin -sm 1000 -outputwig 1 -owtype 3 -odir $odir -norm 0
drompa_peakcall PC_ENRICH $s2 -p $run-PCE-sm1000_${C4}.norm -gt $gt -binsize $bin -sm 1000 -outputwig 1 -owtype 3 -odir $odir -norm 0

drompa_peakcall PC_ENRICH $s2 -p $run-PCE-sm1000_54I2.norm2 -gt $gt -binsize $bin -sm 1000
drompa_peakcall PC_ENRICH $s2 -p $run-PCE-sm1000_rpo21-74I -gt $gt -binsize $bin -sm 1000 -odir PCE -outputwig 1
drompa_draw PC_ENRICH $s1 $s2 $s3 -p $run.40.53.54.sm1000 -gt $gt -sm 1000 -ars $ARS -gftype 3 -gene $genes -ls 100 -rmchr -binsize $bin -scale_ratio 5 -bn 3

rm *~
