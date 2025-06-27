#!/bin/bash

output="/home/gaoyun/poly/data/potato/"
SRR_Acc_List="/home/gaoyun/poly/data/dowloadData/SRR_Acc_List.txt"

prefetch -O $output --max-size 150G --option-file $SRR_Acc_List

SRR_names=`cat $SRR_Acc_List`
# SRR10489263 PromethION Promethion_1 OXFORD_NANOPORE；其他都是MinION MinIon_i 1..18
# SRR10489264 NextSeq 500 Illumina_1
illumina="SRR10489264"
out_dirShort=$output"rawData/shortReads/"
out_dirLong=$output"rawData/longReads/"
[ -d $out_dirShort ] | mkdir -p $out_dirShort
[ -d $out_dirLong ] | mkdir -p $out_dirLong

for SRR_name in ${SRR_names[@]};do
    sra_file=$output$SRR_name"/"$SRR_name".sra"

    if [ "$SRR_name" == "$illumina" ];then
        parallel-fastq-dump -t 32 -O $out_dirShort --split-files --gzip -s $sra_file
    else
        fastq-dump $sra_file --split-3 --gzip --defline-qual '+' -A $SRR_name -O $out_dirLong
    fi
done

# download reference and gff3 file
ref_path=$output"reference/potato/"
[ -d $ref_path ] | mkdir -p $ref_path
#  https://spuddb.uga.edu/dm_v6_1_download.shtml
# DM_1-3_516_R44_potato_genome_assembly.v6.1.fa
# DM_1-3_516_R44_potato.v6.1.working_models.gff3
