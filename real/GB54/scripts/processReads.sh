#!/bin/bash

speciesName="Brettanomyces_bruxellensis"
output="/home/gaoyun/poly/data/GB54/"
mycode="/home/gaoyun/poly/code/phasing/process_reads/GB54/"
threads=32


out_dirShort=$output"rawData/shortReads/"
out_dirLong=$output"rawData/longReads/"
[ -d $out_dirShort ] | mkdir -p $out_dirShort
[ -d $out_dirLong ] | mkdir -p $out_dirLong

ref_path=$output"reference/"$speciesName"/"
[ -d $ref_path ] | mkdir -p $ref_path

ref_file=$ref_path"Brettanomyces_bruxellensis.fasta"

# index ref
if [ ! -f $ref_file".fai" ];then
    samtools faidx $ref_file
fi
if [ ! -f $ref_file".bwt" ];then
    bwa index $ref_file
fi
if [ ! -f $ref_file".dict" ];then
    gatk CreateSequenceDictionary -R $ref_file
fi

mappedLongReads=$output"Mapped/longReads/"
mappedShortReads=$output"Mapped/shortReads/"
[ -d $mappedLongReads ] | mkdir -p $mappedLongReads
[ -d $mappedShortReads ] | mkdir -p $mappedShortReads

shortReadVariants=$output"Variants/shortReads/"
longReadVariants=$output"Variants/longReads/"
[ -d $shortReadVariants ] | mkdir -p $shortReadVariants
[ -d $longReadVariants ] | mkdir -p $longReadVariants

ShortReadsName="ERR4624299"
rawShortReadsPrefix=$out_dirShort$ShortReadsName
shortReadsInitMapPrefix=$mappedShortReads$ShortReadsName".init"
shortReadsInitMap=$shortReadsInitMapPrefix".sorted.bam"
if [ ! -f $shortReadsInitMap ];then
    bwa mem -t $threads -M $ref_file $rawShortReadsPrefix"_1.fastq.gz" $rawShortReadsPrefix"_2.fastq.gz" | samtools view -bS - | samtools sort - -o $shortReadsInitMap
    samtools index $shortReadsInitMap
fi

shortReadsMapPrefix=$mappedShortReads$ShortReadsName
# 需要.bai
bai_file=$shortReadsInitMapPrefix".sorted.bam.bai"
if [ ! -f $bai_file ];then
    # samtools sort $shortReadsMapFile -@ $threads -o $shortReadsMapSortFile
    samtools index -@ $threads $shortReadsMapSortFile
fi

shortReadsMapFinalFile=$shortReadsMapPrefix".sorted.bam"
if [ ! -f $shortReadsMapFinalFile ];then
    shortReadsMapRGFile=$shortReadsMapPrefix".RG.bam"
    shortReadsMapRGSortFile=$shortReadsMapPrefix".RG.sorted.bam"
    shortReadsMapRGSortMDFile=$shortReadsMapPrefix".RG.sorted.MD.bam"
    gatk AddOrReplaceReadGroups -I $shortReadsInitMap -O $shortReadsMapRGFile -RGID "ID_"$speciesName -RGLB "LB_"$speciesName -RGPL "SHORTREADS" -RGPU "PU_"$speciesName -RGSM "SM_"$speciesName
    samtools sort $shortReadsMapRGFile -@ $threads -o $shortReadsMapRGSortFile
    samtools index -@ $threads $shortReadsMapRGSortFile
    gatk MarkDuplicates --REMOVE_DUPLICATES true -I $shortReadsMapRGSortFile -O $shortReadsMapRGSortMDFile -M "/dev/null"
    samtools sort $shortReadsMapRGSortMDFile -@ $threads -o $shortReadsMapFinalFile
    samtools index -@ $threads $shortReadsMapFinalFile
    rm $shortReadsInitMap $shortReadsMapRGFile $shortReadsMapRGSortMDFile $shortReadsMapRGSortFile
fi

shortReadVCF=$shortReadVariants$ShortReadsName".vcf"
ploidy=3
if [ ! -f $shortReadVCF ];then
    gatk HaplotypeCaller -R $ref_file -ploidy $ploidy -I $shortReadsMapFinalFile -O $shortReadVCF --native-pair-hmm-threads $threads
fi


longReadsName="ERR4624298"
# longReadsInitMapPrefix=$mappedLongReads$longReadsName".init"
longReadsInitMapPrefix=$mappedLongReads$longReadsName
longReadsInitMapFile=$longReadsInitMapPrefix".sorted.bam"
if [ ! -f $longReadsInitMapFile ];then
    minimap2 -ax map-ont -t $threads $ref_file $out_dirLong$longReadsName".fastq.gz" | samtools view -bS - | samtools sort - -o $longReadsInitMapFile
    samtools index $longReadsInitMapFile
fi
longReadsMapChrsSam=$mappedLongReads$longReadsName".sorted.sam"
if [ ! -f $longReadsMapChrsSam ];then
    samtools view -h -@ $threads $longReadsInitMapFile -o $longReadsMapChrsSam
fi


# run our phasing
benchmarkingParameters=$mycode"benchmarkingParameters_GB54.tsv"
python $mycode"phase_tool_runners_GB54.py" $benchmarkingParameters ntchap plot Phasing_clean_chr plot_clean_chr


