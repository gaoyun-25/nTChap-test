#!/bin/bash

SRR_Acc_List="/home/gaoyun/poly/code/phasing/process_reads/potato/dowloadData/SRR_Acc_List.txt"
output="/home/gaoyun/poly/data/potato/"
mycode="/home/gaoyun/poly/code/phasing/process_reads/potato/dowloadData/"
threads=16

# SRR10489263 PromethION Promethion_1 OXFORD_NANOPORE；其他都是MinION MinIon_i 1..18
# SRR10489264 NextSeq 500 Illumina_1
promethion="SRR10489263"
illumina="SRR10489264"

# 提取chr2的short bam和long bam; 
# chrs_arr=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12)
chrs_arr=(chr02)
# Reads_name=(longReadsOnt promethion $illumina)
longReads_name=(longReadsOnt)
shortReads_name=($illumina)

SRR_names=`cat $SRR_Acc_List`

out_dirShort=$output"rawData/shortReads/"
out_dirLong=$output"rawData/longReads/"
[ -d $out_dirShort ] | mkdir -p $out_dirShort
[ -d $out_dirLong ] | mkdir -p $out_dirLong

# reference
ref_path=$output"reference/potato/"
ref_full_file=$ref_path"DM_1-3_516_R44_potato_genome_assembly.v6.1.fa.gz"
ref_file=$ref_path"potato.fa"
refChr02_file=$ref_path"potato_chr02.fa"
[ -d $ref_path ] | mkdir -p $ref_path
# 只提取12条染色体 chr01...chr12
if [ ! -f $ref_file ];then
    python $mycode"extract_seq.py" $ref_full_file $ref_file chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12
fi
if [ ! -f $refChr02_file ];then
    python $mycode"extract_seq.py" $ref_full_file $refChr02_file chr02
fi
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
if [ ! -f $refChr02_file".fai" ];then
    samtools faidx $refChr02_file
fi
if [ ! -f $refChr02_file".bwt" ];then
    bwa index $refChr02_file
fi
if [ ! -f $refChr02_file".dict" ];then
    gatk CreateSequenceDictionary -R $refChr02_file
fi

mappedLongReads=$output"Mapped/longReads/"
mappedShortReads=$output"Mapped/shortReads/"
[ -d $mappedLongReads ] | mkdir -p $mappedLongReads
[ -d $mappedShortReads ] | mkdir -p $mappedShortReads

shortReadVariants=$output"Variants/shortReads/"
longReadVariants=$output"Variants/longReads/"
[ -d $shortReadVariants ] | mkdir -p $shortReadVariants
[ -d $longReadVariants ] | mkdir -p $longReadVariants

SRR_list=""
strainName="potato"

# 合并 raw long reads fq
longReadsName=longReadsOnt
longReadsMapPrefix=$mappedLongReads$longReadsName
rawLongReadsFq=$out_dirLong$longReadsName"_1.fastq.gz"
if [ ! -f $rawLongReadsFq ];then
    # zcat $SRR_list | gzip - > $rawLongReadsFq
    cat $SRR_list > $rawLongReadsFq
fi
if [ ! -f $longReadsMapPrefix".sorted.bam" ];then
    minimap2 -ax map-ont -t $threads $ref_file $rawLongReadsFq | samtools view -bS - | samtools sort - -o $longReadsMapPrefix".sorted.bam"
    samtools index $longReadsMapPrefix".sorted.bam"
fi

annotation_path=$ref_path"annotation/"
[ -d $annotation_path ] | mkdir -p $annotation_path
gff3_file=$annotation_path"DM_1-3_516_R44_potato.v6.1.working_models.gff3"
if [ ! -f $gff3_file ];then
    gunzip $gff3_file".gz"
fi

# chr02的short bam和long bam 加tag;
for chr in ${chrs_arr[@]};do
    for reads_name in ${longReads_name[@]};do
        longReadsMapPrefix=$mappedLongReads$reads_name
        longReadsMapChrPrefix=$mappedLongReads$reads_name"."$chr
        # 需要.bai
        bai_file=$longReadsMapPrefix".sorted.bam.bai"
        longReadsMapFile=$longReadsMapPrefix".sorted.bam"
        if [ ! -f $bai_file ];then
            samtools index -@ $threads $longReadsMapFile
        fi

        longReadsMapFinalFile=$longReadsMapChrPrefix".sorted.bam"
        if [ ! -f $longReadsMapFinalFile ];then
            longReadsMapChrPassFile=$longReadsMapChrPrefix".pass.bam"
            # remove unmapped reads and extract alignments mapped to chr02
            samtools view -h -b -F 4 $longReadsMapFile $chr -@ $threads > $longReadsMapChrPassFile
            samtools sort $longReadsMapChrPassFile -@ $threads -o $longReadsMapFinalFile
            samtools index -@ $threads $longReadsMapFinalFile
            rm  $longReadsMapChrPassFile
        fi
        longReadsMapSam=$longReadsMapChrPrefix".sorted.sam"
        if [ ! -f $longReadsMapSam ];then
            samtools view -h -@ $threads $longReadsMapFinalFile -o $longReadsMapSam
        fi

        # convert gff to bed
        gene_file=$annotation_path"all_genes.bed"
        if [ ! -f $gene_file ];then
            grep ^chr $gff3_file | awk '$3 == "gene" {printf $1"\t"$4"\t"$5"\n" } ' - | sortBed -i - > $gene_file
        fi
        # split bed by chromosome
        geneChr_file=$annotation_path"all_genes."$chr".regions"
        if [ ! -f $geneChr_file ];then
            mergeBed -i $gene_file | awk '$1 == "chr02" {printf $1":"$2"-"$3"\n"}' > $geneChr_file
        fi
        # call variants
        vcf=$longReadVariants"variants."$chr".vcf"
        if [ ! -f $vcf ];then
            $mycode"freebayes-parallel" $geneChr_file $threads --ploidy 4 -f $refChr02_file $longReadsMapFinalFile > $vcf
        fi
    done
done


benchmarkingParameters=$mycode"benchmarkingParameters_potato_rawGene.tsv"
python $mycode"phaseToolRunner_gene.py" $benchmarkingParameters ntchap

python $mycode"plot_given_gene_id.py"
python $mycode"count_snps.py"
