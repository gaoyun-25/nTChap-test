import os
import subprocess
from Bio import SeqIO
import gzip
import matplotlib.pyplot as plt
import math
from collections import defaultdict
import pandas as pd
import numpy as np

# Distribution of reduced reads length, reduced reads coverage, Heterozygous SNPs Number of each chromosome

def count_reads_length(long_reads_fq, read_lengths_txt):
    with gzip.open(long_reads_fq, "rt") as f_in, open(read_lengths_txt, "w") as f_out:
        # f_out.write("ReadID\tLength\n")
        for record in SeqIO.parse(f_in, "fastq"):
            read_id = record.id.split()[0]
            read_length = len(record.seq)
            f_out.write(f"{read_id} {read_length}\n")

def plot_reads_lengths_hist(read_lengths_txt, out_path, out_name, bin_width):
    lengths = []
    with open(read_lengths_txt) as f:
        for line in f:
            line = line.strip().split(" ")
            lengths.append(int(line[1]))

    # 绘制直方图（对数刻度）
    max_val = max(lengths)
    min_val = 0
    bins = np.arange(min_val, max_val+bin_width, bin_width)
    reads_lenths_png = os.path.join(out_path, out_name+"_long_reads_length.png")
    reads_lenths_svg = os.path.join(out_path, out_name+"_long_reads_length.svg")
    plt.figure(figsize=(18, 10))
    plt.hist(lengths, bins=bins, log=True, edgecolor="black", alpha=0.7, range=(0, max(lengths)))
    plt.xlabel("Read Length (bp)")
    plt.ylabel("Count (log scale)")
    plt.title("Long Reads Length Distribution")

    # 自动格式化x轴刻度（科学计数法处理大数值）
    plt.ticklabel_format(axis='x', style='plain')  # 禁用科学计数法
    ax = plt.gca()

    # 旋转x轴标签并调整密度
    plt.xticks(rotation=45, ha='right', fontsize=10)  # 45度旋转，右对齐

    # 自动调整布局
    plt.tight_layout(pad=3)  # 增加布局边距

    # 手动调整子图位置（根据需要微调）
    plt.subplots_adjust(bottom=0.18)  

    plt.savefig(reads_lenths_png, dpi=300, bbox_inches="tight")
    plt.savefig(reads_lenths_svg, dpi=300, bbox_inches='tight')
    plt.close()

def plot_chrs_reads_lengths_hist(out_path, out_name, bam_file, chrs, bin_width):

    for chr in chrs:
        out_chr_name = out_name+"."+chr
        chr_reads_id_txt = os.path.join(out_path, out_chr_name+"_reads_id.txt")
        os.system("samtools view "+bam_file+" "+chr+" -F 4 | awk '{print $1}' > "+chr_reads_id_txt)
        # 从BAM文件提取已筛选reads的长度
        chr_read_lengths_txt = os.path.join(out_path, out_chr_name+"_read_lengths.txt")
        os.system("samtools view "+bam_file+" | grep -wFf "+chr_reads_id_txt+" | awk '{print $1, length($10)}' > "+chr_read_lengths_txt)
        plot_reads_lengths_hist(chr_read_lengths_txt, out_path, out_chr_name, bin_width)

def count_chrs_intervals_snps(hetSNPs_positions_file, intervals_length,out_path, out_name):
    # 按照chr长度进行划分，看固定长度区间多少het snps，是否有非常少的区域
    chr_pos =  dict()
    with open(hetSNPs_positions_file) as f:
        for line in f:
            line = line.strip().split("\t")
            chr_pos.setdefault(line[0], set())
            chr_pos[line[0]].add(int(line[1]))

    # chr_length = dict()
    # with open(het_file) as f:
    #     next(f)
    #     for line in f:
    #         line = line.strip().split("\t")
    #         chr_length[line[0]] = line[1]

    # 每个interval有多少个het snps
    intervals_nSNPs = dict()
    out_nSNPs_file = os.path.join(out_path, out_name+"_block_nSNPs.tsv")
    for chr, positions in chr_pos.items():
        intervals_nSNPs.setdefault(chr, dict())
        max_pos = max(positions)+ intervals_length
        for i in range(0, max_pos, intervals_length):
            for pos in positions:
                if pos >= i and pos < i+intervals_length:
                    block = i + intervals_length/2
                    intervals_nSNPs[chr].setdefault(block, 0)
                    intervals_nSNPs[chr][block] += 1

    with open(out_nSNPs_file, "w") as f:
        f.write("chromosome\tcoordinate\tnSNPs\n")
        for chr, blocks_nsnps in intervals_nSNPs.items():
            for block in sorted(blocks_nsnps):
                f.write(chr+"\t"+str(block)+"\t"+str(blocks_nsnps[block])+"\n")

    # plot_bar_snps(out_nSNPs_file,out_path,out_name+"."+chr)


def plot_bar_snps(chrs_nSNPs_file,out_path,out_name):
    df=pd.read_csv(chrs_nSNPs_file,sep="\t")
    for chr, data in df.groupby('chromosome'):
        out_png=os.path.join(out_path,out_name+"."+chr+"_snps_bar.png")
        out_svg=os.path.join(out_path,out_name+"."+chr+"_snps_bar.svg")

        x = data['coordinate'].values
        y = data['nSNPs'].values

        # 确保位置数据按升序排列
        sorted_indices = np.argsort(x)
        x = x[sorted_indices]
        y = y[sorted_indices]

        # 计算每个柱子的宽度（当前位置到下一个位置的间隔）
        if len(x) > 1:
            widths = np.diff(x)
            widths = np.append(widths, widths[-1])  # 最后一个宽度沿用前一个
        else:
            widths = np.array([1])  # 若只有一个数据点，默认宽度为1

        # 绘制柱状图
        plt.figure(figsize=(12, 6))
        plt.bar(x, y, width=widths, align='edge', edgecolor='black')

        # 设置图表标题和轴标签
        plt.xlabel('Chromosome Positions')
        plt.ylabel('Heterozygous SNPs Number')

        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.savefig(out_svg, dpi=300, bbox_inches='tight')
        plt.close()

def plot_chrs_intervals_snps(hetSNPs_positions_file, intervals_length,out_path, out_name):
    chr_pos =  dict()
    with open(hetSNPs_positions_file) as f:
        for line in f:
            line = line.strip().split("\t")
            chr_pos.setdefault(line[0], set())
            chr_pos[line[0]].add(int(line[1]))

    for chr, chr_len in chr_pos.items():
        hist_name = out_name+"."+chr+"_chr_snps_hist"
        plot_chrs_nSNPs_hist(chr_len, out_path, hist_name, intervals_length)

def plot_chrs_nSNPs_hist(data, out_path, out_name, bin_width):
    reads_lenths_png = os.path.join(out_path, out_name+".png")
    reads_lenths_svg = os.path.join(out_path, out_name+".svg")
    plt.figure(figsize=(12, 6))
    max_val = max(data)
    min_val = 0
    bins = np.arange(min_val, max_val+bin_width, bin_width)
    plt.hist(data, bins=bins, edgecolor="black", alpha=0.7, range=(0, max(data)))
    plt.xlabel("Chromosome Positions", fontsize=18)
    plt.ylabel("Heterozygous SNPs Number", fontsize=18)

    # 调整X轴和Y轴刻度标签大小
    plt.xticks(fontsize=14)  # X轴刻度标签字号
    plt.yticks(fontsize=14)  # Y轴刻度标签字号

    # 自动格式化x轴刻度（科学计数法处理大数值）
    # plt.ticklabel_format(axis='x', style='plain')  # 禁用科学计数法
    ax = plt.gca()

    # 自动调整布局
    plt.tight_layout(pad=3)  # 增加布局边距

    # 手动调整子图位置（根据需要微调）
    plt.subplots_adjust(bottom=0.18)  

    plt.savefig(reads_lenths_png, dpi=300, bbox_inches="tight")
    plt.savefig(reads_lenths_svg, dpi=300, bbox_inches='tight')
    plt.close()

def plot_nSNPs_hist(het_file, out_path, out_name,bin_width):
    chr_length = dict()
    with open(het_file) as f:
        next(f)
        for line in f:
            line = line.strip().split("\t")
            chr_length[line[0]] = line[1]

    for chr, chr_len in chr_length.items():
        block_nSNPs_file = os.path.join(out_path, out_name+"."+chr+"_block_nSNPs.tsv")
        positions = []
        with open(block_nSNPs_file) as f:
            next(f)
            for line in f:
                line = line.strip().split("\t")
                positions.append(int(line[1]))
        out_block_name = out_name+"."+chr+"_nSNPs_freq"
        # plot_chrs_nSNPs_hist(positions, out_path, out_block_name, bin_len)
        reads_lenths_png = os.path.join(out_path, out_block_name+".png")
        reads_lenths_svg = os.path.join(out_path, out_block_name+".svg")
        plt.figure(figsize=(18, 10))
        max_val = max(positions)
        min_val = 0
        bins = np.arange(min_val, max_val+bin_width, bin_width)
        plt.hist(positions, bins=bins, edgecolor="black", alpha=0.7, range=(0, max(positions)))
        plt.xlabel("nPos")
        plt.ylabel("nBlocks")
        # plt.title("Long Reads Length Distribution")

        # 自动格式化x轴刻度（科学计数法处理大数值）
        # plt.ticklabel_format(axis='x', style='plain')  # 禁用科学计数法
        ax = plt.gca()

        # # 旋转x轴标签并调整密度
        # plt.xticks(rotation=45, ha='right', fontsize=10)  # 45度旋转，右对齐

        # # 自动调整布局
        # plt.tight_layout(pad=3)  # 增加布局边距

        # # 手动调整子图位置（根据需要微调）
        # plt.subplots_adjust(bottom=0.18)  

        plt.savefig(reads_lenths_png, dpi=300, bbox_inches="tight")
        # plt.savefig(reads_lenths_svg, dpi=300, bbox_inches='tight')
        plt.close()


def snpReads_coverage(snp_reads_file, intervals_length,snpReads_coverage_file):
    snps_reads = dict()
    with open(snp_reads_file) as f:
        for line in f:
            line = line.strip().split("\t")
            if line[0][0] == ">":
                chr = line[0][1:]
            else:
                read_name = line[0]
                snps = line[1:]
                snps_reads.setdefault(chr,dict())
                for snp in snps:
                    pos = int(snp.split("=")[0])
                    snps_reads[chr].setdefault(pos, 0)
                    snps_reads[chr][pos] += 1

    mean_coverages = dict()
    for chr, pos_cov in snps_reads.items():
        mean_coverages.setdefault(chr, dict())
        max_pos = max(pos_cov)+ intervals_length
        for i in range(0, max_pos, intervals_length):
            block_coverages = list()
            for pos, cov in pos_cov.items():
                if pos >= i and pos < i + intervals_length:
                    block_coverages.append(cov)
            if len(block_coverages) == 0:
                mean = 0
            else:
                mean = sum(block_coverages)/len(block_coverages)
            block = i + intervals_length/2
            mean_coverages[chr][block] = mean

    with open(snpReads_coverage_file, "w") as f:
        f.write("chromosome\tcoordinate\tcoverage\n")
        for chr,block_cov in mean_coverages.items():
            block_sort = sorted(block_cov)
            for block in block_sort:
                f.write(chr+"\t"+str(block)+"\t"+str(block_cov[block])+"\n")

def plot_bar_snpReads_coverage(snpReads_coverage_file,out_path,out_name):
    df=pd.read_csv(snpReads_coverage_file,sep="\t")
    for chr_name, data in df.groupby('chromosome'):
        out_png = os.path.join(out_path, out_name+"."+chr_name+".snpReads_coverage_bar.png")
        out_svg = os.path.join(out_path, out_name+"."+chr_name+".snpReads_coverage_bar.svg")

        x = data['coordinate'].values
        y = data['coverage'].values

        # 确保位置数据按升序排列
        sorted_indices = np.argsort(x)
        x = x[sorted_indices]
        y = y[sorted_indices]

        # 计算每个柱子的宽度（当前位置到下一个位置的间隔）
        if len(x) > 1:
            widths = np.diff(x)
            widths = np.append(widths, widths[-1])  # 最后一个宽度沿用前一个
        else:
            widths = np.array([1])  # 若只有一个数据点，默认宽度为1

        # 绘制柱状图
        plt.figure(figsize=(12, 6))
        plt.bar(x, y, width=widths, align='edge', edgecolor='black')

        # 设置图表标题和轴标签
        plt.xlabel('Chromosome Positions', fontsize=18)
        plt.ylabel('Coverage', fontsize=18)
        # 调整X轴和Y轴刻度标签大小
        plt.xticks(fontsize=14)  # X轴刻度标签字号
        plt.yticks(fontsize=14)  # Y轴刻度标签字号

        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.savefig(out_svg, dpi=300, bbox_inches='tight')
        plt.close()


def count_snpReads_length(snp_reads_file, intervals_length, snpReads_length_file):
    snpReads = dict()
    with open(snp_reads_file) as f:
        for line in f:
            line = line.strip().split()
            if line[0][0] == ">":
                chr = line[0][1:]
            else:
                read_name = line[0]
                snps = line[1:]
                snpReads.setdefault(chr,list())
                positions = set()
                for snp in snps:
                    pos = int(snp.split("=")[0])
                    positions.add(pos)
                positions = sorted(positions)
                snpReads[chr].append(positions)
                # chr:[[pos1..],[pos1...]]

    # 取最大pos
    snpReads_length = dict()
    for chr, positions_list in snpReads.items():
        snpReads_length.setdefault(chr,dict())
        max_list = []
        for positions in positions_list:
            max_list.append(positions[-1])
        max_pos = max(max_list) + intervals_length
        for i in range(0, max_pos, intervals_length):
            chr_snpReads_length = list()
            block = i + intervals_length/2
            for positions in positions_list:
                n = 0
                right = i + intervals_length
                for pos in positions:
                    if pos >= i and pos < right:
                        n += 1
                if n > 0:
                    chr_snpReads_length.append(n)
            if len(chr_snpReads_length) == 0:
                mean = 0
            else:
                mean = sum(chr_snpReads_length)/len(chr_snpReads_length)
            snpReads_length[chr][block] = mean

    with open(snpReads_length_file, "w") as f:
        f.write("chromosome\tcoordinate\tlengths\n")
        for chr, block_length in snpReads_length.items():
            order_pos = sorted(block_length)
            for b in order_pos:
                f.write(chr+"\t"+str(b)+"\t"+str(block_length[b])+"\n")


def plot_bar_snpReads_lengths(snpReads_lengths_file,out_path,out_name):
    df=pd.read_csv(snpReads_lengths_file,sep="\t")
    for chr_name, data in df.groupby('chromosome'):
        out_png = os.path.join(out_path, out_name+"."+chr_name+".snpReads_lengths_bar.png")
        out_svg = os.path.join(out_path, out_name+"."+chr_name+".snpReads_lengths_bar.svg")

        x = data['coordinate'].values
        y = data['lengths'].values

        # 确保位置数据按升序排列
        sorted_indices = np.argsort(x)
        x = x[sorted_indices]
        y = y[sorted_indices]

        # 计算每个柱子的宽度（当前位置到下一个位置的间隔）
        if len(x) > 1:
            widths = np.diff(x)
            widths = np.append(widths, widths[-1])  # 最后一个宽度沿用前一个
        else:
            widths = np.array([1])  # 若只有一个数据点，默认宽度为1

        # 绘制柱状图
        plt.figure(figsize=(12, 6))
        plt.bar(x, y, width=widths, align='edge', edgecolor='black')

        # 设置图表标题和轴标签
        plt.xlabel('Chromosome Position', fontsize=18)
        plt.ylabel('Mean of Reduced Reads Lengths', fontsize=18)
        # 调整X轴和Y轴刻度标签大小
        plt.xticks(fontsize=14)  # X轴刻度标签字号
        plt.yticks(fontsize=14)  # Y轴刻度标签字号

        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.savefig(out_svg, dpi=300, bbox_inches='tight')
        plt.close()


if __name__=="__main__":
    main_path = "/home/gaoyun/poly/data/"
    out_species = "GB54"
    species = "Brettanomyces_bruxellensis"
    raw_long_reads_name = "ERR4624298"
    phase_out_name = "Phasing"

    max_ID = 0.05
    min_ovl = 0.3
    min_sim = 0.3
    min_overlap = 30
    min_ovl_len = 5
    min_len = 0
    # out_name = "phased"+"_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)

    phased_name = species+raw_long_reads_name
    out_dir = os.path.join(main_path, out_species)
    phased_dir = os.path.join(out_dir,phase_out_name,phased_name)
    phased_variants = os.path.join(phased_dir,"Variants/VCF")
    phased_reads_path = os.path.join(phased_dir, "Variants/Reads")
    phased_path = os.path.join(phased_dir,"Phased")
    raw_reads_path = os.path.join(out_dir, "rawData/longReads")
    mapped_long_path = os.path.join(out_dir, "Mapped/longReads")

    hetSNPs_positions_file = os.path.join(phased_variants,phased_name+".hetSNPs.positions.vcf.tsv")
    snp_reads_file = os.path.join(phased_reads_path, phased_name+".hetPositions.SNPxLongReads.vcf.chr.tsv")
    reference_file = os.path.join(out_dir,"reference",species,species+".fasta")
    ref_length_file = os.path.join(out_dir,"reference",species,species+"_chrLengths.txt")
    raw_long_reads_file = os.path.join(raw_reads_path,raw_long_reads_name+".fastq.gz")
    long_bam_file = os.path.join(mapped_long_path, raw_long_reads_name+".sorted.bam")

    out_path = os.path.join(phased_dir,"Analysed")
    out_reads_path = os.path.join(out_path,"reads_analysis")
    out_plot_path = os.path.join(phased_path,"Plots")
    all_paths = [out_path,out_plot_path, out_reads_path]
    for path in all_paths:
        os.makedirs(path,exist_ok=True)

    read_lengths_txt = os.path.join(out_reads_path,"read_lengths_txt")
    chrs = {"NC_054682.1","NC_054683.1","NC_054684.1","NC_054685.1","NC_054686.1","NC_054687.1","NC_054688.1","NC_054689.1","NC_054690.1"}

    het_file = os.path.join(out_path, "count_nHetSNPsPositions.tsv")
    intervals_length = 10000

    # 画直方图
    plot_chrs_intervals_snps(hetSNPs_positions_file, intervals_length,out_reads_path, raw_long_reads_name)

    snpReads_coverage_file = os.path.join(out_reads_path, raw_long_reads_name+"_snpReads_coverage.tsv")
    snpReads_coverage(snp_reads_file, intervals_length,snpReads_coverage_file)
    plot_bar_snpReads_coverage(snpReads_coverage_file,out_reads_path,raw_long_reads_name)

    snpReads_length_file = os.path.join(out_reads_path, raw_long_reads_name+"_snpReads_length.tsv")
    count_snpReads_length(snp_reads_file, intervals_length, snpReads_length_file)
    plot_bar_snpReads_lengths(snpReads_length_file,out_reads_path,raw_long_reads_name)