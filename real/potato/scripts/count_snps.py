import sys
import os
import bisect
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np

# fig a, fig b

# 计算chr02 每个gene上call的 het snps个数，即每个基因的长度
# chr-n1-n2 n2-n1 number_het_SNPs
def count_gene_het(hetSNPs_positions_file, gene_bed,out_file,chrs):
    gene_intervals = dict()
    with open(gene_bed) as f:
        for line in f:
            line = line.strip().split("\t")
            chr = line[0]
            if chr in chrs:
                gene_intervals.setdefault(chr,list())
                gene_intervals[chr].append([int(line[1]), int(line[2])]) #bed都是[start,end) bed 0-based; 下面计数把与右端点相同值的点也计算上了

    het_positions = dict()
    # vcf 1_based
    with open(hetSNPs_positions_file) as f:
        for line in f:
            line = line.strip().split("\t")
            chr = line[0]
            if chr in chrs:
                het_positions.setdefault(chr, list())
                het_positions[chr].append(int(line[1])-1)

    with open(out_file, "w") as g:
        g.write("chr-gene-id\tgene length (bp)\tnumber HetSNPs\tHetSNPs rates\n")
        for chr in chrs:
            intervals = gene_intervals[chr]
            positions = het_positions[chr]
            intervals_postions = count_per_interval(intervals, positions)
            for i,interval in enumerate(intervals):
                len_interval = interval[1] - interval[0]+1
                snp_rate = intervals_postions[i]/len_interval
                g.write(chr+"-"+"-".join([str(x) for x in interval])+"\t"+str(len_interval)+"\t"+str(intervals_postions[i])+"\t"+str(snp_rate)+"\n")


def count_per_interval(intervals, numbers):
    # 预处理：将数排序以支持二分查找
    sorted_numbers = sorted(numbers)
    result = list()
    for interval in intervals:
        # 确保区间左右端点有序
        a, b = interval
        left = min(a, b)
        right = max(a, b)
        # 查找第一个 >= left 的位置
        l_pos = bisect.bisect_left(sorted_numbers, left)
        # 查找第一个 > right 的位置
        r_pos = bisect.bisect_right(sorted_numbers, right)
        # 计算落在区间内的数的数量
        count = r_pos - l_pos
        result.append(count)
    return result


def extract_gene_het_by_id(gene_het_file, gene_id_file, out_file):
    gene_id = set()
    with open(gene_id_file) as f:
        for line in f:
            line=line.strip()
            gene_id.add(line)

    with open(gene_het_file) as f, open(out_file, "w") as g:
        first_line = next(f)
        g.write(first_line)
        for line in f:
            line = line.strip().split("\t")
            if line[0] in gene_id:
                g.write("\t".join(line)+"\n")


def plot_hist_snps_nGene(gene_het_file,out_path):
    out_png=os.path.join(out_path,"snps_nGene_hist.png")
    out_svg=os.path.join(out_path,"snps_nGene_hist.svg")

    df=pd.read_csv(gene_het_file,sep="\t")
    # 计算bins范围（从最小值到最大值，步长0.001）
    min_rate = df['HetSNPs rates'].min()
    max_rate = df['HetSNPs rates'].max()
    bins = np.arange(np.floor(min_rate * 1000)/1000,  # 向下取整到千分位
                    np.ceil(max_rate * 1000)/1000 + 0.001,  # 向上取整到千分位
                    0.001)
    plt.figure(figsize=(12.3, 7))  # 调整宽度为更合理的12.3英寸

    # 绘制直方图
    plt.hist(df['HetSNPs rates'], bins=bins, 
            edgecolor='black',color='#426579')

    # 设置坐标轴标签和标题
    plt.xlabel('Heterozygous SNPs Rates', fontsize=18)
    plt.ylabel('Number of Genes', fontsize=18)
    # plt.title('Distribution of HetSNPs Rates Across Genes')
    # 调整X轴和Y轴刻度标签大小
    plt.xticks(fontsize=14)  # X轴刻度标签字号
    plt.yticks(fontsize=14)  # Y轴刻度标签字号
    # 去除顶部和右侧边框
    ax = plt.gca()
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)

    # 关闭网格线
    ax.grid(False)

    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.savefig(out_svg, dpi=300, bbox_inches='tight')

def plot_dot_snps_geneLength(gene_het_file, out_path):
    out_png=os.path.join(out_path,"snps_geneLength_dot.png")
    out_svg=os.path.join(out_path,"snps_geneLength_dot.svg")
    # 读取数据（假设gene_het_file已定义）
    df = pd.read_csv(gene_het_file, sep="\t")

    plt.figure(figsize=(12, 7))

    # 绘制散点图，添加透明度、边缘颜色和点大小控制
    plt.scatter(x=df['number HetSNPs'], 
                y=df['gene length (bp)'],
                alpha=0.9,
                edgecolor='black',
                s=80,  # 点大小
                color='#426579')

    # 设置坐标轴标签和标题
    plt.xlabel('Number of Heterozygous SNPs', fontsize=18)
    plt.ylabel('Gene Length (bp)', fontsize=18)
    # plt.title('Relationship between HetSNPs Count and Gene Length', fontsize=14, pad=20)
    # 调整X轴和Y轴刻度标签大小
    plt.xticks(fontsize=14)  # X轴刻度标签字号
    plt.yticks(fontsize=14)  # Y轴刻度标签字号
    # 去除顶部和右侧边框
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # 添加浅灰色网格线
    # ax.grid(True, linestyle='--', alpha=0.3, color='grey')

    # 优化坐标轴范围
    plt.xlim(left=-5)  # 为左侧留出空隙
    plt.ylim(bottom=-500)  # 为底部留出空隙

    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.savefig(out_svg, dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    main_path = "/home/gaoyun/poly/data/"
    species = "potato"
    bed_name = "all_genes.chr02.bed"
    phase_out_name = "Phasing_gene"
    phased_name = "potatolongReadsOnt.chr02"
    chrs = {"chr02"}

    out_dir = os.path.join(main_path, species)
    gene_bed = os.path.join(out_dir, "reference/potato/annotation",bed_name)
    phased_dir = os.path.join(out_dir,phase_out_name,phased_name)
    phased_variants = os.path.join(phased_dir,"Variants/VCF")
    hetSNPs_positions_file = os.path.join(phased_variants,phased_name+".hetSNPs.positions.vcf.tsv")

    out_path = os.path.join(phased_dir,"Analysed")
    os.makedirs(out_path,exist_ok=True)
    gene_het_file = os.path.join(out_path, "count_gene_nHetSNPsPositions.tsv")

    count_gene_het(hetSNPs_positions_file, gene_bed,gene_het_file,chrs)

    gene_id_file = "/home/gaoyun/poly/code/phasing/final_nTChap/process_data/real/potato/scripts/out_gene_id.txt"
    gene_id_het_file = os.path.join(out_path, "out_gene_id_nHetSNPsPositions.tsv")
    extract_gene_het_by_id(gene_het_file, gene_id_file, gene_id_het_file)

    plot_hist_snps_nGene(gene_het_file,out_path)

    plot_dot_snps_geneLength(gene_id_het_file, out_path)