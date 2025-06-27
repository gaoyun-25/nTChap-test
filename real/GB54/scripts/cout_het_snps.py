import sys
import os
import bisect
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
import subprocess

# Distribution of chromosome length and heterozygous SNPs rate per chromosome of GB54 genome

# 计算chr02 每个chrs上call的 het snps个数，即每个基因的长度
# chr-n1-n2 n2-n1 number_het_SNPs
def count_chrs_het(hetSNPs_positions_file,ref_length_file, out_file):
    ref_length = dict()
    with open(ref_length_file) as f:
        for line in f:
            line = line.strip().split("\t")
            chr = line[0].split(" ")[0]
            ref_length[chr] = int(line[1])

    het_positions = dict()
    # vcf 1_based
    with open(hetSNPs_positions_file) as f:
        for line in f:
            line = line.strip().split("\t")
            chr = line[0]
            het_positions.setdefault(chr, list())
            het_positions[chr].append(int(line[1]))

    sort_chr_name = sorted(het_positions.keys())
    with open(out_file, "w") as g:
        g.write("chr-id\tchrs length (bp)\tnumber HetSNPs\tHetSNPs rates\n")
        for chr in sort_chr_name:
            positions = het_positions[chr]
            snp_rate = len(positions)/ref_length[chr]
            g.write(chr+"\t"+str(ref_length[chr])+"\t"+str(len(positions))+"\t"+str(snp_rate)+"\n")


def extract_chrs_het_by_id(chrs_het_file, chrs_id_file, out_file):
    chrs_id = set()
    with open(chrs_id_file) as f:
        for line in f:
            line=line.strip()
            chrs_id.add(line)

    with open(chrs_het_file) as f, open(out_file, "w") as g:
        first_line = next(f)
        g.write(first_line)
        for line in f:
            line = line.strip().split("\t")
            if line[0] in chrs_id:
                g.write("\t".join(line)+"\n")


def plot_hist_snps_nchrs(chrs_het_file,out_path):
    out_png=os.path.join(out_path,"snps_nchrs_hist.png")
    out_svg=os.path.join(out_path,"snps_nchrs_hist.svg")

    df=pd.read_csv(chrs_het_file,sep="\t")
    # 计算bins范围（从最小值到最大值，步长0.001）
    min_rate = df['HetSNPs rates'].min()
    max_rate = df['HetSNPs rates'].max()
    bins = np.arange(np.floor(min_rate * 1000)/1000,  # 向下取整到千分位
                    np.ceil(max_rate * 1000)/1000 + 0.001,  # 向上取整到千分位
                    0.001)
    plt.figure(figsize=(12.3, 7))  # 调整宽度为更合理的12.3英寸

    # 绘制直方图
    plt.hist(df['HetSNPs rates'], bins=bins, 
            edgecolor='black', alpha=0.7)

    # 设置坐标轴标签和标题
    plt.xlabel('HetSNPs Rates')
    plt.ylabel('Number of chrs')
    # plt.title('Distribution of HetSNPs Rates Across chrss')

    # 去除顶部和右侧边框
    ax = plt.gca()
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)

    # 关闭网格线
    ax.grid(False)

    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.savefig(out_svg, dpi=300, bbox_inches='tight')

def plot_dot_snps_chrsLength(chrs_het_file, out_path):
    out_png=os.path.join(out_path,"snps_chrsLength_dot.png")
    out_svg=os.path.join(out_path,"snps_chrsLength_dot.svg")
    # 读取数据（假设chrs_het_file已定义）
    df = pd.read_csv(chrs_het_file, sep="\t")

    plt.figure(figsize=(12, 7))

    # 绘制散点图，添加透明度、边缘颜色和点大小控制
    plt.scatter(x=df['HetSNPs rates'], 
                y=df['chrs length (bp)'],
                alpha=0.9,
                edgecolor='black',
                s=80,  # 点大小
                color='#426579')

    # 调整X轴和Y轴刻度标签大小
    plt.xticks(fontsize=14)  # X轴刻度标签字号
    plt.yticks(fontsize=14)  # Y轴刻度标签字号

    # 设置坐标轴标签和标题
    plt.xlabel('Heterozygous SNP rates', fontsize=18)
    plt.ylabel('Chromosome Length (bp)', fontsize=18)
    # plt.title('Relationship between HetSNPs Count and chrs Length', fontsize=14, pad=20)

    # 去除顶部和右侧边框
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # 添加浅灰色网格线
    # ax.grid(True, linestyle='--', alpha=0.3, color='grey')

    # # 优化坐标轴范围
    # plt.xlim(left=-5)  # 为左侧留出空隙
    # plt.ylim(bottom=-500)  # 为底部留出空隙

    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.savefig(out_svg, dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    main_path = "/home/gaoyun/poly/data"
    out_species = "GB54"
    species = "Brettanomyces_bruxellensis"
    phase_out_name = "Phasing"
    phased_name = "Brettanomyces_bruxellensisERR4624298"

    max_ID = 0.05
    min_ovl = 0.3
    min_sim = 0.3
    min_overlap = 30
    min_ovl_len = 5
    min_len = 0
    out_name = "phased"+"_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)

    out_dir = os.path.join(main_path, out_species)
    phased_dir = os.path.join(out_dir,phase_out_name,phased_name)
    phased_variants = os.path.join(phased_dir,"Variants/VCF")
    phased_path = os.path.join(phased_dir,"Phased")


    hetSNPs_positions_file = os.path.join(phased_variants,phased_name+".hetSNPs.positions.vcf.tsv")
    reference_file = os.path.join(out_dir,"reference",species,species+".fasta")
    ref_length_file = os.path.join(out_dir,"reference",species,species+"_chrLengths.txt")
    statistics_ref_len_command = "/home/gaoyun/poly/code/phasing/process_reads/GB54/statistics_chr_lengths.py"
    p = subprocess.run(["python",statistics_ref_len_command,reference_file,ref_length_file],stderr=subprocess.PIPE,stdout=subprocess.PIPE,universal_newlines=True)

    out_path = os.path.join(phased_dir,"Analysed")
    out_plot_path = os.path.join(phased_path,"Plots")
    all_paths = [out_path,out_plot_path]
    for path in all_paths:
        os.makedirs(path,exist_ok=True)

    het_file = os.path.join(out_path, "count_nHetSNPsPositions.tsv")
    count_chrs_het(hetSNPs_positions_file, ref_length_file,het_file)

    # chrs_id_file = "/home/gaoyun/poly/code/phasing/process_reads/CBS1483/dowloadData/out_chrs_id.txt"
    # chrs_id_het_file = os.path.join(out_path, "out_chrs_id_nHetSNPsPositions.tsv")
    # # extract_chrs_het_by_id(het_file, chrs_id_file, chrs_id_het_file)

    # plot_hist_snps_nchrs(het_file,out_path)
    plot_dot_snps_chrsLength(het_file, out_path)

    # plot_dot_snps_chrsLength(chrs_id_het_file, out_path)

    # consensusClusterPath = os.path.join(phased_path, out_name+"_variants.tsv")
    # readClusterPath = os.path.join(phased_path, out_name+"_clustered_read_name.tsv")
    # longReadFastQFilesPath = ""
    # plot_command = "/home/gaoyun/poly/code/phasing/table/analyze_GB54/plot.py"
    # p = subprocess.run(["python",plot_command,consensusClusterPath,readClusterPath,out_plot_path,longReadFastQFilesPath,out_name],stderr=subprocess.PIPE,stdout=subprocess.PIPE,universal_newlines=True)

