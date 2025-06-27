import sys
import os
import bisect
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
import subprocess

# Distribution of phased heterozygous SNPs number covered by clusters

def hetSNPs_nclusters(consensus_snps_file,het_cid_file,het_ncid_file,het_positions_tsv):
    clusters = set()
    het_cid = dict()
    with open(consensus_snps_file) as f:
        for line in f:
            line = line.strip().split("\t")
            clusters.add(line[0])
            chr_pos = line[1]+":"+line[2]
            het_cid.setdefault(chr_pos,set())
            het_cid[chr_pos].add(line[0])

    het_positions = set()
    with open(het_positions_tsv) as f:
        for line in f:
            line = line.strip().split("\t")
            chr_pos = line[0]+":"+line[1]
            het_positions.add(chr_pos)

    phased_het_positions = set(het_cid.keys())
    unphased_het_positions = het_positions - phased_het_positions

    with open(het_cid_file, "w") as g:
        g.write("unphased pos: "+str(len(unphased_het_positions))+"\t"+",".join(unphased_het_positions)+"\n")
        g.write("chr:pos\tncids\tcids\n")
        for pos,cids in het_cid.items():
            g.write(pos+"\t"+str(len(cids))+"\t"+",".join(cids)+"\n")

    het_ncid = dict()
    for pos,cids in het_cid.items():
        ncid = len(cids)
        het_ncid.setdefault(ncid,0)
        het_ncid[ncid] += 1
    total_pos = len(unphased_het_positions)
    ncid234 = 0
    n = {2,3,4}
    n3 = {3}
    ncid3 = 0 
    with open(het_ncid_file, "w") as g:
        g.write("Covered Clusters Number\tHetSNPs position Number\n")
        for ncids, npos in het_ncid.items():
            g.write(str(ncids)+"\t"+str(npos)+"\n")
            total_pos += npos
            if ncids in n:
                ncid234 += npos
            if ncids in n3:
                ncid3 = npos
        g.write(str(0)+"\t"+str(len(unphased_het_positions))+"\n")

    with open(note, "a") as f:
        f.write(consensus_snps_file+"\n")
        f.write("total pos: "+str(total_pos)+"\n")
        f.write("unphased_het_positions: "+str(len(unphased_het_positions))+"\n")
        f.write("unphased rate(%) "+str((len(unphased_het_positions)/total_pos)*100)+"\n")
        f.write("phased_rates234(%): "+str((ncid234/total_pos)*100)+"\n")
        f.write("phased_rates3(%): "+str((ncid3/total_pos)*100)+"\n")
        f.write("cids numbers per chrs(len(clusters)/9): "+str(len(clusters)/9)+"\n")
        f.write("cids numbers per haplotype(len(clusters)/27): "+str(len(clusters)/27)+"\n")


def plot_hist_snps_nclusters(het_ncid_file,out_path,out_name):
    out_png=os.path.join(out_path,out_name+"_snps_nclusters_hist.png")
    out_svg=os.path.join(out_path,out_name+"_snps_nclusters_hist.svg")

    df=pd.read_csv(het_ncid_file,sep="\t")

    # 提取数据列（假设列名为 'x' 和 'y'，按实际列名修改）
    x = df['Covered Clusters Number'].tolist()
    y = df['HetSNPs position Number'].tolist()

    # 创建画布和坐标轴
    fig, ax = plt.subplots(figsize=(10, 6))

    # 绘制条形图（宽度为1，对齐到刻度边缘）
    bars = ax.bar(
        x, 
        y, 
        width=0.6, 
        # align='edge',  # 条形左侧对齐刻度线
        edgecolor='black', 
        color='#426579',
        # linewidth=0.5
    )
    # 在条形顶部标注数值
    for bar in bars:
        height = bar.get_height()
        ax.text(
            x=bar.get_x() + bar.get_width() / 2,  # 水平居中
            y=height + 0.3,                       # 垂直位置（数值上方）
            s=f'{int(height)}',                   # 显示整数
            ha='center', va='bottom',             # 对齐方式
            fontsize=14
        )

    ax.set_xticks(x)  # 确定刻度位置
    ax.set_xticklabels([str(i) for i in x], fontsize=14)
    # 设置Y轴刻度和标签 - 增大字体
    ax.tick_params(axis='y', labelsize=14)  # 设置Y轴刻度标签大小

    # 添加标签和标题
    ax.set_xlabel('Covered Clusters Number per Position', fontsize=18)
    ax.set_ylabel('Number of Heterozygous SNP position', fontsize=18)

    # 关闭网格线
    ax.grid(False)

    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.savefig(out_svg, dpi=300, bbox_inches='tight')


def count_hetSNPs_per_chrs(het_cid_file,out_file):
    chrs_hetPos = dict()
    chrs_hetPos["unphased"] = dict()
    chrs_hetPos["phased"] = dict()
    with open(het_cid_file) as f:
        for line in f:
            line = line.strip().split("\t")
            if "unphased pos" in line[0]:
                unphased_pos = line[1].split(",")
                for chr_pos in unphased_pos:
                    chr_pos_split = chr_pos.split(":")
                    chrs_hetPos["unphased"].setdefault(chr_pos_split[0], set())
                    chrs_hetPos["unphased"][chr_pos_split[0]].add(chr_pos_split[1])
            elif "chr:pos" in line[0]:
                continue
            else:
                chr_pos_split = line[0].split(":")
                chrs_hetPos["phased"].setdefault(chr_pos_split[0], set())
                chrs_hetPos["phased"][chr_pos_split[0]].add(chr_pos_split[1])

    with open(out_file,"w") as g:
        g.write("unphased\n")
        unphased_chr_pos = sorted(chrs_hetPos["unphased"].keys())
        for chr in unphased_chr_pos:
            positions = chrs_hetPos["unphased"][chr]
            g.write(chr+"\t"+str(len(positions))+"\n")
        g.write("phased\n")
        for chr, positions in chrs_hetPos["phased"].items():
            g.write(chr+"\t"+str(len(positions))+"\n")

def count_chr_nclusters(cosnensus_path, note_txt):
    chr_clusters = dict()
    with open(cosnensus_path) as f:
        for line in f:
            line = line.strip().split("\t")
            chr = line[1]
            chr_clusters.setdefault(chr, set())
            chr_clusters[chr].add(line[0])
    
    with open(note_txt, "a") as g:
        g.write(cosnensus_path+"\n")
        for chr, cids in chr_clusters.items():
            g.write(chr+"\t"+str(len(cids))+"\n")

if __name__ == "__main__":
    main_path = "/home/gaoyun/poly/data/"
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

    out_path = os.path.join(phased_dir,"Analysed")
    out_plot_path = os.path.join(phased_path,"Plots")
    all_paths = [out_path,out_plot_path]
    for path in all_paths:
        os.makedirs(path,exist_ok=True)

    consensusClusterPath = os.path.join(phased_path, out_name+"_variants.tsv")
    readClusterPath = os.path.join(phased_path, out_name+"_clustered_read_name.tsv")

    note = os.path.join(out_path,"note.txt")
    het_cid_file = os.path.join(out_path,"hetSNPs_positions_covered_cid.tsv")
    het_ncid_file = os.path.join(out_path,"hetSNPs_positions_number_of_covered_cid.tsv")
    # hetSNPs_nclusters(consensusClusterPath,het_cid_file,het_ncid_file,hetSNPs_positions_file)
    # out_hist_name = "phased"
    # plot_hist_snps_nclusters(het_ncid_file,out_path,out_hist_name)
    # out_chr_phased_pos_file = os.path.join(out_path,"chr_phased_positions.tsv")
    # count_hetSNPs_per_chrs(het_cid_file,out_chr_phased_pos_file)
    # # count_chr_nclusters(consensusClusterPath, note)

    clean_ways = ["dup"]
    out_clean_name = "Cleaned_chr_f_d_99"
    clean_path = os.path.join(phased_path, out_clean_name)
    for clean_way in clean_ways:
        clean_consensusClusterPath = os.path.join(clean_path, out_name+"_cleaned_"+clean_way+"_variants.tsv")
        clean_readClusterPath = os.path.join(clean_path, out_name+"_cleaned_"+clean_way+"_clustered_read_name.tsv")
        het_cid_file = os.path.join(out_path,out_name+"_cleaned_"+clean_way+"_"+out_clean_name+"_hetSNPs_positions_covered_cid.tsv")
        het_ncid_file = os.path.join(out_path,out_name+"_cleaned_"+clean_way+"_"+out_clean_name+"_hetSNPs_positions_number_of_covered_cid.tsv")
        hetSNPs_nclusters(clean_consensusClusterPath,het_cid_file,het_ncid_file,hetSNPs_positions_file)
        out_hist_name = out_name+"_cleaned_"+clean_way+"_"+out_clean_name
        plot_hist_snps_nclusters(het_ncid_file,out_path,out_hist_name)
        out_chr_phased_pos_file = os.path.join(out_path,out_name+"_cleaned_"+clean_way+"_"+out_clean_name+"_chr_phased_positions.tsv")
        count_hetSNPs_per_chrs(het_cid_file,out_chr_phased_pos_file)
        count_chr_nclusters(clean_consensusClusterPath, note)