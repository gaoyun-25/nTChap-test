import os
import matplotlib
matplotlib.use('Agg')
from plotnine import *
import pysam
import random
import sys
import math
import subprocess

# gene_cid file -- clusters are present in each gene
# tag bam with XC=cluster id

def load_consensus_cluster(file):
    clusters_consensus = dict()
    with open(file) as f:
        for line in f:
            line = line.strip("\n").split("\t")
            clusters_consensus.setdefault(line[0], set())
            snp = line[1]+":"+line[2]+"="+line[3]
            clusters_consensus[line[0]].add(snp)

    return clusters_consensus


def getCidGeneVariants(clustered_file,gene_cid_file,clustered_gene_file,clustered_reads_file,clustered_reads_gene_file,out_gene_id):
    with open(clustered_reads_file) as f:
        with open(clustered_reads_gene_file,"w") as g:
            for line in f:
                g.write(line)

    gene_clusters = dict()
    clusters_consensus = load_consensus_cluster(clustered_file)
    with open(gene_cid_file,"r") as f:
        for line in f:
            line = line.strip().split("\t")
            gene_clusters.setdefault(line[5],{"length":0, "clusters":set()})
            gene_clusters[line[5]]["length"] = int(line[4])
            gene_clusters[line[5]]["clusters"].add(line[0])

    with open(clustered_gene_file,"w") as f:
        for gene_id in out_gene_id:
            if gene_id in gene_clusters.keys():
                for cid in gene_clusters[gene_id]["clusters"]:
                    for snp in clusters_consensus[cid]:
                        pos_base = snp.split(":")[1].split("=")
                        f.write(cid+"\t"+gene_id+"\t"+pos_base[0]+"\t"+pos_base[1]+"\n")

def getCidGene(clustered_file,bed_file, chrom, gene_cid_file):
    # 一个cluster可能属于多个gene regions
    # cid chr from to gene(;)
    cid_positions = dict()
    with open(clustered_file) as f:
        for line in f:
            line = line.strip("\n").split("\t")
            cid_chr = line[0]+"+"+line[1]
            cid_positions.setdefault(cid_chr,set())
            cid_positions[cid_chr].add(int(line[2]))

    cid_positions_from_to = dict()
    # with open(out_file, "w") as g:
    for cid_chr, positions in cid_positions.items():
        min_pos = min(positions)
        max_pos = max(positions)
            # g.write(cid_chr.split("+")[0]+"\t"+cid_chr.split("+")[1]+"\t"+str(min_pos)+"\t"+str(max_pos)+"\n")
        cid_positions_from_to[cid_chr] = (min_pos,max_pos)
    print("cid_positions_from_to", len(cid_positions_from_to))

    # load bed file
    gene_regions = dict()
    with open(bed_file) as f:
        for line in f:
            line = line.strip().split("\t")
            chr = line[0]
            if chr == chrom:
                gene_regions.setdefault(chr, set())
                gene_regions[chr].add((int(line[1]),int(line[2])))

    # 比较gene区域和cid的区域是否有重叠的部分，找出gene区域有哪些cid cov
    gene_cid = dict()
    for gene_interval in gene_regions[chrom]:
        gene_region_string = chrom+"-"+str(gene_interval[0])+"-"+str(gene_interval[1])
        interval_length = gene_interval[1] - gene_interval[0] + 1
        gene_cid.setdefault(gene_region_string,set())
        for cid_chr, position_interval in cid_positions_from_to.items():
            start_interval = max(gene_interval[0], position_interval[0])
            end_interval = min(gene_interval[1], position_interval[1])
            if start_interval <= end_interval:
                # 有交集
                gene_cid[gene_region_string].add((cid_chr.split("+")[0],cid_chr.split("+")[1],str(position_interval[0]),str(position_interval[1]),str(interval_length)))

    with open(gene_cid_file, "w") as f:
        for gene, cid_sets in gene_cid.items():
            for cid_tuple in cid_sets:
                f.write("\t".join(cid_tuple)+"\t"+gene+"\n")


def loadClusteredReads(clustered_reads_file):
    clustered_reads = dict()
    with open(clustered_reads_file) as f:
        for line in f:
            line = line.strip().split("\t")
            clustered_reads.setdefault(line[0],set())
            clustered_reads[line[0]].add(line[1])

    return clustered_reads


def loadBamFile(bam_file):
    alignments = dict()
    if bam_file.endswith("bam"):
        with pysam.AlignmentFile(bam_file,"rb") as f:
            header_dict = f.header.to_dict()
            for line in f:
                read_name = line.query_name
                alignments.setdefault(read_name, set())
                alignments[read_name].add(line)
    elif bam_file.endswith("sam"):
        with pysam.AlignmentFile(bam_file,"r") as f:
            header_dict = f.header.to_dict()
            for line in f:
                read_name = line.query_name
                alignments.setdefault(read_name, set())
                alignments[read_name].add(line)

    return alignments,header_dict


def generate_random_color():
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))


def tagBam(bam_file, gene_cid_file, clustered_reads_file,num_top,out_path, out_name):
    # bed文件中选择chr02 gene中长度前num_top的regions
    # 输出前num_top gene的bam，每个bam中加tag，属于哪个cid
    split_bam_path = os.path.join(out_path,"Bam")
    os.makedirs(split_bam_path,exist_ok=True)

    clustered_reads = loadClusteredReads(clustered_reads_file)

    gene_clusters_info = dict()
    with open(gene_cid_file,"r") as f:
        for line in f:
            line = line.strip().split("\t")
            gene_clusters_info.setdefault(line[5],{"length":0, "clusters":dict()})
            gene_clusters_info[line[5]]["length"] = int(line[4])
            gene_clusters_info[line[5]]["clusters"][line[0]] = clustered_reads[line[0]]

    # sort gene length
    sorted_gene = sorted(gene_clusters_info, key=lambda k:gene_clusters_info[k]["length"], reverse=True)
    i = 1
    # 打开bam
    alignments,header_dict = loadBamFile(bam_file)

    for gene_id in sorted_gene:
        if i > num_top:
            break
        out_gene_bam = os.path.join(split_bam_path,out_name+"_"+gene_id+".bam")
        used_colors = set()
        with pysam.AlignmentFile(out_gene_bam,"wb",header=header_dict) as f:
            for cid, reads_name in gene_clusters_info[gene_id]["clusters"].items():
                cid_algns = list()
                color = generate_random_color()
                while True:
                    if color not in used_colors:
                        used_colors.add(color)
                        break
                    else:
                        color = generate_random_color()
                for read_name in reads_name:
                    for algn in alignments[read_name]:
                        algn.set_tag("XC",cid,value_type="Z")
                        algn.set_tag("YC","175,238,238",value_type="Z")
                        cid_algns.append(algn)
                        f.write(algn)
        i += 1
        out_gene_sort_bam = os.path.join(split_bam_path,out_name+"_"+gene_id+".sort.bam")
        pysam.sort("-o",out_gene_sort_bam,out_gene_bam)
        pysam.index(out_gene_sort_bam)

    return None


def given_gene_id(file):
    gene_id = list()
    with open(file) as f:
        for line in f:
            line = line.strip()
            gene_id.append(line)

    return gene_id


if __name__ == "__main__":

    max_ID = 0.05
    min_ovl = 0.3
    min_sim = 0.3
    min_overlap = 30
    min_ovl_len = 5
    min_len = 0
    chrom = "chr02"
    out_name = "phased_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap)+"_"+str(min_ovl_len)
    num_top = 20

    # outdir = sys.argv[1]
    # phased_dir = sys.argv[2]
    # bed_file = sys.argv[3]
    # bam_file = sys.argv[4]
    # chrom = sys.argv[5]
    # num_top = int(sys.argv[6])
    # out_name = sys.argv[7]

    outdir = "/home/gaoyun/poly/data/potato/Phasing_gene_chr_nx_vcf/potatolongReadsOnt.chr02/Phased"
    phased_dir = outdir
    bed_file = "/home/gaoyun/poly/data/potato/reference/potato/annotation/all_genes.bed"
    bam_file = "/home/gaoyun/poly/data/potato/Mapped/longReads/longReadsOnt.chr02.sorted.sam"

    out_path = os.path.join(outdir,"PlotsGene/")
    all_path = [out_path]
    for path in all_path:
        os.makedirs(path, exist_ok=True)

    clustered_file = os.path.join(phased_dir, out_name+"_variants.tsv")
    clustered_reads_file = os.path.join(phased_dir, out_name+"_clustered_read_name.tsv")
    gene_cid_file = os.path.join(out_path, out_name+"_blocks_gene.tsv")
    getCidGene(clustered_file,bed_file, chrom, gene_cid_file)

    clustered_gene_file = os.path.join(out_path,out_name+"_given_gene_variants.tsv")
    clustered_reads_gene_file = os.path.join(out_path, out_name+"_given_gene_clustered_read_name.tsv")
    given_cid_txt = "/home/gaoyun/poly/code/phasing/final_nTChap/process_data/real/potato/scripts/out_gene_id.txt"
    out_gene_id = given_gene_id(given_cid_txt)
    getCidGeneVariants(clustered_file,gene_cid_file,clustered_gene_file,clustered_reads_file,clustered_reads_gene_file,out_gene_id)

    tagBam(bam_file, gene_cid_file, clustered_reads_gene_file,num_top,out_path,out_name)

    command_py = "/home/gaoyun/poly/code/phasing/final_nTChap/process_data/real/potato/scripts/plot.py"
    longReadFastQFilesPath = ""
    outName = out_name+"_given"
    p=subprocess.run(["python3", command_py, clustered_gene_file, clustered_reads_gene_file,out_path, longReadFastQFilesPath,outName],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)

