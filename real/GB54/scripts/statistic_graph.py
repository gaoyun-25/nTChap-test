import sys
import os
from itertools import combinations
import multiprocessing
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# 计算加权全局稠密度
def weighted_global_density(G):
    n = G.number_of_nodes()
    if n < 2:
        return 0.0  # 避免除以零
    
    # 最大可能边数（无向图）
    max_possible_edges = n * (n - 1) / 2
    
    # 考虑所有可能的节点对
    total_weight = 0.0
    nodes = list(G.nodes)
    
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            u, v = nodes[i], nodes[j]
            # 如果边存在，获取权重；否则权重为0
            if G.has_edge(u, v):
                total_weight += G[u][v]['weight']
    
    density = total_weight / max_possible_edges
    return density

# 计算权重方差
def weight_variance(G):
    # 获取所有实际存在的边的权重
    weights = [data['weight'] for u, v, data in G.edges(data=True)]
    
    if len(weights) < 2:
        return 0.0  # 方差需要至少2个值
    
    return np.var(weights)

# 计算平均加权聚集系数
def average_weighted_clustering_coefficient(G):
    # 计算每个节点的加权聚集系数
    clustering_coeffs = []
    
    for node in G.nodes:
        neighbors = list(G.neighbors(node))
        k = len(neighbors)
        
        if k < 2:
            clustering_coeffs.append(0.0)
            continue
        
        # 计算三角形可能数
        triangle_sum = 0.0
        possible_triangles = 0
        
        # 检查所有邻居对
        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                u, v = neighbors[i], neighbors[j]
                
                # 只有u和v相连时才形成三角形
                if G.has_edge(u, v):
                    w_ij = G[u][v]['weight']
                    w_iu = G[node][u]['weight']
                    w_iv = G[node][v]['weight']
                    
                    # 使用几何平均值作为权重测量
                    triangle_weight = (w_iu * w_iv * w_ij) ** (1/3)
                    triangle_sum += triangle_weight
                
                possible_triangles += 1
        
        # 如果没有任何可能的三角形，系数为0
        if possible_triangles == 0:
            clustering_coeffs.append(0.0)
        else:
            clustering_coeffs.append(triangle_sum / possible_triangles)
    
    return np.mean(clustering_coeffs)

# 计算平均节点强度
def average_node_strength(G):
    strengths = []
    
    for node in G.nodes:
        # 计算节点的强度（所有相邻边权重之和）
        strength = sum(G[node][neighbor]['weight'] for neighbor in G.neighbors(node))
        strengths.append(strength)
    
    return np.mean(strengths)


def load_reduced_reads(reduced_read_file):
    # load reduced reads
    reads_dict = dict()
    with open(reduced_read_file, "r") as f:
        for line in f:
            if line[0] == ">":
                line=line.strip()
                ctg_name = line[1:]
                reads_dict[ctg_name] = dict()
            else:
                line = line.strip("\n").split("\t")
                reads_dict[ctg_name][line[0]] = frozenset(line[1:])

    return reads_dict


# calculate similarity between reads
def get_similarity(common_pos, consensusA, consensusB, min_ovl, min_sim, min_overlap, min_ovl_len):
    # common_pos = self.positions & other.positions
    if len(consensusA) < len(consensusB):
        setA = consensusA
        setB = consensusB
    else:
        setA = consensusB
        setB = consensusA
    if len(common_pos)/len(setA) >= min_ovl and len(common_pos) >= min_ovl_len or len(common_pos) >= min_overlap:
        common_set = setA & setB
        local_similarity = len(common_set)/len(common_pos)
        if local_similarity < min_sim:
            local_similarity = 0
    else:
        local_similarity = 0

    return local_similarity

# initialize G and cluster_info
def initializeG(queue, out_queue, min_ovl, min_sim, min_overlap, min_ovl_len):
    local_infos = dict()
    while True:
        combination_batch = queue.get()
        if combination_batch == "end":
            out_queue.put(local_infos)
            return None
        for comb in combination_batch:
            v_i = comb[0]
            v_j = comb[1]
            cluster_i = clusters_infos[v_i]
            cluster_j = clusters_infos[v_j]
            local_infos.setdefault(v_i, {})
            local_infos.setdefault(v_j, {})
            common_positions = cluster_i["positions"] & cluster_j["positions"]
            if len(common_positions) >= min_ovl_len:
                similarity = get_similarity(common_positions, cluster_i["consensus"], cluster_j["consensus"], min_ovl, min_sim, min_overlap, min_ovl_len)
                if similarity == 0:
                    continue
                local_infos[v_i][v_j] = similarity
                local_infos[v_j][v_i] = similarity


def start_initializeG(queue, out_queue, clusters_infos_queue, clusters_infos, min_ovl, min_sim, min_overlap, min_ovl_len, thread_number):
    all_process = list()
    for thread in range(0, thread_number-1):
        initializeG_p = multiprocessing.Process(target=initializeG, args=((queue),(out_queue),(min_ovl),(min_sim),(min_overlap),(min_ovl_len)))
        initializeG_p.daemon = True
        initializeG_p.start()
        all_process.append(initializeG_p)
    initialize_merge_p = multiprocessing.Process(target=initialize_merge, args=((out_queue),(clusters_infos_queue),(clusters_infos),))
    initialize_merge_p.daemon = True
    initialize_merge_p.start()
    all_process.append(initialize_merge_p)

    return all_process


def initialize_merge(out_queue, clusters_infos_queue, clusters_infos):
    while True:
        local_dict = out_queue.get()
        if local_dict == "end":
            clusters_infos_queue.put(clusters_infos)
            return None
        for nameA in local_dict.keys():
            for nameB, similarity in local_dict[nameA].items():
                clusters_infos[nameA]["similarities"][nameB] = {'weight':similarity}

def compute_reads_similarity(reduced_reads, min_ovl, min_sim, min_overlap, min_ovl_len, thread_number):
    global clusters_infos
    clusters_infos = dict()
    for read_name, var_seq in reduced_reads.items():
        positions = set([var.split("=")[0] for var in var_seq])
        clusters_infos[read_name] = {"positions": positions, "consensus": var_seq, "similarities": dict()}

    print("compute similarity between reduced reads")
    batch_size = 50000

    #Create a queue for all the combinations
    combination_queue = multiprocessing.Queue(maxsize=3*thread_number)
    local_similarity_dict_queue = multiprocessing.Queue()
    clusters_infos_queue = multiprocessing.Queue()

    # Start the cache filling processes up before filling up the queue
    all_process = start_initializeG(combination_queue, local_similarity_dict_queue, clusters_infos_queue, clusters_infos, min_ovl, min_sim, min_overlap, min_ovl_len, thread_number)
    all_initG = all_process[0:-1]
    init_merge = all_process[-1]

    # Fill up the queue with combinations
    current_batch = []
    i = 0
    for comb in combinations(clusters_infos.keys(), 2):
        current_batch.append(comb)
        if len(current_batch) > batch_size:
            combination_queue.put(current_batch)
            i += 1
            current_batch = []
    #  Ensure that when current_batch is less than 50000, it can still be successfully added.
    if len(current_batch) > 0:
        combination_queue.put(current_batch)
        current_batch=[]

    # Fill the queue with end signals
    for thread in range(0, thread_number):
        combination_queue.put("end")

    for initG_p in all_initG:
        initG_p.join()

    local_similarity_dict_queue.put("end")

    clusters_infos = clusters_infos_queue.get()

    init_merge.join()

    G = nx.Graph()
    for cidA, cluster in clusters_infos.items():
        for cidB, similarity in cluster["similarities"].items():
            G.add_edge(cidA, cidB, weight=similarity['weight'])
    del clusters_infos

    return G

def main(reads_dict, out_path, out_name, thread_number,min_ovl, min_sim, min_overlap, min_ovl_len):
    statistical_magnitude_file = os.path.join(out_path, out_name+"_"+"statistical_magnitude.tsv")

    statistical_magnitude_text = "chromosomes\tstatistical_magnitude\tvalue\n"
    for ctg_name, reduced_reads in reads_dict.items():
        G = compute_reads_similarity(reduced_reads, min_ovl, min_sim, min_overlap, min_ovl_len, thread_number)

        # 计算并显示统计量
        density = weighted_global_density(G)
        variance = weight_variance(G)
        clustering = average_weighted_clustering_coefficient(G)
        strength = average_node_strength(G)
        node_count = G.number_of_nodes()
        edge_count = G.number_of_edges()
        statistical_magnitude_text += ctg_name+"\t"+"Global Weighted Density\t"+str(density)+"\n"
        statistical_magnitude_text += ctg_name+"\t"+"Weight Variance\t"+str(variance)+"\n"
        statistical_magnitude_text += ctg_name+"\t"+"Avg. Weighted Clustering Coefficient\t"+str(clustering)+"\n"
        statistical_magnitude_text += ctg_name+"\t"+"Avg. Node Strength\t"+str(strength)+"\n"
        statistical_magnitude_text += ctg_name+"\t"+"Number of Nodes\t"+str(node_count)+"\n"
        statistical_magnitude_text += ctg_name+"\t"+"Number of Edges\t"+str(edge_count)+"\n"

    with open(statistical_magnitude_file,"w") as f:
        f.write(statistical_magnitude_text)
    
    # # 添加解释
    # print("\n结果解释:")
    # print(f"- 加权全局稠密度({density:.4f})接近0表示整体相似度低，接近1表示整体相似度高")
    # print(f"- 权重方差({variance:.4f})高暗示权重分布不均匀，可能表示存在社区结构")
    # print(f"- 平均加权聚集系数({clustering:.4f})高表示图的局部聚类紧密")
    # print(f"- 平均节点强度({strength:.4f})表示每个节点的平均加权连接度")

def plot_statistical_magnitude(statistical_magnitude_file, out_path, out_name):
    # 读取数据文件
    data = pd.read_csv(statistical_magnitude_file, sep='\t')

    statistical_magnitude_png = os.path.join(out_path, out_name+"_statistical_magnitude.png")
    statistical_magnitude_svg = os.path.join(out_path, out_name+"_statistical_magnitude.svg")

    # 确保染色体顺序保留文件中的原始顺序
    chromosomes_ordered = data['chromosomes'].unique()

    # 数据重组：转为宽格式，每一行对应一个染色体
    pivot_df = data.pivot(index='chromosomes', columns='statistical_magnitude', values='value')

    # 创建组合标签：将"Number of Nodes"和"Number of Edges"合并为注释
    # 注意：确保这些统计量在文件中存在
    pivot_df['node_edge_label'] = pivot_df.apply(
        lambda row: f"{row.name}\n({int(row['Number of Nodes'])}, {int(row['Number of Edges'])})", 
        axis=1
    )

    # 设置图形和双y轴
    fig, ax_left = plt.subplots(figsize=(14, 7))
    ax_right = ax_left.twinx()

    # 左侧y轴的统计量（0-1范围）
    left_stats = [
        'Global Weighted Density',
        'Weight Variance',
        'Avg. Weighted Clustering Coefficient'
    ]
    # 右侧y轴的统计量
    right_stats = ['Avg. Node Strength']

    # 绘制左侧y轴的统计量
    for stat in left_stats:
        ax_left.plot(pivot_df.index, pivot_df[stat], marker='o', label=stat)

    # 绘制右侧y轴的统计量
    for stat in right_stats:
        ax_right.plot(pivot_df.index, pivot_df[stat], marker='s', linestyle='--', label=stat, color='purple')

    # 设置坐标轴标签
    ax_left.set_xlabel('Chromosomes and (number of nodes, number of edges)',fontsize=14)
    ax_left.set_ylabel('Statistical Magnitudes (0-1 scale)', color="#083046", fontsize=14)
    ax_right.set_ylabel('Avg. Node Strength (60-80 scale)', color='purple', fontsize=14)

    # 设置y轴范围
    ax_left.set_ylim(0, 1)
    ax_right.set_ylim(60, 90)

    # 设置坐标轴刻度文字大小
    ax_left.tick_params(axis='both', which='major', labelsize=12)
    ax_right.tick_params(axis='y', which='major', labelsize=12)
    # 设置x轴刻度位置和标签
    ax_left.set_xticks(range(len(pivot_df)))
    ax_left.set_xticklabels(pivot_df['node_edge_label'], rotation=45, ha='right', fontsize=14)

    # 添加图例
    lines_left, labels_left = ax_left.get_legend_handles_labels()
    lines_right, labels_right = ax_right.get_legend_handles_labels()
    ax_left.legend(lines_left + lines_right, labels_left + labels_right, loc='best', fontsize=14)

    # 添加标题和网格
    ax_left.grid(True, alpha=0.3)

    # 调整布局
    plt.tight_layout()
    plt.savefig(statistical_magnitude_png, dpi=300, bbox_inches="tight")
    plt.savefig(statistical_magnitude_svg, dpi=300, bbox_inches='tight')




if __name__ == "__main__":
    main_path = "/home/gaoyun/poly/data"
    out_species = "GB54"
    species = "Brettanomyces_bruxellensis"
    phase_out_name = "Phasing"
    phased_name = "Brettanomyces_bruxellensisERR4624298"
    thread_number = 16

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

    out_path = os.path.join(phased_dir,"Analysed")
    out_plot_path = os.path.join(phased_path,"Plots")
    all_paths = [out_path,out_plot_path]
    for path in all_paths:
        os.makedirs(path,exist_ok=True)

    # reduced_read_file = os.path.join(phased_dir,"Variants/Reads/",phased_name+".hetPositions.SNPxLongReads.vcf.chr.tsv")
    # reads_dict = load_reduced_reads(reduced_read_file)

    # main(reads_dict, out_path, out_name, thread_number,min_ovl, min_sim, min_overlap, min_ovl_len)
    statistical_magnitude_file = os.path.join(out_path, out_name+"_"+"statistical_magnitude.tsv")
    plot_statistical_magnitude(statistical_magnitude_file, out_path, out_name)
