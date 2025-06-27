import pandas as pd
import matplotlib.pyplot as plt


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
    ax_left.set_xlabel('Chromosomes and (number of nodes, number of edges)',fontsize=16)
    ax_left.set_ylabel('Statistical Magnitudes (0-1 scale)', color="#083046", fontsize=16)
    ax_right.set_ylabel('Avg. Node Strength (60-80 scale)', color='purple', fontsize=16)

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
    ax_left.legend(lines_left + lines_right, labels_left + labels_right, loc='best')

    # 添加标题和网格
    ax_left.grid(True, alpha=0.3)

    # 调整布局
    plt.tight_layout()
    plt.savefig(statistical_magnitude_png, dpi=300, bbox_inches="tight")
    plt.savefig(statistical_magnitude_svg, dpi=300, bbox_inches='tight')
