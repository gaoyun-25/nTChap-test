import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba

# To validate the post-processing step
# plot statistics results figure for potato chr1 datasets and their post-processing results

# 计算四分位数统计量（保持原逻辑）
def calculate_quantiles(x):
    return pd.Series({
        'q1': np.quantile(x, 0.25),
        'median': np.median(x),
        'q3': np.quantile(x, 0.75)
    })


def dot_fig(benchmark_file,outPath,out_name,bench_keys,x_groups,colors_list,group_name):

    outputSVG = os.path.join(outPath,out_name+".svg")
    outputPNG=os.path.join(outPath,out_name+".png")

    tbl=pd.read_csv(benchmark_file,sep="\t")

    df = tbl
    for group in x_groups:
        condition = (tbl[group_name] == group)
        value = (x_groups.index(group)-1)*0.2
        tbl.loc[condition, "x_pos"] = value

    # 创建大画布 (figsize根据需求调整)
    fig = plt.figure(figsize=(15, 6), facecolor='white')  # 画布尺寸和背景色
    fig.subplots_adjust(hspace=0.3, wspace=0.3)  # 控制子图间距

    # 创建 2行×3列 的子图网格
    axes = []
    for i in range(len(bench_keys)):  # 6个子图
        ax = fig.add_subplot(2, 3, i+1)  # 位置索引从1开始
        axes.append(ax)

        bench_key = bench_keys[i]
        if bench_key == "nHaplotigs per haplotype":
            ylabs = bench_key
        elif bench_key == "time":
            ylabs = "wall time (s)"
        else:
            ylabs = bench_key + "  " +"(%)"

        stats = df.groupby([group_name])[bench_key].apply(calculate_quantiles).unstack().reset_index()

        colors = dict()
        for i, group in enumerate(x_groups):
            colors[group] = colors_list[i]

        box_height = 0.15 # 箱体宽度

        # 绘制散点图（保持原逻辑）
        for group, group_df in df.groupby(group_name):
            ax.scatter(
                y=group_df['x_pos'],  # 原x_pos改为y轴
                x=group_df[bench_key], 
                color=colors[group],
                s=40,
                edgecolor='white',
                alpha=0.9,
                label=group,
                zorder=3  # 确保散点在上层
            )
        ax.yaxis.set_ticklabels([])

        # 绘制统计箱型图（修改以下部分）
        for group, stats_row in stats.groupby(group_name):
            # 获取对应散点的坐标参数
            # x_label_idx = x_labels.index(x_label)
            group_idx = x_groups.index(group)
            y_center = (group_idx - 1)*0.2  # 复用散点偏移逻辑
            
            # 颜色参数设置
            edge_color = colors[group]
            face_color = to_rgba(edge_color, alpha=0.4)  # 添加透明度
            
            # 绘制IQR箱体
            rect = plt.Rectangle(
                (stats_row['q1'].values[0], y_center - box_height/2),  # x起始=q1, y位置调整
                stats_row['q3'].values[0] - stats_row['q1'].values[0],  # 宽度=q3-q1
                box_height, 
                # edgecolor=edge_color,
                facecolor=face_color,  # 使用带透明度的填充色
                linewidth=2,
                zorder=2
            )
            ax.add_patch(rect)
            
            # 绘制中位线（保持原颜色不透明）
            ax.vlines(
                stats_row['median'].values[0],
                y_center - box_height/2,
                y_center + box_height/2,
                colors=edge_color,
                linewidths=2,
                zorder=2
            )

        ax.set_xlabel(ylabs, fontsize=12)

        ax.yaxis.grid(False)

        ax.xaxis.grid(
            True,
            color="white",    # 白色
            linestyle="-",     # 实线
            linewidth=0.8,     # 线宽
            alpha=0.9,         # 不透明度（0=透明，1=不透明）
            zorder=0           # 确保网格在数据点下方
        )

        # 去除四周边框线
        for spine in ax.spines.values():
            spine.set_visible(False)

        # 可选：去除刻度线（保留标签文字）
        ax.tick_params(axis='both', 
                    which='both',
                    length=0,        # 刻度线长度设为0
                    width=0)   

        # # 图例设置（保持原逻辑）
        handles, labels = ax.get_legend_handles_labels()
        legend = ax.legend(
            handles[:4], labels[:4],
            title="",
            bbox_to_anchor=(1.02, 1),
            loc='upper left',
            borderaxespad=0.
        )

        ax.set_facecolor('#eeeeee')

    handles_all, labels_all = [], []

    for ax in axes:
        # 禁止显示子图图例
        ax.legend_.remove() if ax.legend_ else None  # 清除已有图例

        # 收集句柄和标签
        handles, labels = ax.get_legend_handles_labels()
        handles_all.extend(handles)
        labels_all.extend(labels)

    # 去重处理（如果多个子图有相同标签）
    unique_labels = pd.unique(labels_all)
    unique_handles = [handles_all[labels_all.index(label)] for label in unique_labels]
    # 在大图顶部添加全局图例（可选）
    fig.legend(
        unique_handles, unique_labels,
        loc='upper center',
        ncol=4,
        fontsize=16,
        frameon=False,
        bbox_to_anchor=(0.5, 0.96)  # 调整位置
    )

    # 保存完整大图
    plt.savefig(
        outputPNG,
        dpi=300,
        bbox_inches="tight",
        facecolor='white'  # 与画布背景一致
    )
    plt.savefig(
        outputSVG,
        bbox_inches="tight",
        facecolor='white'  # 与画布背景一致
    )

if __name__ == "__main__":

    bench_keys = ["accuracy","error rate","missing rate","nHaplotigs per haplotype","reads accuracy rate"]
    inputdir = "/home/gaoyun/poly/code/phasing/process_reads/GB54/"
    benchmark_file = os.path.join(inputdir,"potato_chr1_acc_lofreq_vcfchrnx_round.tsv")
    out_path = "/home/gaoyun/poly/data/GB54/Phasing_vcf/Brettanomyces_bruxellensisERR4624298/Analysed"
    os.makedirs(out_path,exist_ok=True)
    mutations = [0.05,0.1,0.3,1.0,3.0]
    covs = [10,15,20]
    ploidys = [3,4,5,6,7,8,9,10,11,12]
    methods = ["nTChap","nTChap post-processing"]
    colors_list = ["#006699", "#FF6B35"]

    out_name = "All_potato_chr1"
    cov = None
    mutations = None
    dot_fig(benchmark_file,out_path,out_name,bench_keys,methods,colors_list,"tools")
