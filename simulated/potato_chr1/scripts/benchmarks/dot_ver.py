import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba

# To validate coverage and heterozygosity levels on simulated potato chr1 datasets
# dot graphs for coverage and heterozygosity levels

def calculate_quantiles(x):
    return pd.Series({
        'q1': np.quantile(x, 0.25),
        'median': np.median(x),
        'q3': np.quantile(x, 0.75)
    })

def dot_fig(benchmark_file,outPath,out_name,bench_keys,filter_labels,filter_name,x_labels,x_groups,colors_list,label_name, group_name):

    outputSVG = os.path.join(outPath,out_name+".svg")
    outputPNG=os.path.join(outPath,out_name+".png")

    tbl=pd.read_csv(benchmark_file,sep="\t")
    if filter_labels == None:
        df = tbl
    else:
        df=tbl[(tbl[filter_name] == filter_labels)]

    x_lable1 = x_labels[0]

    for x_label_idx, x_label in enumerate(x_labels):
        for group in x_groups:
            if filter_labels == None:
                condition = (tbl[label_name] == x_label) & (tbl[group_name] == group)
            else:
                condition = (tbl[label_name] == x_label) & (tbl[group_name] == group) & (tbl[filter_name] == filter_labels)
            value = x_label_idx + (x_groups.index(group)-1)*0.2
            tbl.loc[condition, "x_pos"] = value

    # 创建大画布 (figsize根据需求调整)
    fig = plt.figure(figsize=(18, 10), facecolor='white')  # 画布尺寸和背景色
    fig.subplots_adjust(hspace=0.2, wspace=0.3)  # 控制子图间距

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

        stats = df.groupby([label_name, group_name])[bench_key].apply(calculate_quantiles).unstack().reset_index()

        colors = dict()
        for i, group in enumerate(x_groups):
            colors[group] = colors_list[i]

        box_height = 0.2 # 箱体宽度

        # 绘制散点图（保持原逻辑）
        for (x_label, group), group_df in df.groupby([label_name, group_name]):
            ax.scatter(
                y=group_df['x_pos'],  # 原x_pos改为y轴
                x=group_df[bench_key], 
                color=colors[group],
                s=40,
                edgecolor='white',
                alpha=0.9,
                label=group if x_label == x_lable1 else "",
                zorder=3  # 确保散点在上层
            )

        # 绘制统计箱型图（修改以下部分）
        for (x_label, group), stats_row in stats.groupby([label_name, group_name]):
            # 获取对应散点的坐标参数
            x_label_idx = x_labels.index(x_label)
            group_idx = x_groups.index(group)
            y_center = x_label_idx + (group_idx - 1)*0.2  # 复用散点偏移逻辑
            
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
        # 坐标轴优化
        ax.set_yticks(np.arange(len(x_labels))+0.1)  # 方法标签居中
        ax.set_yticklabels(x_labels, fontsize=12)
        ax.set_ylim(-0.5, len(x_labels)-0.5 + 0.2)  # 扩展边界包含所有箱型图
        ax.set_xlabel(ylabs, fontsize=16)
        ax.yaxis.grid(False)

        # 步骤2：获取原始刻度位置并计算新位置
        original_yticks = np.arange(len(x_labels)) + 0.1  # 原代码中的y刻度（已上调0.1）
        adjusted_ticks = original_yticks + 0.5  # 继续上调0.5

        # 步骤3：在调整后的位置绘制水平线
        for y in adjusted_ticks:
            ax.axhline(
                y, 
                color='white', 
                linestyle='--', 
                linewidth=0.8, 
                alpha=0.9, 
                zorder=0
            )
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

    # 假设每个子图通过 label 参数定义标签
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
        bbox_to_anchor=(0.5, 0.93)  # 调整位置
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

    bench_keys = ["accuracy","error rate","missing rate","nHaplotigs per haplotype","reads accuracy rate","time"]
    inputdir = "/home/gaoyun/poly/data/"
    benchmark_file = os.path.join(inputdir,"potato_chr1/accuracyMetrics/acc_time_lofreq_vcfchrnx_round.tsv")
    out_path = os.path.join(inputdir,"potato_chr1/accuracyMetrics/table/")
    os.makedirs(out_path,exist_ok=True)
    mutations = [0.05,0.1,0.3,1.0,3.0]
    covs = [10,15,20]
    ploidys = [3,4,5,6,7,8,9,10,11,12]
    methods = ["flopp","nPhase","nTChap","whatshap polyphase"]
    colors_list = ["#006633", "#006699", "#CC0033", "#FF9900"]

    out_name = "CovPloidy"
    mut = None
    dot_fig(benchmark_file,out_path,out_name,bench_keys,mut,"heterozygosityRates",covs,methods,colors_list,"coverage","tools")

    out_name = "HetPloidy"
    cov = None
    dot_fig(benchmark_file,out_path,out_name,bench_keys,cov,"coverage",mutations,methods,colors_list,"heterozygosityRates","tools")
