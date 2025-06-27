import os
from plotnine import *
import pandas as pd
import patchworklib as pw
from functools import reduce
from plotnine import ggplot, aes, geom_point, theme, guides, element_blank, element_text

# line graph of SaCere

def loadAccuracyBenchmarkMerge(benchmark_file,outPath,out_name,bench_keys):
    outputSVG = os.path.join(outPath,out_name+".svg")
    outputPNG=os.path.join(outPath,out_name+".png")

    gg_list = []
    tbl=pd.read_csv(benchmark_file,sep="\t")
    data=tbl
    for bench_key in bench_keys:
        if bench_key == "nHaplotigs per haplotype":
            ylabs = bench_key
        elif bench_key == "time":
            ylabs = "wall time (s)"
        else:
            ylabs = bench_key + "  " +"(%)"

        g = (ggplot(data,aes(x="ploidy_het_cov",y=bench_key,group="tools",shape="tools",color="tools")) +
            #  labs(title="mutation="+str(mut)+"%   "+"coverage="+str(cov)+"X")+ theme(plot_title=element_text(hjust=0.5))+
            labs(y=ylabs) +
            theme_classic() +
            labs(x="(ploidy, heterozygosityRates (%), coverage)") +
            scale_color_manual(values=["#006633", "#006699","#CC0033",  "#FF9900"]) +
            theme(legend_position="none",
                  axis_title_x=element_text(size=14),
                  axis_title_y=element_text(size=14),
                  axis_text_x=element_text(size=9, angle=45, hjust=1),
                  axis_text_y=element_text(size=11),)+
            geom_point(size=5) + geom_line(size=1))

        figsize = (6, 3)
        gg_list.append(pw.load_ggplot(g,figsize=figsize))
    
    n_cols = 2
    total_plots = len(gg_list)
    remainder = total_plots % n_cols
    if remainder != 0:
        blanks_needed = n_cols - remainder
        for _ in range(blanks_needed):
            blank_brick = ggplot_blank()
            gg_list.append(pw.load_ggplot(blank_brick,figsize=figsize))

    # 按行切割并拼接
    rows = []
    for i in range(0, len(gg_list), n_cols):
        row_plots = gg_list[i:i+n_cols]
        row = reduce(lambda a, b: a | b, row_plots)  # 横向拼接
        rows.append(row)

    # 纵向拼接所有行
    combined = reduce(lambda a, b: a / b, rows)

    # 保存或显示结果
    combined.savefig(outputPNG)
    combined.savefig(outputSVG)  

    # 生成图例专用图（作为标题）
    outputLegendSVG = os.path.join(outPath,out_name+"_legend.svg")
    outputLegendPNG=os.path.join(outPath,out_name+"_legend.png")
    legend_plot = generate_legend(data, "accuracy")
    ggsave(legend_plot,filename=outputLegendSVG)
    ggsave(legend_plot,filename=outputLegendPNG)

    return None

def ggplot_blank():
    return ggplot(pd.DataFrame()) + theme_void()

def generate_legend(data, bench_key):
    # 生成图例专用图（作为标题）
    return (
        ggplot(data,aes(x="ploidy_het",y=bench_key,group="tools",shape="tools",color="tools"))
        + geom_point(size=5)
        + geom_line(size=1)
        + scale_color_manual(values=["#006633", "#006699","#CC0033",  "#FF9900"])
        + theme_classic()
        # + guides(color=guide_legend(title="", title_position="top"))
        + theme(
            # 隐藏所有非图例元素
            panel_background=element_blank(),
            axis_line=element_blank(),
            axis_text=element_blank(),
            axis_ticks=element_blank(),
            plot_title=element_blank(),
            # 控制图例样式
            legend_position="top",
            # legend_title=element_text(size=14, weight="bold"),  # 标题样式
            legend_text=element_text(size=14),                   # 图例文字样式
            legend_background=element_blank(),                  # 透明背景
            legend_margin=0,                                     # 无边距
            plot_margin=0                                        # 全屏显示
        )
    )

if __name__ == "__main__":

    inputdir = "/home/gaoyun/poly/data/"
    benchmark_file = os.path.join(inputdir,"SaCere/accuracyMetrics/acc_time_vcfchrnx_round.tsv")
    bench_keys = ["accuracy","error rate","missing rate","nHaplotigs per haplotype","reads accuracy rate","time"]
    out_path = os.path.join(inputdir,"SaCere/accuracyMetrics/table/")
    os.makedirs(out_path,exist_ok=True)

    out_name = "all_SaCere_accuracy"
    loadAccuracyBenchmarkMerge(benchmark_file,out_path,out_name,bench_keys)

