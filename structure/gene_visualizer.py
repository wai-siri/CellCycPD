import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib import font_manager
from matplotlib.path import Path
import os

def read_gene_data(exon_file, intron_file, promoter_file):
    """读取基因数据，包括外显子、内含子和启动子信息"""
    exon_data = pd.read_csv(exon_file, sep='\t')
    intron_data = pd.read_csv(intron_file, sep='\t')
    promoter_data = pd.read_csv(promoter_file, sep='\t')
    
    # 移除内含子长度为负的记录
    if not intron_data.empty:
        intron_data = intron_data[
            intron_data['intron_chrom_end'] > intron_data['intron_chrom_start']
        ]
    
    return exon_data, intron_data, promoter_data

def format_bp(bp):
    """格式化碱基对长度显示"""
    if bp >= 1000000:
        return f'{bp/1000000:.1f}Mb'
    elif bp >= 1000:
        return f'{bp/1000:.1f}kb'
    return f'{bp}bp'

def create_exon_shape(x, y, width, height):
    """创建外显子形状（矩形）"""
    vertices = [
        (x, y),
        (x, y + height),
        (x + width, y + height),
        (x + width, y),
        (x, y),
    ]
    codes = [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
    return Path(vertices, codes)

def plot_gene_structure(exon_data, intron_data, promoter_data, output_dir='./gene_structures'):
    """绘制基因结构图
    
    参数说明：
    exon_data: 外显子数据，包含位置和基因信息
    intron_data: 内含子数据，由相邻外显子推算
    promoter_data: 启动子数据，包含TSS位置
    output_dir: 输出目录
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 设置图形样式
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # 获取所有基因名称
    genes = exon_data['external_gene_name'].unique()
    
    for gene_name in genes:
        # 获取该基因的数据
        gene_exons = exon_data[exon_data['external_gene_name'] == gene_name].sort_values('exon_chrom_start')
        gene_introns = intron_data[intron_data['external_gene_name'] == gene_name].sort_values('intron_chrom_start')
        gene_promoters = promoter_data[promoter_data['external_gene_name'] == gene_name]
        
        if len(gene_exons) == 0:
            continue
            
        # 创建图形
        fig, ax = plt.subplots(figsize=(15, 3))
        fig.patch.set_facecolor('white')
        
        # 获取基因信息
        strand = gene_exons.iloc[0]['strand']
        gene_start = min(gene_promoters['promoter_start'].min(), gene_exons['exon_chrom_start'].min())
        gene_end = max(gene_promoters['promoter_end'].max(), gene_exons['exon_chrom_end'].max())
        gene_length = gene_end - gene_start
        
        # 计算压缩后的相对位置
        def compress_position(pos):
            """压缩长内含子，返回相对位置"""
            total_length = 0
            compressed_pos = pos - gene_start
            
            # 检查每个内含子
            for _, intron in gene_introns.iterrows():
                if pos > intron['intron_chrom_end']:
                    if intron['intron_chrom_end'] - intron['intron_chrom_start'] > 5000:
                        total_length += intron['intron_chrom_end'] - intron['intron_chrom_start'] - 300
                elif pos > intron['intron_chrom_start']:
                    if intron['intron_chrom_end'] - intron['intron_chrom_start'] > 5000:
                        ratio = (pos - intron['intron_chrom_start']) / (intron['intron_chrom_end'] - intron['intron_chrom_start'])
                        total_length += (pos - intron['intron_chrom_start']) - (ratio * 300)
            
            return (compressed_pos - total_length) / (gene_length - total_length)
        
        # 绘制基因主干（细线）
        ax.plot([0, 1], [0.5, 0.5], '#D3D3D3', linewidth=1, zorder=1)
        
        # 绘制启动子区域
        for _, promoter in gene_promoters.iterrows():
            # 计算启动子范围
            if strand == '+':
                tss_pos = promoter['promoter_start']  # 正链TSS在启动子起始位置
                promoter_start = tss_pos - 200  # 上游200bp
                promoter_end = tss_pos + 40    # 下游40bp
            else:
                tss_pos = promoter['promoter_end']    # 负链TSS在启动子末端位置
                promoter_start = tss_pos - 40   # 上游40bp（相对于基因方向）
                promoter_end = tss_pos + 200    # 下游200bp（相对于基因方向）
            
            # 转换为相对位置
            rel_tss_pos = compress_position(tss_pos)
            rel_start = compress_position(promoter_start)
            rel_end = compress_position(promoter_end)
            
            # 绘制启动子区域
            ax.axvspan(rel_start, rel_end, ymin=0.4, ymax=0.6, 
                      color='#FFB6C1', alpha=0.3, zorder=2)
            
            # 绘制TSS位点
            marker = '>' if strand == '+' else '<'
            ax.plot(rel_tss_pos, 0.5, marker=marker, color='#FF4500', 
                   markersize=8, zorder=4)
        
        # 绘制外显子（矩形）
        for i, exon in gene_exons.iterrows():
            start_pos = compress_position(exon['exon_chrom_start'])
            end_pos = compress_position(exon['exon_chrom_end'])
            width = end_pos - start_pos
            
            # 创建外显子形状
            exon_path = create_exon_shape(start_pos, 0.35, width, 0.3)
            exon_patch = patches.PathPatch(
                exon_path, facecolor='#4472C4', alpha=0.8,
                edgecolor='#2F528F', linewidth=0.5, zorder=3
            )
            ax.add_patch(exon_patch)
            
            # 添加长度标注（只标注较长的外显子）
            exon_length = exon['exon_chrom_end'] - exon['exon_chrom_start']
            if exon_length >= 1000:
                center = (start_pos + end_pos) / 2
                y_offset = 0.2 if i % 2 == 0 else -0.2
                y_text = 0.5 + y_offset
                
                ax.annotate(
                    format_bp(exon_length),
                    xy=(center, 0.5),
                    xytext=(center, y_text),
                    ha='center', va='center',
                    fontsize=6, color='#666666',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1),
                    arrowprops=dict(
                        arrowstyle='->', color='#666666',
                        connectionstyle=f'arc3,rad={0.2 if y_offset > 0 else -0.2}',
                        alpha=0.6
                    )
                )
        
        # 生成均匀分布的刻度位置
        num_ticks = 8  # 期望的刻度数量
        total_bp = gene_end - gene_start
        
        # 根据基因长度动态确定刻度间隔
        if total_bp > 1000000:
            base_step = 100000  # 100kb步长
            unit = 'Mb'
            divider = 1000000
        elif total_bp > 100000:
            base_step = 10000   # 10kb步长
            unit = 'kb'
            divider = 1000
        else:
            base_step = 1000    # 1kb步长
            unit = 'kb'
            divider = 1000
        
        step = max(base_step, total_bp // num_ticks)
        step = round(step / base_step) * base_step
        
        # 计算起始位置（四舍五入到步长的倍数）
        start_pos = (gene_start // step) * step
        
        # 生成刻度位置
        tick_positions = np.arange(start_pos, gene_end + step, step)
        tick_positions = tick_positions[tick_positions >= gene_start]
        tick_positions = tick_positions[tick_positions <= gene_end]
        
        # 计算压缩后的刻度位置
        compressed_ticks = [compress_position(pos) for pos in tick_positions]
        tick_labels = [f'{(x-gene_start)/divider:.1f}' for x in tick_positions]
        
        # 移除重叠的刻度
        final_ticks = []
        final_labels = []
        last_pos = float('-inf')
        min_distance = 0.1  # 最小间距（相对坐标）
        
        for pos, label in zip(compressed_ticks, tick_labels):
            if pos - last_pos >= min_distance:
                final_ticks.append(pos)
                final_labels.append(label)
                last_pos = pos
        
        # 添加刻度和标签
        ax.set_xticks(final_ticks)
        ax.set_xticklabels(final_labels, fontsize=8)
        ax.set_xlabel(f'Position ({unit})', fontsize=9)
        
        # 设置图形属性
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        
        # 设置标题
        title = f'{gene_name}'
        subtitle = f'Chromosome {gene_exons.iloc[0]["chromosome_name"]}, Strand {strand}'
        ax.text(0.5, 1.15, title, ha='center', va='bottom', fontsize=12, fontweight='bold')
        ax.text(0.5, 1.05, subtitle, ha='center', va='bottom', fontsize=9, color='gray')
        
        # 添加图例
        legend_elements = [
            patches.Patch(facecolor='#4472C4', alpha=0.8, edgecolor='#2F528F',
                        label='Exon'),
            plt.Line2D([0], [0], color='#D3D3D3', linewidth=1, label='Intron'),
            patches.Patch(facecolor='#FFB6C1', alpha=0.3, label='Promoter Region'),
            plt.Line2D([0], [0], marker='>' if strand == '+' else '<',
                      color='#FF4500', label='TSS', markersize=8, linestyle='none')
        ]
        
        ax.legend(
            handles=legend_elements,
            loc='upper right',
            frameon=True,
            facecolor='white',
            edgecolor='none',
            fontsize=8,
            bbox_to_anchor=(1, 1.2)
        )
        
        # 保存图片
        plt.savefig(os.path.join(output_dir, f'{gene_name}.png'),
                   bbox_inches='tight', dpi=300)
        plt.close()

def main():
    """主函数"""
    print("\nReading gene data...")
    # 输入文件
    exon_file = "D:\\R_Project\\exon_data.tsv"
    intron_file = "D:\\R_Project\\intron_data.tsv"
    promoter_file = "D:\\R_Project\\promoter_data.tsv"
    
    # 输出目录
    output_dir = "code/tool/gene_structures"
    
    # 读取数据
    exon_data, intron_data, promoter_data = read_gene_data(exon_file, intron_file, promoter_file)
    
    # 绘制基因结构图
    print("Generating gene structure visualizations...")
    plot_gene_structure(exon_data, intron_data, promoter_data, output_dir)
    print("Visualizations have been saved to ./gene_structures")

if __name__ == "__main__":
    main()
