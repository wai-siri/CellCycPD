# #!/usr/bin/env python
# # -*- coding: utf-8 -*-
# """
# 基因结构和增强子可视化工具

# 该工具可视化基因结构，包括外显子、内含子、启动子和增强子。
# 增强子数据来自两个来源：近端增强子(pELS)和远端增强子(dELS)。
# """

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.path import Path
import os
from typing import Dict, List, Tuple
import logging

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 重置matplotlib字体设置为默认值
plt.rcParams.update(plt.rcParamsDefault)


def read_gene_data(exon_file: str, intron_file: str, promoter_file: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """读取基因数据，包括外显子、内含子和启动子"""
    logger.info(f"Reading gene data: {exon_file}")
    exon_data = pd.read_csv(exon_file, sep='\t')
    
    logger.info(f"Reading intron data: {intron_file}")
    intron_data = pd.read_csv(intron_file, sep='\t')
    
    logger.info(f"Reading promoter data: {promoter_file}")
    promoter_data = pd.read_csv(promoter_file, sep='\t')
    
    # 移除长度为负的内含子
    if not intron_data.empty:
        intron_data = intron_data[
            intron_data['intron_chrom_end'] > intron_data['intron_chrom_start']
        ]
    
    return exon_data, intron_data, promoter_data


def read_enhancer_data(proximal_enhancer_file: str, distal_enhancer_file: str) -> Tuple[pd.DataFrame, pd.DataFrame]:

    logger.info(f"Reading proximal enhancer data: {proximal_enhancer_file}")
    # BED格式：染色体、起始位置、结束位置、名称、得分、类型
    proximal_enhancers = pd.read_csv(proximal_enhancer_file, sep='\t', header=None,
                                  names=['chromosome', 'start', 'end', 'id1', 'id2', 'type'])
    
    logger.info(f"Reading distal enhancer data: {distal_enhancer_file}")
    distal_enhancers = pd.read_csv(distal_enhancer_file, sep='\t', header=None,
                                names=['chromosome', 'start', 'end', 'id1', 'id2', 'type'])
    
    # 转换染色体格式以匹配基因数据
    proximal_enhancers['chromosome'] = proximal_enhancers['chromosome'].str.replace('chr', '')
    distal_enhancers['chromosome'] = distal_enhancers['chromosome'].str.replace('chr', '')
    
    return proximal_enhancers, distal_enhancers


def format_bp(bp: int) -> str:
    """格式化碱基对长度显示"""
    if bp >= 1000000:
        return f'{bp/1000000:.1f}Mb'
    elif bp >= 1000:
        return f'{bp/1000:.1f}kb'
    return f'{bp}bp'


def create_exon_shape(x: float, y: float, width: float, height: float) -> Path:
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


def find_gene_enhancers(gene_chrom: str, gene_start: int, gene_end: int, 
                       proximal_enhancers: pd.DataFrame, distal_enhancers: pd.DataFrame, 
                       window_size: int = 500) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """查找基因周围的增强子
    
    参数:
        gene_chrom: 基因所在的染色体
        gene_start: 基因起始位置
        gene_end: 基因结束位置
        proximal_enhancers: 近端增强子数据
        distal_enhancers: 远端增强子数据
        window_size: 基因上下游的窗口大小，默认500bp
        
    返回:
        (gene_proximal_enhancers, gene_distal_enhancers): 与基因相关的近端和远端增强子
    """
    # 扩展基因区域
    extended_start = max(0, gene_start - window_size)
    extended_end = gene_end + window_size
    
    # 筛选近端增强子
    gene_proximal_enhancers = proximal_enhancers[
        (proximal_enhancers['chromosome'] == gene_chrom) &
        (proximal_enhancers['start'] >= extended_start) &
        (proximal_enhancers['end'] <= extended_end)
    ]
    
    # 筛选远端增强子
    gene_distal_enhancers = distal_enhancers[
        (distal_enhancers['chromosome'] == gene_chrom) &
        (distal_enhancers['start'] >= extended_start) &
        (distal_enhancers['end'] <= extended_end)
    ]
    
    return gene_proximal_enhancers, gene_distal_enhancers


def plot_gene_structure_with_enhancers(exon_data: pd.DataFrame, intron_data: pd.DataFrame, 
                                     promoter_data: pd.DataFrame, proximal_enhancers: pd.DataFrame, 
                                     distal_enhancers: pd.DataFrame, output_dir: str = './gene_enhancer_structures'):
    """绘制带有增强子的基因结构图
    
    参数:
    exon_data: 包含位置和基因信息的外显子数据
    intron_data: 从相邻外显子派生的内含子数据
    promoter_data: 包含TSS位置的启动子数据
    proximal_enhancers: 近端增强子数据
    distal_enhancers: 远端增强子数据
    output_dir: 输出目录
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 设置图形样式
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # 获取所有基因名称
    genes = exon_data['external_gene_name'].unique()
    logger.info(f"Found {len(genes)} genes")
    
    for gene_name in genes:
        try:
            # 获取此基因的数据
            gene_exons = exon_data[exon_data['external_gene_name'] == gene_name].sort_values('exon_chrom_start')
            gene_introns = intron_data[intron_data['external_gene_name'] == gene_name].sort_values('intron_chrom_start')
            gene_promoters = promoter_data[promoter_data['external_gene_name'] == gene_name]
            
            if len(gene_exons) == 0:
                logger.warning(f"Skipping gene {gene_name}: no exon data")
                continue
                
            # 创建图形
            fig, ax = plt.subplots(figsize=(15, 4))
            fig.patch.set_facecolor('white')
            
            # 获取基因信息
            strand = gene_exons.iloc[0]['strand']
            chromosome = gene_exons.iloc[0]['chromosome_name']
            gene_start = min(gene_promoters['promoter_start'].min() if not gene_promoters.empty else float('inf'), 
                            gene_exons['exon_chrom_start'].min())
            gene_end = max(gene_promoters['promoter_end'].max() if not gene_promoters.empty else float('-inf'), 
                          gene_exons['exon_chrom_end'].max())
            gene_length = gene_end - gene_start
            
            # 查找基因周围的增强子
            gene_proximal_enhancers, gene_distal_enhancers = find_gene_enhancers(
                chromosome, gene_start, gene_end, proximal_enhancers, distal_enhancers, window_size=500)
            
            logger.info(f"Processing gene {gene_name} (Chromosome {chromosome}, {gene_start}-{gene_end}), " +
                       f"found {len(gene_proximal_enhancers)} proximal enhancers and {len(gene_distal_enhancers)} distal enhancers")
            
            # 计算压缩的相对位置
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
                
                return (compressed_pos - total_length) / (gene_length - total_length) if (gene_length - total_length) > 0 else 0
            
            # 绘制基因骨架（细线）
            ax.plot([0, 1], [0.5, 0.5], '#D3D3D3', linewidth=1, zorder=1)
            
            # 绘制启动子区域
            for _, promoter in gene_promoters.iterrows():
                # 计算启动子范围
                if strand == '+':
                    tss_pos = promoter['promoter_start']  # 正链上TSS在启动子起始位置
                    promoter_start = tss_pos - 200  # 上游200bp
                    promoter_end = tss_pos + 40    # 下游40bp
                else:
                    tss_pos = promoter['promoter_end']    # 负链上TSS在启动子结束位置
                    promoter_start = tss_pos - 40   # 上游40bp（相对于基因方向）
                    promoter_end = tss_pos + 200    # 下游200bp（相对于基因方向）
                
                # 转换为相对位置
                rel_tss_pos = compress_position(tss_pos)
                rel_start = compress_position(promoter_start)
                rel_end = compress_position(promoter_end)
                
                # Draw promoter region
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
                
                # 添加长度标注（仅对较长的外显子）
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
            
            # 绘制增强子
            # 近端增强子 - 绘制在上方
            for _, enhancer in gene_proximal_enhancers.iterrows():
                start_pos = compress_position(enhancer['start'])
                end_pos = compress_position(enhancer['end'])
                width = max(0.005, end_pos - start_pos)  # 确保可见性
                
                # 绘制近端增强子（上方）
                rect = patches.Rectangle(
                    (start_pos, 0.7), width, 0.15,
                    facecolor='#FFA500', edgecolor='#FF8C00',
                    linewidth=0.5, alpha=0.7, zorder=3
                )
                ax.add_patch(rect)
            
            # 远端增强子 - 绘制在下方
            for _, enhancer in gene_distal_enhancers.iterrows():
                start_pos = compress_position(enhancer['start'])
                end_pos = compress_position(enhancer['end'])
                width = max(0.005, end_pos - start_pos)  # 确保可见性
                
                # 绘制远端增强子（下方）
                rect = patches.Rectangle(
                    (start_pos, 0.15), width, 0.15,
                    facecolor='#32CD32', edgecolor='#228B22',
                    linewidth=0.5, alpha=0.7, zorder=3
                )
                ax.add_patch(rect)
            
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
            
            # 计算压缩的刻度位置
            compressed_ticks = [compress_position(pos) for pos in tick_positions]
            tick_labels = [f'{(x-gene_start)/divider:.1f}' for x in tick_positions]
            
            # 移除重叠的刻度
            final_ticks = []
            final_labels = []
            last_pos = float('-inf')
            min_distance = 0.1  # 最小距离（相对坐标）
            
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
            subtitle = f'Chromosome {chromosome}, Strand {strand}'
            ax.text(0.5, 1.15, title, ha='center', va='bottom', fontsize=12, fontweight='bold')
            ax.text(0.5, 1.05, subtitle, ha='center', va='bottom', fontsize=9, color='gray')
            
            # 添加图例
            legend_elements = [
                patches.Patch(facecolor='#4472C4', alpha=0.8, edgecolor='#2F528F', label='Exon'),
                plt.Line2D([0], [0], color='#D3D3D3', linewidth=1, label='Intron'),
                patches.Patch(facecolor='#FFB6C1', alpha=0.3, label='Promoter Region'),
                plt.Line2D([0], [0], marker='>' if strand == '+' else '<', color='#FF4500', 
                          label='TSS', markersize=8, linestyle='none'),
                patches.Patch(facecolor='#FFA500', alpha=0.7, label='Proximal Enhancer'),
                patches.Patch(facecolor='#32CD32', alpha=0.7, label='Distal Enhancer')
            ]
            
            ax.legend(
                handles=legend_elements,
                loc='upper center',
                frameon=True,
                facecolor='white',
                edgecolor='none',
                fontsize=8,
                bbox_to_anchor=(0.5, -0.15),
                ncol=3
            )
            
            # 保存图形
            plt.tight_layout()
            plt.subplots_adjust(bottom=0.2)  # 在底部添加额外空间用于图例
            plt.savefig(os.path.join(output_dir, f"{gene_name}.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            logger.error(f"Error processing gene {gene_name}: {str(e)}")
            plt.close()


def main():
    """主函数"""
    # 输入文件
    exon_file = "D:\\R_Project\\exon_data.tsv"
    intron_file = "D:\\R_Project\\intron_data.tsv"
    promoter_file = "D:\\R_Project\\promoter_data.tsv"
    proximal_enhancer_file = "D:\\R_Project\\cellcycle\\code\\tool\\gene_stru_code\\enhancer\\GRCh38-cCREs.pELS.bed"
    distal_enhancer_file = "D:\\R_Project\\cellcycle\\code\\tool\\gene_stru_code\\enhancer\\GRCh38-cCREs.dELS.bed"
    
    # 输出目录
    output_dir = "D:\\R_Project\\cellcycle\\code\\tool\\gene_enhancer_structures"
    
    # 读取数据
    logger.info("Starting data reading...")
    exon_data, intron_data, promoter_data = read_gene_data(exon_file, intron_file, promoter_file)
    proximal_enhancers, distal_enhancers = read_enhancer_data(proximal_enhancer_file, distal_enhancer_file)
    
    # 绘制基因结构图
    logger.info("Starting gene structure visualization...")
    plot_gene_structure_with_enhancers(
        exon_data, intron_data, promoter_data, 
        proximal_enhancers, distal_enhancers, output_dir
    )
    logger.info(f"Visualizations have been saved to {output_dir}")


if __name__ == "__main__":
    main()
