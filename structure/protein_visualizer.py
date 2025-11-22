"""蛋白质特征可视化模块"""

import json
import logging
import numpy as np
from pathlib import Path
from typing import Dict, List
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from collections import defaultdict

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ProteinVisualizer:
    """蛋白质特征可视化器"""
    
    def __init__(self):
        self.feature_data_dir = Path("data/protein_features")
        self.output_dir = Path("data/protein_visualizations")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 定义颜色方案
        self.colors = {
            'Domain': '#FFB6C1',      # 结构域（粉红色）
            'Region': '#90EE90',      # 区域（浅绿色）
            'Helix': '#87CEEB',       # α螺旋（天蓝色）
            'Beta strand': '#FFA07A',  # β折叠（浅鲑色）
            'Turn': '#DDA0DD'         # 转角（梅红色）
        }
        
        # 定义修饰类型的颜色
        self.modification_colors = {
            'Acetylation': '#FF0000',    # 红色
            'Glycosylation': '#00FF00',  # 绿色
            'Methylation': '#0000FF',    # 蓝色
            'Phosphorylation': '#FFA500', # 橙色
            'SUMOylation': '#800080',    # 紫色
            'Ubiquitylation': '#008080', # 青色
            'Others': '#808080'          # 灰色
        }

        # 定义二级结构的显示位置
        self.structure_positions = {
            'Helix': 0.6,        # α螺旋显示在上层
            'Beta strand': 0.4,  # β折叠显示在中层
            'Turn': 0.2         # 转角显示在下层
        }

        # 定义二级结构的图例标签
        self.structure_labels = {
            'Helix': 'α-Helix',
            'Beta strand': 'β-Strand',
            'Turn': 'Turn'
        }

    def _get_modification_type(self, description: str) -> str:
        """根据描述确定修饰类型"""
        description = description.lower()
        if 'acetyl' in description:
            return 'Acetylation'
        elif 'glyco' in description:
            return 'Glycosylation'
        elif 'methyl' in description:
            return 'Methylation'
        elif 'phospho' in description or 'phosphoserine' in description or 'phosphothreonine' in description or 'phosphotyrosine' in description:
            return 'Phosphorylation'
        elif 'sumo' in description:
            return 'SUMOylation'
        elif 'ubiquit' in description:
            return 'Ubiquitylation'
        else:
            return 'Others'

    def _create_legend_elements(self, data: Dict) -> List:
        """创建图例元素"""
        legend_elements = []
        
        # 添加结构域图例
        if data.get('domains'):
            legend_elements.append(Patch(facecolor=self.colors['Domain'],
                                      edgecolor='black',
                                      label='Domain'))
        
        # 添加区域图例
        if data.get('regions'):
            legend_elements.append(Patch(facecolor=self.colors['Region'],
                                      edgecolor='black',
                                      label='Region'))
        
        # 添加二级结构图例
        if data.get('secondary_structure'):
            # 获取数据中存在的所有二级结构类型
            struct_types = {struct['type'] for struct in data['secondary_structure']}
            # 按固定顺序添加图例
            for struct_type in ['Helix', 'Beta strand', 'Turn']:
                if struct_type in struct_types:
                    legend_elements.append(Patch(facecolor=self.colors[struct_type],
                                              edgecolor='black',
                                              label=self.structure_labels[struct_type]))
        
        # 添加修饰类型图例
        mod_types = {self._get_modification_type(mod['description'])
                    for mod in data.get('modifications', [])}
        for mod_type in sorted(mod_types):
            legend_elements.append(Line2D([0], [0], marker='o', color='w',
                                       markerfacecolor=self.modification_colors[mod_type],
                                       label=mod_type, markersize=8))
        
        return legend_elements

    def visualize_protein(self, feature_file: Path) -> None:
        """可视化蛋白质特征"""
        # 读取特征数据
        with open(feature_file, 'r') as f:
            data = json.load(f)
            
        # 创建图形
        fig, ax = plt.subplots(figsize=(15, 3))
        
        # 绘制主链
        sequence_length = data['sequence_length']
        ax.plot([0, sequence_length], [0, 0], color='black', linewidth=2)
        
        # 绘制结构域
        if data.get('domains'):
            for domain in data['domains']:
                start = domain['start']
                end = domain['end']
                rect = Rectangle((start, -0.3), end-start, 0.6,
                               facecolor=self.colors['Domain'],
                               edgecolor='black',
                               alpha=0.7)
                ax.add_patch(rect)
        
        # 绘制区域
        if data.get('regions'):
            for region in data['regions']:
                start = region['start']
                end = region['end']
                rect = Rectangle((start, -0.2), end-start, 0.4,
                               facecolor=self.colors['Region'],
                               edgecolor='black',
                               alpha=0.7)
                ax.add_patch(rect)
        
        # 绘制二级结构
        if data.get('secondary_structure'):
            height = 0.15  # 二级结构的高度
            for struct in data['secondary_structure']:
                start = struct['start']
                end = struct['end']
                struct_type = struct.get('type', '')
                if struct_type in self.structure_positions:
                    y_pos = self.structure_positions[struct_type]
                    rect = Rectangle((start, y_pos), end-start, height,
                                   facecolor=self.colors[struct_type],
                                   edgecolor='black',
                                   alpha=0.8)
                    ax.add_patch(rect)
        
        # 绘制修饰位点
        for mod in data.get('modifications', []):
            mod_type = self._get_modification_type(mod['description'])
            pos = mod['start']
            color = self.modification_colors[mod_type]
            ax.plot(pos, 0, 'o', color=color, markersize=6)
        
        # 设置坐标轴
        ax.set_xlim(-50, sequence_length + 50)
        ax.set_ylim(-0.5, 1.0)
        ax.set_yticks([])
        ax.set_xlabel('Amino Acid Position')
        ax.set_title(f"{data['gene_name']} ({sequence_length} aa)")
        
        # 添加图例
        legend_elements = self._create_legend_elements(data)
        ax.legend(handles=legend_elements, loc='upper center',
                 bbox_to_anchor=(0.5, -0.15), ncol=5)
        
        # 保存图像
        output_file = self.output_dir / f"{data['gene_name']}.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"保存可视化结果到: {output_file}")

    def process_all_proteins(self) -> None:
        """处理所有蛋白质的特征文件"""
        feature_files = list(self.feature_data_dir.glob("*_features.json"))
        for feature_file in feature_files:
            logging.info(f"Processing {feature_file.stem}")
            self.visualize_protein(feature_file)

if __name__ == "__main__":
    # 示例使用
    visualizer = ProteinVisualizer()
    visualizer.process_all_proteins()
