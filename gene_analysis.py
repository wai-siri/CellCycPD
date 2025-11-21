import pandas as pd
from typing import List, Dict, Set
import os
from itertools import combinations

class GeneAnalyzer:
    def __init__(self, file_paths: Dict[str, str]):
        """
        初始化基因分析器
        
        Args:
            file_paths (Dict[str, str]): 数据文件路径字典，键为数据库名称，值为文件路径
        """
        self.file_paths = file_paths
        self.gene_sets: Dict[str, Set[str]] = {}
        self.output_dir = os.path.dirname(list(file_paths.values())[0])  # 使用输入文件的目录作为输出目录
        self.load_data()

    def load_data(self):
        """加载所有数据集的基因符号"""
        for db_name, file_path in self.file_paths.items():
            try:
                df = pd.read_csv(file_path)
                # 获取第一列（gene_symbol）的值，转换为集合，排除表头
                genes = set(df.iloc[:, 0].dropna().str.strip().unique())
                self.gene_sets[db_name] = genes
                print(f"从 {db_name} 加载了 {len(genes)} 个基因")
            except Exception as e:
                print(f"加载 {db_name} 数据时出错: {str(e)}")

    def find_intersections(self) -> Set[str]:
        """
        找出所有数据集的交集基因
        
        Returns:
            Set[str]: 交集基因集合
        """
        if not self.gene_sets:
            return set()
        
        # 从所有集合中获取交集
        intersection = set.intersection(*self.gene_sets.values())
        return intersection

    def find_three_dataset_intersections(self) -> Dict[str, Set[str]]:
        """
        找出在任意三个数据集中共同出现的基因
        
        Returns:
            Dict[str, Set[str]]: 数据集组合到基因集合的映射
        """
        results = {}
        # 获取所有可能的三个数据集的组合
        for combo in combinations(self.gene_sets.keys(), 3):
            # 计算这三个数据集的交集
            intersection = set.intersection(*(self.gene_sets[db] for db in combo))
            # 排除在所有四个数据集中都出现的基因
            all_intersection = self.find_intersections()
            unique_to_three = intersection - all_intersection
            if unique_to_three:
                results[' & '.join(combo)] = unique_to_three
        return results

    def get_high_confidence_genes(self) -> pd.DataFrame:
        """
        获取在至少三个数据集中出现的所有基因
        
        Returns:
            pd.DataFrame: 包含基因符号和来源数据集的DataFrame
        """
        # 获取四数据集交集
        all_intersection = self.find_intersections()
        # 获取三数据集交集
        three_set_intersections = self.find_three_dataset_intersections()
        
        # 合并结果
        result_data = []
        
        # 添加四数据集交集的基因
        for gene in all_intersection:
            result_data.append({
                'gene_symbol': gene,
                'databases': 'All Databases',
                'intersection_count': 4
            })
        
        # 添加三数据集交集的基因
        for combo, genes in three_set_intersections.items():
            for gene in genes:
                result_data.append({
                    'gene_symbol': gene,
                    'databases': combo,
                    'intersection_count': 3
                })
        
        # 转换为DataFrame并按基因符号排序
        df = pd.DataFrame(result_data)
        return df.sort_values(['intersection_count', 'gene_symbol'], ascending=[False, True])

    def save_results(self):
        """保存分析结果到文件"""
        # 保存高置信度基因（在至少三个数据集中出现的基因）
        high_confidence_df = self.get_high_confidence_genes()
        intersection_file = os.path.join(self.output_dir, 'intersection.csv')
        high_confidence_df.to_csv(intersection_file, index=False)
        print(f"\n已将高置信度基因保存到: {intersection_file}")

def main():
    # 定义数据文件路径
    file_paths = {
        'Reactome': r"D:\R_Project\reatome.csv",
        'WikiPathway': r"D:\R_Project\wikipathway.csv",
        'Liang': r"D:\R_Project\liang.csv",
        'KEGG': r"D:\R_Project\KEGG.csv"
    }

    # 初始化分析器并执行分析
    analyzer = GeneAnalyzer(file_paths)
    analyzer.save_results()

    # 打印分析结果摘要
    intersection = analyzer.find_intersections()
    three_set_intersections = analyzer.find_three_dataset_intersections()

    print("\n分析结果摘要:")
    print(f"在所有数据集中共同出现的基因数量: {len(intersection)}")
    print("\n在三个数据集中共同出现的基因数量:")
    for combo, genes in three_set_intersections.items():
        print(f"数据集 {combo}: {len(genes)} 个基因")

if __name__ == "__main__":
    main()
