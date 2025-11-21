import pandas as pd

def compare_gene_lists(intersection_file: str, reactome_file: str) -> None:
    """
    比对intersection.csv中的基因是否都在reactome.csv中出现
    
    Args:
        intersection_file (str): intersection.csv的文件路径
        reactome_file (str): reactome.csv的文件路径
    """
    # 读取文件
    try:
        intersection_df = pd.read_csv(intersection_file)
        reactome_df = pd.read_csv(reactome_file)
        
        # 获取基因列表
        intersection_genes = set(intersection_df['gene_symbol'])
        reactome_genes = set(reactome_df.iloc[:, 0])  # 第一列为gene symbol
        
        # 计算统计信息
        total_intersection_genes = len(intersection_genes)
        genes_in_reactome = intersection_genes.intersection(reactome_genes)
        genes_not_in_reactome = intersection_genes - reactome_genes
        
        # 打印结果
        print(f"\n基因比对结果:")
        print(f"intersection.csv中共有 {total_intersection_genes} 个基因")
        print(f"其中 {len(genes_in_reactome)} 个基因在reactome.csv中出现")
        print(f"有 {len(genes_not_in_reactome)} 个基因未在reactome.csv中出现")
        
        # 如果有未出现的基因，将它们保存到文件中
        if genes_not_in_reactome:
            print("\n未在reactome.csv中出现的基因:")
            for gene in sorted(genes_not_in_reactome):
                # 获取该基因在intersection.csv中的数据库来源信息
                source_info = intersection_df[intersection_df['gene_symbol'] == gene]['databases'].iloc[0]
                print(f"基因: {gene}, 来源: {source_info}")
            
            # 保存未出现的基因到文件
            missing_genes_df = intersection_df[intersection_df['gene_symbol'].isin(genes_not_in_reactome)]
            output_file = r"D:\R_Project\genes_not_in_reactome.csv"
            missing_genes_df.to_csv(output_file, index=False)
            print(f"\n未出现的基因详细信息已保存到: {output_file}")
            
    except Exception as e:
        print(f"处理文件时出错: {str(e)}")

def main():
    # 文件路径
    intersection_file = r"D:\R_Project\intersection.csv"
    reactome_file = r"D:\R_Project\reatome.csv"
    
    # 执行比对
    compare_gene_lists(intersection_file, reactome_file)

if __name__ == "__main__":
    main()
