import pandas as pd
import numpy as np

def analyze_gene_data(exon_file, intron_file):
    """分析基因数据中的重复和重叠问题"""
    print("读取数据文件...")
    exon_data = pd.read_csv(exon_file, sep='\t')
    intron_data = pd.read_csv(intron_file, sep='\t')
    
    print("\n基本统计信息:")
    print(f"外显子数量: {len(exon_data)}")
    print(f"内含子数量: {len(intron_data)}")
    print(f"基因数量: {len(exon_data['external_gene_name'].unique())}")
    
    # 检查重复行
    print("\n检查重复行:")
    exon_duplicates = exon_data[exon_data.duplicated()]
    if len(exon_duplicates) > 0:
        print(f"发现 {len(exon_duplicates)} 个重复的外显子记录")
        print(exon_duplicates)
    else:
        print("未发现完全重复的外显子记录")
    
    # 检查每个基因的外显子重叠情况
    print("\n检查外显子重叠:")
    overlaps = []
    for gene in exon_data['external_gene_name'].unique():
        gene_exons = exon_data[exon_data['external_gene_name'] == gene].sort_values('exon_chrom_start')
        
        if len(gene_exons) > 1:
            for i in range(len(gene_exons)-1):
                current = gene_exons.iloc[i]
                next_exon = gene_exons.iloc[i+1]
                
                if current['exon_chrom_end'] > next_exon['exon_chrom_start']:
                    overlaps.append({
                        'gene': gene,
                        'exon1_id': current['ensembl_exon_id'],
                        'exon2_id': next_exon['ensembl_exon_id'],
                        'overlap_size': current['exon_chrom_end'] - next_exon['exon_chrom_start']
                    })
    
    if overlaps:
        print(f"发现 {len(overlaps)} 个重叠的外显子对:")
        overlaps_df = pd.DataFrame(overlaps)
        print(overlaps_df)
    else:
        print("未发现重叠的外显子")
    
    # 检查内含子长度分布
    print("\n内含子长度分布:")
    intron_lengths = intron_data['intron_chrom_end'] - intron_data['intron_chrom_start']
    print(f"最短内含子: {intron_lengths.min():,} bp")
    print(f"最长内含子: {intron_lengths.max():,} bp")
    print(f"中位数长度: {intron_lengths.median():,} bp")
    print(f"平均长度: {intron_lengths.mean():,.0f} bp")
    print("\n长度分位数:")
    print(intron_lengths.describe())
    
    # 输出建议的压缩阈值
    q75 = intron_lengths.quantile(0.75)
    suggested_threshold = max(10000, q75)
    print(f"\n建议的内含子压缩阈值: {suggested_threshold:,.0f} bp")
    print(f"超过该阈值的内含子数量: {sum(intron_lengths > suggested_threshold)}")

if __name__ == "__main__":
    exon_file = "D:\\R_Project\\exon_data.tsv"
    intron_file = "D:\\R_Project\\intron_data.tsv"
    analyze_gene_data(exon_file, intron_file)
