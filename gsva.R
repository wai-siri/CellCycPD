# 安装必需的R包（如果尚未安装）
#install.packages("BiocManager")
#BiocManager::install(c("GSVA", "GSEABase", "limma", "Biobase"))

# 加载包
library(GSVA)
library(GSEABase)
library(limma)
library(Biobase)

# 设置工作目录（可选，根据需要修改）
# setwd("d:/R_Project/cellcycle/code/tool/GSVA")

# 从txt文件加载自定义基因集
# 读取基因列表文件
gene_list <- readLines("D:/R_Project/cellcycle/code/tool/GSVA/1.txt")
# 移除空行
gene_list <- gene_list[gene_list != ""]
# 创建基因集列表
gene_set_list <- list(CellCycle_Genes = gene_list)

# 将基因集列表转换为GSVA可用的GeneSetCollection对象
gene_sets <- lapply(names(gene_set_list), function(set_name) {
  genes <- gene_set_list[[set_name]]
  new("GeneSet", setName = set_name, geneIds = genes)
})
gene_sets <- GeneSetCollection(gene_sets)

# 查看基因集的内容
print(gene_sets)

# 读取TPM格式的表达数据
expression_data <- read.table("D:/R_Project/cellcycle/code/tool/GSVA/mRNA_exp_TPM_only_TCGA/TCGA-UVM.gene_expression_TPM.tsv", header = TRUE, sep = "\t", row.names = 1)

# 查看数据
print(paste("表达数据维度：", nrow(expression_data), "基因 x", ncol(expression_data), "样本"))
print("表达数据前几行：")
head(expression_data[, 1:5])  # 只显示前5个样本的数据

# TPM数据需要进行log2转换（添加一个小值避免log2(0)）
log2_expression_data <- log2(expression_data + 0.01)

# 计算基因集评分
print("开始计算基因集评分...")
gsva_scores <- gsva(as.matrix(log2_expression_data), gene_sets, method = "gsva", kcdf = "Gaussian", abs.ranking = TRUE)

# 查看结果（得到每个样本的基因集评分）
print("GSVA评分结果前几个样本：")
print(gsva_scores[, 1:5])

# 提取CellCycle基因集的评分
cellcycle_score <- gsva_scores["CellCycle_Genes", ]

# 计算评分的中位数
median_score <- median(cellcycle_score)

# 根据中位数切分样本
high_score_group <- which(cellcycle_score > median_score)
low_score_group <- which(cellcycle_score <= median_score)

# 查看高低评分组的样本
cat("高评分组：", length(high_score_group), "个样本\n")
cat("低评分组：", length(low_score_group), "个样本\n")

# 创建分组信息
groups <- rep("Low", length(cellcycle_score))
groups[high_score_group] <- "High"

# 绘制箱线图
pdf("D:/R_Project/cellcycle/code/tool/GSVA/33cancerGSVA_Boxplot_result/UVM_Boxplot.pdf", width = 6, height = 5)
boxplot(cellcycle_score ~ factor(groups), col = c("lightblue", "lightgreen"), 
        names = c("Low Score", "High Score"), main = "Cell Cycle Pathway Score by Group")
dev.off()

# 保存评分结果
result_df <- data.frame(Sample = colnames(expression_data), 
                       CellCycle_Score = cellcycle_score,
                       Group = groups)
write.table(result_df, file = "D:/R_Project/cellcycle/code/tool/GSVA/33cancerGSVA_Boxplot_result/UVM_Results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# 保存原始GSVA评分矩阵
write.table(t(gsva_scores), file = "D:/R_Project/cellcycle/code/tool/GSVA/33cancerGSVA/UVM.txt", sep = "\t", quote = FALSE)
