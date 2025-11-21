# 载入需要的包
library(treeio)
library(ggtree)
library(tidytree)
library(stringr)
library(ggplot2)

# 设置工作目录
setwd("D:/R_Project/cellcycle/code/tool/tree")

# 导入Newick树文件与多序列比对FASTA文件
# 使用相对路径，因为已经设置了工作目录
tree <- read.tree("nwk_data/ABL1.nwk")
msa <- read.fasta("fas_data/ABL1_CDS_all.fas")

# 从文件名中提取基因名称作为输出文件名的一部分
gene_name <- basename(tools::file_path_sans_ext("nwk_data/ABL1.nwk"))

# 创建输出目录
output_dir <- "res/png"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 处理空格问题，替换为空格的基因名
names(msa) <- str_replace(names(msa), " ", "_")

# 确保只使用同时在树和序列中都出现的物种
# 获取树中的物种名称
tree_taxa <- tree$tip.label
# 获取序列中的物种名称
msa_taxa <- names(msa)

# 找出同时在树和序列中都出现的物种
common_taxa <- intersect(tree_taxa, msa_taxa)

# 如果有不匹配的物种，输出警告信息
if (length(common_taxa) < length(tree_taxa) || length(common_taxa) < length(msa_taxa)) {
  cat("\n\n警告: 树和序列中的物种不完全匹配\n")
  cat("  树中物种数量: ", length(tree_taxa), "\n")
  cat("  序列中物种数量: ", length(msa_taxa), "\n")
  cat("  共有物种数量: ", length(common_taxa), "\n")
  
  # 显示不在序列中的树物种
  tree_only <- setdiff(tree_taxa, msa_taxa)
  if (length(tree_only) > 0) {
    cat("  只在树中出现的物种: ", paste(tree_only, collapse=", "), "\n")
  }
  
  # 显示不在树中的序列物种
  msa_only <- setdiff(msa_taxa, tree_taxa)
  if (length(msa_only) > 0) {
    cat("  只在序列中出现的物种: ", paste(msa_only, collapse=", "), "\n")
  }
  
  # 过滤树，只保留共有物种
  tree <- drop.tip(tree, setdiff(tree_taxa, common_taxa))
  
  # 过滤序列，只保留共有物种
  msa <- msa[names(msa) %in% common_taxa]
  
  cat("  已自动过滤不匹配的物种\n\n")
}

# 先绘制进化树
p1 <- ggtree(tree, branch.length = "none", layout = "rectangular") + 
  geom_tiplab(size = 3, color = "blue", fontface = "bold", align = TRUE) +  # 设置标签大小和颜色
  theme_tree() +  # 使用ggtree的主题
  theme(legend.position = "none")  # 去除图例

# 美化树的分支：给树枝加上颜色渐变效果
p1 <- p1 + 
  geom_tree(aes(color = branch.length), size = 1) +  # 使用geom_tree来设置颜色
  scale_color_gradient(low = "gray", high = "red")  # 为树枝加渐变色

# 然后使用msaplot将多序列比对与树形合并
# 增加offset参数值，以增加树和多序列比对之间的间距
p2 <- msaplot(p1, msa, offset = 8, width = 2)

# 可以进一步调整多序列比对的可视化效果
p2 <- p2 + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"),
        axis.title.x = element_text(size = 10, face = "bold", color = "black")) + 
  ggtitle("Phylogenetic Tree and MSA") +  # 添加标题
  theme(plot.title = element_text(size = 15, hjust = 0.5))  # 调整标题样式

# 打印最终结果
print(p2)

# 生成高清晰度PNG图片
# 设置输出文件路径
output_file <- file.path(output_dir, paste0(gene_name, "_tree_msa.png"))

# 使用ggsave保存高清晰度PNG图片
# dpi设置为300，这是高清晰度的标准
# 根据序列数量自动调整高度
# 计算最长的物种名称长度，用于动态调整图片宽度
max_label_length <- max(nchar(tree$tip.label))
# 根据最长物种名称长度调整宽度，确保有足够空间显示物种名称
width_factor <- max(1, max_label_length / 15)  # 每15个字符增加1英寸宽度

ggsave(
  filename = output_file,
  plot = p2,
  width = 14 + width_factor,  # 增加基础宽度，并根据物种名称长度动态调整
  height = max(8, length(tree$tip.label) * 0.2),  # 根据物种数量动态调整高度
  dpi = 300,  # 高清晰度
  bg = "white",  # 白色背景
  limitsize = FALSE  # 允许生成大尺寸图片
)

# 输出成功信息
cat("\n\n高清晰度PNG图片已生成: ", output_file, "\n")

# 添加参数化功能，可以通过命令行参数指定基因名称
run_with_gene <- function(gene_name) {
  # 构建文件路径
  tree_file <- paste0("nwk_data/", gene_name, ".nwk")
  msa_file <- paste0("fas_data/", gene_name, "_CDS_all.fas")
  
  # 检查文件是否存在
  if (!file.exists(tree_file) || !file.exists(msa_file)) {
    cat("\n错误: 无法找到文件 ", tree_file, " 或 ", msa_file, "\n")
    return(FALSE)
  }
  
  # 读取文件
  tree <- read.tree(tree_file)
  msa <- read.fasta(msa_file)
  
  # 处理空格问题
  names(msa) <- str_replace(names(msa), " ", "_")
  
  # 确保只使用同时在树和序列中都出现的物种
  # 获取树中的物种名称
  tree_taxa <- tree$tip.label
  # 获取序列中的物种名称
  msa_taxa <- names(msa)
  
  # 找出同时在树和序列中都出现的物种
  common_taxa <- intersect(tree_taxa, msa_taxa)
  
  # 如果有不匹配的物种，输出警告信息
  if (length(common_taxa) < length(tree_taxa) || length(common_taxa) < length(msa_taxa)) {
    cat("\n\n警告: ", gene_name, " 基因的树和序列中的物种不完全匹配\n")
    cat("  树中物种数量: ", length(tree_taxa), "\n")
    cat("  序列中物种数量: ", length(msa_taxa), "\n")
    cat("  共有物种数量: ", length(common_taxa), "\n")
    
    # 显示不在序列中的树物种
    tree_only <- setdiff(tree_taxa, msa_taxa)
    if (length(tree_only) > 0) {
      cat("  只在树中出现的物种: ", paste(tree_only, collapse=", "), "\n")
    }
    
    # 显示不在树中的序列物种
    msa_only <- setdiff(msa_taxa, tree_taxa)
    if (length(msa_only) > 0) {
      cat("  只在序列中出现的物种: ", paste(msa_only, collapse=", "), "\n")
    }
    
    # 过滤树，只保留共有物种
    tree <- drop.tip(tree, setdiff(tree_taxa, common_taxa))
    
    # 过滤序列，只保留共有物种
    msa <- msa[names(msa) %in% common_taxa]
    
    cat("  已自动过滤不匹配的物种\n\n")
  }
  
  # 绘制进化树
  p1 <- ggtree(tree, branch.length = "none", layout = "rectangular") + 
    geom_tiplab(size = 3, color = "blue", fontface = "bold", align = TRUE) +
    theme_tree() +
    theme(legend.position = "none")
  
  # 美化树的分支
  p1 <- p1 + 
    geom_tree(aes(color = branch.length), size = 1) +
    scale_color_gradient(low = "gray", high = "red")
  
  # 合并多序列比对与树形
  # 使用更大的offset值，确保物种名称不被遮挡
  p2 <- msaplot(p1, msa, offset = 8, width = 2)
  
  # 调整多序列比对的可视化效果
  p2 <- p2 + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"),
          axis.title.x = element_text(size = 10, face = "bold", color = "black")) + 
    ggtitle(paste0(gene_name, " Phylogenetic Tree and MSA")) +
    theme(plot.title = element_text(size = 15, hjust = 0.5))
  
  # 设置输出文件路径
  output_file <- file.path(output_dir, paste0(gene_name, "_tree_msa.png"))
  
  # 计算最长的物种名称长度，用于动态调整图片宽度
  max_label_length <- max(nchar(tree$tip.label))
  # 根据最长物种名称长度调整宽度，确保有足够空间显示物种名称
  width_factor <- max(1, max_label_length / 15)  # 每15个字符增加1英寸宽度
  
  # 保存高清晰度PNG图片
  ggsave(
    filename = output_file,
    plot = p2,
    width = 14 + width_factor,  # 增加基础宽度，并根据物种名称长度动态调整
    height = max(8, length(tree$tip.label) * 0.2),  # 根据物种数量动态调整高度
    dpi = 300,  # 高清晰度
    bg = "white",  # 白色背景
    limitsize = FALSE  # 允许生成大尺寸图片
  )
  
  cat("\n高清晰度PNG图片已生成: ", output_file, "\n")
  
  # 返回绘图对象，便于显示
  return(p2)
}

# 如果有命令行参数，则使用参数指定的基因名称
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  p2 <- run_with_gene(args[1])
  if (!is.logical(p2)) {
    print(p2)
  }
}
