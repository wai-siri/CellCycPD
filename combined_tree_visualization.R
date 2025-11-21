
# 加载必要的包
library(Biostrings)
library(DECIPHER)
library(ggtree)
library(treeio)
library(ggsci)
library(ggtreeExtra)
library(ggnewscale)
library(ggplot2)
library(stringr)
library(tidytree)
library(ape)

# 设置工作目录
setwd("d:/R_Project/cellcycle/code/tool/tree")

# 函数：从命令行参数或默认值获取基因名称
get_gene_name <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    return(args[1])
  } else {
    return("CDK1")  # 默认基因名称
  }
}

# 获取基因名称
gene_name <- get_gene_name()

# 构建文件路径
tree_file <- paste0("nwk_data/", gene_name, ".nwk")
msa_file <- paste0("fas_data/", gene_name, "_CDS_all.fas")
orig_file <- paste0("data1/", gene_name, "_CDS_all.fasta")

# 检查文件是否存在
if (!file.exists(tree_file) || !file.exists(msa_file) || !file.exists(orig_file)) {
  stop(paste("错误: 无法找到必要的文件。请确认以下文件存在:\n", 
             tree_file, "\n", msa_file, "\n", orig_file))
}

# 读取数据
tree <- read.tree(tree_file)
seqs <- readDNAStringSet(msa_file)
seqs_orig <- readDNAStringSet(orig_file)

# 处理序列名称中的空格
names(seqs) <- gsub(" ", "_", names(seqs))

# 确保树和序列中的物种匹配
tree_taxa <- tree$tip.label
msa_taxa <- names(seqs)
common_taxa <- intersect(tree_taxa, msa_taxa)

# 如果有不匹配的物种，输出警告信息并过滤
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
  seqs <- seqs[names(seqs) %in% common_taxa]
  
  cat("  已自动过滤不匹配的物种\n\n")
}


# 计算与人类参考序列的保守率
calculate_conservation_rate <- function(aligned_seqs) {
  # 获取比对长度和密码子数量
  align_length <- width(aligned_seqs)[1]
  codon_count <- floor(align_length/3)
  
  # 初始化氨基酸矩阵
  aa_matrix <- matrix("", nrow=length(aligned_seqs), ncol=codon_count)
  
  # 将DNA序列转换为氨基酸序列
  for(i in seq_along(aligned_seqs)) {
    seq_str <- as.character(aligned_seqs[[i]])
    for(j in seq_len(codon_count)) {
      pos <- (j-1)*3 + 1
      if(pos+2 <= align_length) {
        codon <- substr(seq_str, pos, pos+2)
        if(grepl("-", codon) || !all(strsplit(codon, "")[[1]] %in% c("A", "C", "G", "T"))) {
          aa_matrix[i,j] <- "-"
        } else {
          aa_matrix[i,j] <- as.character(translate(DNAString(codon)))
        }
      } else {
        aa_matrix[i,j] <- "-"
      }
    }
  }
  
  # 找到人类参考序列
  human_idx <- which(names(aligned_seqs) == "Homo_sapiens")
  if(length(human_idx) == 0) stop("无法找到人类参考序列")
  
  # 获取人类参考氨基酸序列
  human_aa <- aa_matrix[human_idx,]
  conservation_rates <- numeric(nrow(aa_matrix))
  
  # 计算保守率
  for(i in seq_len(nrow(aa_matrix))) {
    current_aa <- aa_matrix[i,]
    valid_positions <- which(human_aa != "-" & current_aa != "-")
    total_valid <- length(valid_positions)
    
    conservation_rates[i] <- if(total_valid > 0) {
      sum(human_aa[valid_positions] == current_aa[valid_positions]) / total_valid
    } else {
      0
    }
  }
  
  return(conservation_rates)
}

# 计算每个序列与人类参考序列的保守率
seq_conservation <- calculate_conservation_rate(seqs)

# 计算序列信息
seq_info <- data.frame(
  label = gsub(" ", "_", names(seqs)),
  gc_content = rowSums(letterFrequency(seqs, letters=c("G", "C"))) / width(seqs),
  conservation = seq_conservation,
  length = width(seqs_orig)
)

# 定义物种分类
taxonomy <- data.frame(
  label = c(
    # Mammals
    "Balaena_mysticetus", "Heterocephalus_glaber", "Loxodonta_africana",
    "Myotis_brandtii", "Homo_sapiens", "Mus_musculus", "Macaca_mulatta",
    "Canis_lupus_familiaris", "Rattus_norvegicus", "Macropus_eugenii",
    "Bos_taurus", "Tursiops_truncatus", "Equus_caballus", "Felis_catus",
    "Sus_scrofa", "Oryctolagus_cuniculus", "Peromyscus_maniculatus",
    "Saimiri_sciureus", "Gorilla_gorilla", "Equus_asinus", "Ochotona_curzoniae",
    "Erinaceus_europaeus", "Monodelphis_domestica",
    # Birds
    "Gallus_gallus", "Taeniopygia_guttata", "Nipponia_nippon",
    "Falco_peregrinus", "Melopsittacus_undulatus", "Aptenodytes_forsteri",
    "Struthio_camelus", "Haemorhous_mexicanus", "Pezoporus_wallicus",
    "Astur_gentilis", "Cuculus_canorus", "Tyto_alba",
    # Reptiles
    "Chelonia_mydas", "Alligator_mississippiensis", "Chamaeleo_calyptratus",
    "Python_bivittatus", "Anolis_carolinensis", "Crocodylus_porosus",
    # Amphibians
    "Xenopus_laevis", "Ambystoma_mexicanum", "Bombina_orientalis",
    "Rana_temporaria", "Hyla_cinerea",
    # Osteichthyes
    "Danio_rerio", "Oryzias_latipes", "Gasterosteus_aculeatus",
    "Oreochromis_niloticus", "Salmo_salar", "Takifugu_rubripes",
    "Carassius_auratus", "Tetraodon_nigroviridis", "Mugil_cephalus",
    "Paralichthys_olivaceus", "Ictalurus_punctatus",
    # Chondrichthyes
    "Somniosus_microcephalus", "Stegostoma_fasciatum", "Squalus_acanthias",
    "Leucoraja_erinacea", "Rhincodon_typus", "Sphyrna_lewini", "Callorhinchus_milii",
    # Outgroup
    "Ciona_intestinalis", "Branchiostoma_floridae", "Strongylocentrotus_purpuratus"
  ),
  Class = c(
    rep("Mammals", 23),
    rep("Birds", 12),
    rep("Reptiles", 6),
    rep("Amphibians", 5),
    rep("Osteichthyes", 11),
    rep("Chondrichthyes", 7),
    rep("Outgroup", 3)
  )
)

# 计算引导值点的位置
node_info <- fortify(tree)
bootstrap_data <- node_info[!is.na(as.numeric(node_info$label)) & as.numeric(node_info$label) > 0.8,]

# 计算每个叶子节点的角度
leaf_info <- node_info[node_info$isTip,]

# 过滤taxonomy数据框，只保留存在于树中的物种
taxonomy <- taxonomy[taxonomy$label %in% tree$tip.label,]

# 计算文本角度
taxonomy$text_angle <- sapply(taxonomy$label, function(lbl) {
  node <- leaf_info[leaf_info$label == lbl,]
  if(nrow(node) > 0) {
    angle <- (atan2(node$y, node$x) * 180/pi) %% 360
    if(angle > 180) angle <- angle - 180
    return(angle)
  }
  return(0)
})

# 创建输出目录
output_dir <- "res/png"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#---------- 第一部分：创建圆形树 ----------#

# 创建树的基本可视化
p_circular <- ggtree(tree, 
                     layout="circular",   # 使用圆形布局
                     branch.length="none") +  # 使用统一树枝长度，不根据branch.length绘制
  
  geom_tree(aes(color=branch.length), size=0.8, lineend="round", linejoin="round", alpha=0.9) +  # 颜色映射要用树枝长度
  scale_color_viridis_c(name="Evolutionary\nDistance", 
                       option="turbo",   
                       begin=0.2, end=0.8,  # 调整颜色范围
                       guide=guide_colorbar(title.position="top")) +
  # 添加引导值点
  geom_nodepoint(aes(subset=!isTip & as.numeric(label) > 0.8, size=as.numeric(label)),
                color="#333333",
                alpha=0.8) +  # 增加透明度
  # 调整点的大小范围
  scale_size_continuous(name="Bootstrap Value",
                        range=c(0.8, 2.5),  # 点的大小范围
                        limits=c(0.8, 1.0)) +
  # 物种分类圈层
  new_scale_fill() +
  geom_fruit(data=taxonomy,
             geom=geom_tile,
             mapping=aes(y=label, fill=Class),
             width=12,     # 分类圈层宽度
             offset=0.5) +  # 向外移动
  scale_fill_npg(alpha=0.5) +  # 透明
  # 物种名称注释在分类圈层上
  geom_fruit(data=taxonomy,
             geom=geom_text,
             mapping=aes(y=label, label=label, angle=text_angle),
             size=2.5,
             offset=0.05,     # 与分类圈层使用相同的offset
             hjust=0.5,     # 文本居中
             vjust=0.5) +   # 垂直居中
  # GC含量圈层
  new_scale_fill() +
  geom_fruit(data=seq_info,
             geom=geom_tile,
             mapping=aes(y=label, fill=gc_content),
             width=2,     # GC含量圈层宽度
             offset=0.55) +  # 分类圈层外侧
  scale_fill_viridis_c(name="GC Content",
                       option="viridis") +
  
  # 保守率层
  new_scale_fill() +
  geom_fruit(data=seq_info,
             geom=geom_tile,
             mapping=aes(y=label, fill=conservation),
             width=3,
             offset=0.3) +  # 外移距离
  scale_fill_gradient(name="Conservation\nRate",
                      low="yellow", high="red",
                      limits=c(0, 1),
                      guide=guide_colorbar(title.position="top")) +
  # 保守率数值标签
  geom_fruit(data=seq_info,
             geom=geom_text,
             mapping=aes(y=label, label=sprintf("%.2f", conservation)),
             size=3.0,  # 字体大小
             offset=0,  # 保守率圈层相同的偏移量
             color="#ffffff",  
             fontface="bold") +
  # 树形图主题
  theme_tree2() +
  # 去掉坐标轴和背景
  theme_void() +
  # 图例位置和样式
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(6, 6, 6, 6),
    legend.box.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(fill="white", color=NA),
    legend.title = element_text(size=10, face="bold"),
    legend.text = element_text(size=8),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    plot.background = element_rect(fill="white", color=NA)
  )

#---------- 第二部分：创建矩形树+多序列比对 ----------#

# 绘制矩形进化树
p_rectangular <- ggtree(tree, branch.length = "none", layout = "rectangular") + 
  geom_tiplab(size = 3, color = "#3366CC", fontface = "bold", align = TRUE) +  
  theme_tree() +  
  theme(legend.position = "none")  # 去除图例

p_rectangular <- p_rectangular + 
  geom_tree(aes(color = branch.length), size = 1) +  
  scale_color_viridis_c(option = "turbo", begin = 0.2, end = 0.7) +  
  # 添加bootstrap值标签
  geom_nodelab(aes(label=label, subset=!isTip & as.numeric(label) > 0.1), 
              size=2.5, color="#333333", hjust=1.2, vjust=-0.2)  

# 将DNAStringSet转换为ape包可以处理的DNAbin格式
try({
  dna_bin <- as.DNAbin(seqs)
  
  # 然后使用msaplot将多序列比对与树形合并
  # 增加offset参数值，增加树和多序列比对之间的间距
  p_msa <- msaplot(p_rectangular, dna_bin, offset = 8, width = 2)
}, silent = FALSE)

p_msa <- p_msa + 
  # 移除底部数轴标记
  theme(axis.text.x = element_blank(),  # 移除x轴文本
        axis.ticks.x = element_blank(),  # 移除x轴刻度线
        axis.line.x = element_blank(),   # 移除x轴线
        axis.title.x = element_blank())  # 调整标题样式

#---------- 第三部分：分别保存圆形树和矩形树 ----------#

# 保存圆形树 PDF 格式
ggsave(
  filename = file.path(output_dir, paste0(gene_name, "_circular_tree.pdf")),
  plot = p_circular,
  width = 12,
  height = 12,
  limitsize = FALSE
)

# 保存圆形树 PNG 格式
ggsave(
  filename = file.path(output_dir, paste0(gene_name, "_circular_tree.png")),
  plot = p_circular,
  width = 12,
  height = 12,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

# 保存矩形树 PDF 格式
ggsave(
  filename = file.path(output_dir, paste0(gene_name, "_rectangular_tree_msa.pdf")),
  plot = p_msa,
  width = 15,
  height = 12,
  limitsize = FALSE
)

# 保存矩形树 PNG 格式
ggsave(
  filename = file.path(output_dir, paste0(gene_name, "_rectangular_tree_msa.png")),
  plot = p_msa,
  width = 15,
  height = 12,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

cat("\n\n完成！已分别生成以下文件:\n")
cat(file.path(output_dir, paste0(gene_name, "_circular_tree.pdf")), "\n")
cat(file.path(output_dir, paste0(gene_name, "_circular_tree.png")), "\n")
cat(file.path(output_dir, paste0(gene_name, "_rectangular_tree_msa.pdf")), "\n")
cat(file.path(output_dir, paste0(gene_name, "_rectangular_tree_msa.png")), "\n")

# 显示圆形树
print(p_circular)

# 显示矩形树
print(p_msa)
