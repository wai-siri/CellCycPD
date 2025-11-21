# 设置工作目录
setwd("d:/R_Project/cellcycle/code/tool/tree")

# 加载必要的包
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("treeio")
BiocManager::install("Biostrings")
BiocManager::install("ggsci")  # 美化配色
BiocManager::install("ggtreeExtra")  # 添加外圈注释

# 设置CRAN镜像
repos <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
if (!requireNamespace("ggnewscale", quietly=TRUE))
  install.packages("ggnewscale", repos=repos)  # 添加多个颜色映射

library(ggtree)
library(treeio)
library(ggplot2)
library(Biostrings)
library(ggsci)
library(ggtreeExtra)
library(ggnewscale)
library(DECIPHER)
library(stringr)
library(tidytree)

# 加载BLOSUM62矩阵
data("BLOSUM62", package = "Biostrings")
# 创建本地副本以避免绑定问题
blosum62_matrix <- BLOSUM62

# 读取树文件
tree <- read.tree("nwk_data/ABL1.nwk")

# 读取序列数据
seqs <- readDNAStringSet("fas_data/ABL1_CDS_all.fas")
seqs_orig <- readDNAStringSet("data1/ABL1_CDS_all.fasta")

# 加载需要的包
if (!requireNamespace("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")
if (!requireNamespace("DECIPHER", quietly = TRUE))
  BiocManager::install("DECIPHER")

library(Biostrings)
library(DECIPHER)

# 定义BLOSUM62矩阵 - 用于氨基酸保守性评分
data(BLOSUM62)

# 将DNA序列翻译成蛋白质序列
translate_dna <- function(dna_seq) {
  # 移除gap后再翻译
  clean_seq <- gsub("-", "", as.character(dna_seq))
  # 确保长度是3的倍数
  remainder <- nchar(clean_seq) %% 3
  if(remainder > 0) {
    clean_seq <- substr(clean_seq, 1, nchar(clean_seq) - remainder)
  }
  # 翻译成蛋白质序列
  if(nchar(clean_seq) >= 3) {
    prot_seq <- suppressWarnings(translate(DNAString(clean_seq)))
    return(as.character(prot_seq))
  } else {
    return("")
  }
}

# 计算位点级别的保守性得分
calculate_position_conservation <- function(aligned_seqs) {
  # 获取比对长度
  align_length <- width(aligned_seqs)[1]
  
  # 将DNA序列转换为氨基酸序列矩阵
  aa_matrix <- matrix("", nrow=length(aligned_seqs), ncol=align_length/3)
  
  # 对每个序列，按密码子(3个碱基)提取并翻译
  for(i in seq_along(aligned_seqs)) {
    seq_str <- as.character(aligned_seqs[[i]])
    for(j in 1:(align_length/3)) {
      pos <- (j-1)*3 + 1
      codon <- substr(seq_str, pos, pos+2)
      # 如果密码子包含gap或非标准碱基，标记为gap
      if(grepl("-", codon) || !all(strsplit(codon, "")[[1]] %in% c("A", "C", "G", "T"))) {
        aa_matrix[i,j] <- "-"
      } else {
        aa <- translate(DNAString(codon))
        aa_matrix[i,j] <- as.character(aa)
      }
    }
  }
  
  # 初始化保守性得分向量
  cons_scores <- numeric(ncol(aa_matrix))
  
  # 对每个位置计算保守性得分
  for(pos in seq_len(ncol(aa_matrix))) {
    # 提取该位置的所有氨基酸
    aas <- aa_matrix[,pos]
    # 去除gap
    valid_aas <- aas[aas != "-"]
    
    if(length(valid_aas) <= 1) {
      # 如果只有一个或没有有效氨基酸，设为0
      cons_scores[pos] <- 0
    } else {
      # 计算该位置的平均BLOSUM62得分
      total_score <- 0
      pair_count <- 0
      
      for(i in seq_len(max(0, length(valid_aas)-1))) {
        for(j in seq(i+1, length(valid_aas))) {
          aa1 <- valid_aas[i]
          aa2 <- valid_aas[j]
          
          # 获取BLOSUM62得分
          if(aa1 %in% rownames(blosum62_matrix) && aa2 %in% colnames(blosum62_matrix)) {
            score <- blosum62_matrix[aa1, aa2]
            total_score <- total_score + score
            pair_count <- pair_count + 1
          }
        }
      }
      
      if(pair_count > 0) {
        # 平均得分
        avg_score <- total_score / pair_count
        # 归一化到0-1范围 (BLOSUM62得分范围大约从-4到11)
        norm_score <- (avg_score + 4) / 15
        cons_scores[pos] <- norm_score
      } else {
        cons_scores[pos] <- 0
      }
    }
  }
  
  return(cons_scores)
}

# 计算每个序列的平均保守性得分
position_cons_scores <- calculate_position_conservation(seqs)

# 为每个序列计算平均保守性得分
seq_conservation <- numeric(length(seqs))
for(i in seq_along(seqs)) {
  # 获取序列中的有效位置
  seq_str <- as.character(seqs[[i]])
  valid_positions <- numeric(length(position_cons_scores))
  
  for(j in seq_along(position_cons_scores)) {
    pos <- (j-1)*3 + 1
    codon <- substr(seq_str, pos, pos+2)
    # 如果密码子不包含gap且只有标准碱基，则该位置有效
    if(!grepl("-", codon) && all(strsplit(codon, "")[[1]] %in% c("A", "C", "G", "T"))) {
      valid_positions[j] <- 1
    } else {
      valid_positions[j] <- 0
    }
  }
  
  # 计算该序列的平均保守性得分
  if(sum(valid_positions) > 0) {
    seq_conservation[i] <- sum(position_cons_scores * valid_positions) / sum(valid_positions)
  } else {
    seq_conservation[i] <- 0
  }
}

# 标准化保存性得分
seq_conservation_z <- scale(seq_conservation)[,1]

# 计算序列信息
seq_info <- data.frame(
    label = gsub(" ", "_", names(seqs)),
    gc_content = rowSums(letterFrequency(seqs, letters=c("G", "C"))) / width(seqs),
    seq_length = width(seqs),
    conservation = seq_conservation
)

# 标准化序列长度和保存性得分
seq_info$seq_length_z <- scale(seq_info$seq_length)[,1]
seq_info$conservation_z <- seq_conservation_z

# 计算序列信息
seq_info$gc_content <- rowSums(letterFrequency(seqs, letters=c("G", "C"))) / width(seqs)
seq_info$length <- width(seqs_orig)

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
taxonomy$text_angle <- sapply(seq_len(nrow(taxonomy)), function(i) {
  node <- leaf_info[leaf_info$label == taxonomy$label[i],]
  if(nrow(node) > 0) {
    # 计算角度，使文本沿着树枝方向
    angle <- (atan2(node$y, node$x) * 180/pi) %% 360
    # 对于180-360度的角度，调整文本方向以保持可读性
    if(angle > 180) {
      angle <- angle - 180
    }
    angle
  } else {
    0
  }
})

# 创建树的基本可视化
p <- ggtree(tree, 
            layout="circular",   # 使用圆形布局
            size=0.9,           # 增加树枝线条粗细
            color="darkgray") +  # 设置基本树枪颜色
     # 添加引导值点
     geom_point(data=bootstrap_data,
                aes(x=x, y=y, size=as.numeric(label)),
                color="black",
                alpha=0.5) +
     # 调整点的大小范围
     scale_size_continuous(name="Bootstrap Value",
                          range=c(0.5, 1.5),
                          limits=c(0.8, 1.0)) +
     # 添加物种分类圈层，使用透明度来表示
     new_scale_fill() +
     geom_fruit(data=taxonomy,
                geom=geom_tile,
                mapping=aes(y=label, fill=Class),
                width=0.4,     # 增加分类圈层宽度
                offset=0.38) +  # 向外移动，避免与树重叠
     scale_fill_npg(alpha=0.5) +  # 使用半透明的颜色
     # 添加物种名称注释，直接写在分类圈层上
     geom_fruit(data=taxonomy,
                geom=geom_text,
                mapping=aes(y=label, label=label, angle=text_angle),
                size=2.5,
                offset=0.05,     # 与分类圈层使用相同的offset
                hjust=0.5,     # 文本居中对齐
                vjust=0.5) +   # 垂直居中对齐
     # 添加GC含量圈层
     new_scale_fill() +
     geom_fruit(data=seq_info,
                geom=geom_tile,
                mapping=aes(y=label, fill=gc_content),
                width=0.08,     # 减GC含量圈层宽度
                offset=0.35) +  # 放在分类圈层外侧
     scale_fill_viridis_c(name="GC Content",
                         option="viridis") +

     # 添加保存性得分层
     new_scale_fill() +
     geom_fruit(data=seq_info,
                geom=geom_tile,
                mapping=aes(y=label, fill=conservation_z),
                width=0.08,
                offset=0.3) +  # 紧贴GC含量圈层
     scale_fill_gradient2(name="Conservation\n(Z-score)",
                         low="yellow", mid="orange", high="red",
                         midpoint=0,
                         guide=guide_colorbar(title.position="top"))
     # 使用树形图主题
     theme_tree2() +
     # 去掉坐标轴和背景
     theme_void() +
     # 调整图例位置和样式
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

# 保存图形
ggsave("ABL1_phylogenetic_tree.pdf", p, width=15, height=15, bg="white")
