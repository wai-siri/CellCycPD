###PANCER 基因表达小提琴图和箱式图
# 读取基因列表
gene <- read.delim("gene.txt", header=F, stringsAsFactors=F)$V1

# 加载TCGA和GTEx的组合数据
load("self/TCGA_GTEx_pancancer_mrna_pheno.rdata")

# 创建新的输出目录
output_dir <- "output/new"
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(latex2exp)
library(rstatix)
library(RColorBrewer)

gene <- gsub(" ","",gene)

scientific_colors <- c("#2166AC", "#B2182B")  # 蓝色和红色，分别代表Normal和Tumor
names(scientific_colors) <- c("Normal", "Tumor")

# 设置面板背景色
panel_bg_color <- "#E6F3FF"  

# 自定义wilcox检验函数
my_wilcox_test <- function(data) {
  # 检查每组至少有3个观测值
  n1 <- sum(data$Group == "Normal")
  n2 <- sum(data$Group == "Tumor")
  
  if (n1 >= 3 && n2 >= 3) {
    test <- wilcox.test(gene ~ Group, data = data)
    return(test$p.value)
  } else {
    return(NA)
  }
}

for (i in gene) {
  print(i)
  plot_df <- tcga_gtex_mrna_pheno %>%
    select(1:4,all_of(i)) %>%
    filter(sample_type %in% c("GTEx_normal","TCGA_normal","TCGA_tumor"))
  
  colnames(plot_df)[5] <- "gene"
  head(plot_df)
  # 按中位数由高到低排列
  data_new <- plot_df %>% 
    group_by(project) %>% 
    mutate(median = median(gene), group_max = max(gene)) %>% 
    arrange(desc(median))
  data_new$sample_type <- ifelse(data_new$sample_type != "TCGA_tumor","Normal",data_new$sample_type )
  data_new$sample_type <- ifelse(data_new$sample_type == "TCGA_tumor","Tumor",data_new$sample_type)
  data_new$gene <- ifelse(data_new$gene < -5,5, data_new$gene)
  # 调整因子顺序
  data_new$project <- factor(data_new$project, levels = unique(data_new$project))
  
  data_new <- data.frame(data_new)
  data_new$sample_type <- as.factor(data_new$sample_type)
  
  head(data_new)
  colnames(data_new)[4] <- "Group"
  
  data_new[1:6,]
  data_new<- as.data.frame(data_new)
  
  if (!is.factor(data_new$Group))  { 
    data_new$Group <- as.factor(data_new$Group)  
  } 
  
  # 计算每个project的p值
  p_values <- data_new %>%
    group_by(project) %>%
    summarise(
      p.value = my_wilcox_test(pick(everything())),  # 替换cur_data()以避免警告
      y.position = max(gene) * 1.1,  # 调整高度使其在顶部显示
      x.position = 1  # 将x位置设在中央
    )
  
  # 调整p值
  p_values$p.adj <- p.adjust(p_values$p.value, method = "bonferroni")
  
  # 创建标签，显示所有癌症的具体p值
  p_values$p.label <- sapply(p_values$p.adj, function(p) {
    if (is.na(p)) return("p = NA")
    # 格式化p值，使用科学计数法显示非常小的p值
    if (p < 0.001) {
      return(sprintf("p = %.2e", p))
    } else {
      return(sprintf("p = %.4f", p))
    }
  })
  
  p <- ggplot(data = data_new, aes(x = project, y = gene)) + 
    # 小提琴图设置
    geom_violin(aes(fill = Group), alpha = 0.7) +
    # 箱线图设置
    geom_boxplot(aes(colour = Group), fill="white", 
                outlier.shape = NA, width = 0.2,
                outlier.color = "black",
                position = position_dodge(0.9)) + 
    # 颜色设置
    scale_fill_manual(values = scientific_colors) +
    scale_color_manual(values = c("black", "black")) +
    # 主题设置
    theme_bw() + 
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 10),
      legend.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 12, face = "bold.italic"),
      axis.title.y = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 10, face = "bold", color = "black"),
      strip.background = element_rect(fill = panel_bg_color, color = "black"),  # 添加边框
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_blank()
    ) + 
    # 标签设置
    ylab(paste0(i, " expression: log2(TPM+0.01)")) + 
    xlab("Cancer") + 
    # 分面设置
    facet_wrap(. ~ project, scales = "free_x", nrow = 3) + 
    scale_x_discrete(breaks = unique(data_new$project))
  
  # 添加p值标签，在每个的顶部中央显示
  p1 <- p + geom_text(
    data = p_values,
    aes(x = x.position, y = y.position, label = p.label),
    size = 3,
    hjust = 0.5,  # 水平居中对齐
    vjust = 0,    # 垂直对齐到顶部
    fontface = "bold"  # 加粗显示
  )
  
  # 保存图片
  ggsave(file.path(output_dir, paste0(i, ".png")), 
         plot = p1, width = 12, height = 8, dpi = 400)
  
  print(paste("Completed:", i))
}
