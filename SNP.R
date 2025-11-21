# 加载所需包  
library(maftools)  # MAF 可视化  
library(data.table)  # 高效数据操作  
library(pbapply)

# 加载目标基因文件  
gene_file <- "D:/R_Project/intersection.csv"  
gene_info <- fread(gene_file)  

# 提取基因符号和基因范围  
gene_symbols <- gene_info[[1]]  # 第一列是基因符号  
gene_ranges <- gene_info[, .(  
  chromosome_name = gene_info[[6]],  # 第6列：染色体名  
  start_position = gene_info[[7]],  # 第7列：起始位置  
  end_position = gene_info[[8]]     # 第8列：终止位置  
)]  

# 获取癌症类型数据  
data_dir <- "D:/R_Project/TCGA_SNP_Data/"  
data_files <- list.files(data_dir, pattern = "_SNP\\.Rdata$", full.names = TRUE)  


# 设置保存路径  
output_dir <- "D:/R_Project/Visualization_Results/"  
dir.create(output_dir, showWarnings = FALSE)  


# 主函数：处理癌症类型数据并完成可视化  
process_cancer_data <- function(data_file) {
  load(data_file)
  cancer_type <- gsub(pattern = "_SNP\\.Rdata$", replacement = "", x = basename(data_file))
  message(paste("正在处理癌症类型：", cancer_type))
  
  maf_data <- data
  if (!"Hugo_Symbol" %in% colnames(maf_data)) {
    stop(paste("Error: 未找到 'Hugo_Symbol' 列，请确认数据列名。"))
  }
  
  # 筛选目标基因突变
  target_mutations <- maf_data[maf_data$Hugo_Symbol %in% gene_symbols, ]
  for (i in 1:nrow(gene_ranges)) {
    current_range <- gene_ranges[i, ]
    filtered <- maf_data[
      maf_data$Chromosome == current_range$chromosome_name &
        maf_data$Start_Position >= current_range$start_position &
        maf_data$End_Position <= current_range$end_position,
    ]
    target_mutations <- rbind(target_mutations, filtered)
  }
  target_mutations <- unique(target_mutations)
  
  if (nrow(target_mutations) == 0) {
    message(paste("癌症类型", cancer_type, "未发现符合的目标突变"))
    return(NULL)
  }
  
  maf <- read.maf(maf = target_mutations)
  
  png(file.path(output_dir, paste0(cancer_type, "_Summary_Plot.png")), width = 2400, height = 1800, res = 300)
  plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
  dev.off()
  
  png(file.path(output_dir, paste0(cancer_type, "_Oncoplot.png")), width = 2400, height = 1800, res = 300)
  oncoplot(maf = maf, top = 20, draw_titv = TRUE)
  dev.off()
  
  titv_res <- titv(maf = maf, plot = FALSE)
  png(file.path(output_dir, paste0(cancer_type, "_TiTv_Plot.png")), width = 2400, height = 1800, res = 300)
  plotTiTv(res = titv_res)
  dev.off()
  
  png(file.path(output_dir, paste0(cancer_type, "_Gene_Interactions.png")), width = 2400, height = 1800, res = 300)
  somaticInteractions(maf = maf, top = 20, pvalue = c(0.05, 0.1))
  dev.off()
  
  png(file.path(output_dir, paste0(cancer_type, "_Rainfall_Plot.png")), width = 2400, height = 1800, res = 300)
  rainfallPlot(maf = maf, detectChangePoints = TRUE, pointSize = 0.4)
  dev.off()
  
  message(paste("已完成癌症类型的分析：", cancer_type))
  return(cancer_type)
}

# 处理所有癌症数据  
results <- pblapply(data_files, process_cancer_data)

# 输出结果  
processed_cancers <- unlist(results)  
message("以下癌症类型的数据分析已完成：")  
print(processed_cancers)