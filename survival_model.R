# 基于GSVA通路评分结果构建简单的Cox比例风险回归模型

library(survival)
library(survminer)
library(dplyr)
library(tidyr)

# 设置工作目录
# setwd("d:/R_Project/cellcycle/code/tool/GSVA")

# 1. 读取GSVA评分矩阵
cat("读取GSVA评分矩阵...\n")
gsva_scores <- read.table("D:/R_Project/cellcycle/code/tool/GSVA/GSVA_Scores_Matrix.txt", 
                          header = TRUE, sep = "\t", check.names = FALSE)

# 2. 读取数据
cat("读取数据...\n")

# 读取GSVA评分矩阵
# 由于GSVA_Scores_Matrix.txt的第一列是样本ID，但没有列名，我们需要特殊处理
gsva_data <- read.table("D:/R_Project/cellcycle/code/tool/GSVA/GSVA_Scores_Matrix.txt", header = TRUE, sep = "\t", check.names = FALSE)

# 创建正确的数据框，第一列为样本ID，第二列为CellCycle_Genes评分
gsva_scores <- data.frame(
  sample = rownames(gsva_data),
  CellCycle_Genes = gsva_data$CellCycle_Genes,
  stringsAsFactors = FALSE
)

# 读取临床数据
clinical_data <- read.table("D:/R_Project/cellcycle/code/tool/GSVA/BRCA_clinicalMatrix", header = TRUE, sep = "\t", check.names = FALSE)

# 输出数据的基本信息
cat("GSVA评分数据维度:", dim(gsva_scores), "\n")
cat("临床数据维度:", dim(clinical_data), "\n")
cat("GSVA评分数据列名:", paste(colnames(gsva_scores)[seq_len(min(3, ncol(gsva_scores)))], collapse=", "), if(ncol(gsva_scores)>3) "..." else "", "\n")
cat("临床数据前10列列名:", paste(colnames(clinical_data)[seq_len(min(10, ncol(clinical_data)))], collapse=", "), "...\n")

# 3. 读取样本注释信息
cat("读取样本注释信息...\n")
sample_info <- read.table("D:/R_Project/cellcycle/code/tool/GSVA/TCGA-BRCA.sample.info.tsv", 
                         header = TRUE, sep = "\t", check.names = FALSE)

# 4. 数据预处理
cat("处理数据...\n")

# 检查GSVA评分矩阵的结构
cat("GSVA评分矩阵维度:", dim(gsva_scores), "\n")
cat("GSVA评分矩阵列名:", colnames(gsva_scores), "\n")

# 处理GSVA评分矩阵
# 我们已经创建了正确的数据框，确保sample列存在
if(!"sample" %in% colnames(gsva_scores)) {
  stop("GSVA评分矩阵中找不到sample列")
}

# 统一样本ID格式 - 将GSVA评分中的点替换为连字符
gsva_scores$sample <- gsub("\\.", "-", gsva_scores$sample)

# 首先检查临床数据的列名
cat("临床数据列名：", paste(colnames(clinical_data)[seq_len(min(10, ncol(clinical_data)))], collapse=", "), "...\n")

# 提取临床数据中的生存信息
clinical_subset <- clinical_data

# 检查必要的列是否存在
if("sampleID" %in% colnames(clinical_data)) {
  clinical_subset <- clinical_subset %>% rename(sample = sampleID)
} else if("sample" %in% colnames(clinical_data)) {
  # 已经有sample列，不需要重命名
} else {
  # 如果没有样本ID列，使用第一列作为样本ID
  clinical_subset <- clinical_subset %>% 
    mutate(sample = clinical_subset[[1]]) %>%
    select(sample, everything())
}

# 检查并重命名必要的生存列
if("OS_Time_nature2012" %in% colnames(clinical_subset)) {
  clinical_subset <- clinical_subset %>% rename(OS_time = OS_Time_nature2012)
} else if(!"OS_time" %in% colnames(clinical_subset)) {
  stop("找不到生存时间列")
}

if("OS_event_nature2012" %in% colnames(clinical_subset)) {
  clinical_subset <- clinical_subset %>% rename(OS_event = OS_event_nature2012)
} else if(!"OS_event" %in% colnames(clinical_subset)) {
  stop("找不到生存事件列")
}

# 重命名其他可能的列
if("Age_at_Initial_Pathologic_Diagnosis_nature2012" %in% colnames(clinical_subset)) {
  clinical_subset <- clinical_subset %>% rename(Age = Age_at_Initial_Pathologic_Diagnosis_nature2012)
}
if("pathologic_stage" %in% colnames(clinical_subset)) {
  clinical_subset <- clinical_subset %>% rename(Stage = pathologic_stage)
}
if("pathologic_T" %in% colnames(clinical_subset)) {
  clinical_subset <- clinical_subset %>% rename(T_stage = pathologic_T)
}
if("pathologic_N" %in% colnames(clinical_subset)) {
  clinical_subset <- clinical_subset %>% rename(N_stage = pathologic_N)
}
if("pathologic_M" %in% colnames(clinical_subset)) {
  clinical_subset <- clinical_subset %>% rename(M_stage = pathologic_M)
}
if("ER_Status_nature2012" %in% colnames(clinical_subset)) {
  clinical_subset <- clinical_subset %>% rename(ER_status = ER_Status_nature2012)
}
if("PR_Status_nature2012" %in% colnames(clinical_subset)) {
  clinical_subset <- clinical_subset %>% rename(PR_status = PR_Status_nature2012)
}
if("HER2_Final_Status_nature2012" %in% colnames(clinical_subset)) {
  clinical_subset <- clinical_subset %>% rename(HER2_status = HER2_Final_Status_nature2012)
}

# 只保留样本ID的前15个字符（TCGA-XX-XXXX-XX格式）以匹配临床数据
gsva_scores$sample_short <- substr(gsva_scores$sample, 1, 15)
clinical_subset$sample_short <- substr(clinical_subset$sample, 1, 15)









# 5. 合并数据集
cat("合并数据集...\n")

# 输出合并前的数据集大小
cat("GSVA评分数据集大小:", nrow(gsva_scores), "行\n")
cat("临床数据子集大小:", nrow(clinical_subset), "行\n")

# 合并GSVA评分和临床数据（使用简化的样本ID）
merged_data <- merge(clinical_subset, gsva_scores, by = "sample_short", all = FALSE)

# 确保合并后的数据有一个明确的sample列
# 使用clinical_subset中的sample列作为最终的sample列
if("sample.x" %in% colnames(merged_data)) {
  # 如果合并后有sample.x和sample.y，使用sample.x作为sample
  merged_data <- merged_data %>% 
    rename(sample = sample.x) %>%
    select(-sample.y)
} else if(!"sample" %in% colnames(merged_data)) {
  # 如果没有sample列，使用sample_short作为sample
  merged_data$sample <- merged_data$sample_short
}

# 输出合并后的数据集大小和列名
cat("合并后的数据集大小:", nrow(merged_data), "行\n")
cat("合并后的数据集列名:", paste(colnames(merged_data)[seq_len(min(10, ncol(merged_data)))], collapse=", "), "...\n")

# 6. 数据清洗
cat("清洗数据...\n")
# 检查关键变量的缺失情况
cat("OS_time缺失值数量:", sum(is.na(merged_data$OS_time)), "\n")
cat("OS_event缺失值数量:", sum(is.na(merged_data$OS_event)), "\n")
cat("CellCycle_Genes缺失值数量:", sum(is.na(merged_data$CellCycle_Genes)), "\n")

# 移除缺失生存数据的样本
merged_data <- merged_data %>%
  filter(!is.na(OS_time) & !is.na(OS_event) & !is.na(CellCycle_Genes))

# 输出清洗后的数据集大小
cat("清洗后的数据集大小:", nrow(merged_data), "行\n")

# 将生存事件转换为数值型
merged_data$OS_event <- as.numeric(merged_data$OS_event)

# 确保有足够的事件样本
cat("事件样本数量:", sum(merged_data$OS_event), "\n")







# 7. 构建Cox比例风险回归模型
cat("构建Cox比例风险回归模型...\n")

# 检查数据是否足够进行分析
if(nrow(merged_data) < 10 || sum(merged_data$OS_event) < 5) {
  stop("数据量不足或事件数量太少，无法构建可靠的Cox模型")
}

# 尝试不同的模型，从简单到复杂
cat("尝试构建单变量模型...\n")
simple_model <- tryCatch({
  coxph(Surv(OS_time, OS_event) ~ CellCycle_Genes, data = merged_data)
}, error = function(e) {
  cat("单变量模型构建失败:", e$message, "\n")
  return(NULL)
})

# 初始化cox_model变量以避免未定义错误
cox_model <- NULL

# 如果简单模型成功，尝试添加年龄
if(!is.null(simple_model)) {
  cat("尝试添加年龄变量...\n")
  age_model <- tryCatch({
    coxph(Surv(OS_time, OS_event) ~ CellCycle_Genes + Age, data = merged_data)
  }, error = function(e) {
    cat("添加年龄变量失败:", e$message, "\n")
    return(simple_model)
  })
  
  # 如果添加年龄成功，尝试添加分期
  if(!identical(age_model, simple_model)) {
    cat("尝试添加分期变量...\n")
    full_model <- tryCatch({
      coxph(Surv(OS_time, OS_event) ~ CellCycle_Genes + Age + Stage, data = merged_data)
    }, error = function(e) {
      cat("添加分期变量失败:", e$message, "\n")
      return(age_model)
    })
    cox_model <- full_model
  } else {
    cox_model <- age_model
  }
} else {
  stop("无法构建任何Cox模型，请检查数据")
}

# 确保cox_model已定义
if(is.null(cox_model)) {
  stop("无法构建 Cox 模型，请检查数据")
}

# 8. 输出模型结果
cat("输出模型结果...\n")
model_summary <- summary(cox_model)
print(model_summary)

# 9. 保存模型结果
cat("保存模型结果...\n")
sink("D:/R_Project/cellcycle/code/tool/GSVA/cox_model_results.txt")
print(model_summary)
sink()

# 10. 可视化生存曲线
cat("绘制生存曲线...\n")

# 根据CellCycle_Genes评分的中位数将样本分为高低两组
median_score <- median(merged_data$CellCycle_Genes, na.rm = TRUE)
merged_data$score_group <- ifelse(merged_data$CellCycle_Genes > median_score, "High", "Low")

# 绘制Kaplan-Meier生存曲线
km_fit <- survfit(Surv(OS_time, OS_event) ~ score_group, data = merged_data)

# 保存生存曲线图
pdf("D:/R_Project/cellcycle/code/tool/GSVA/survival_curve.pdf", width = 8, height = 6)
ggsurvplot(
  km_fit,
  data = merged_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  linetype = "strata",
  surv.median.line = "hv",
  ggtheme = theme_bw(),
  palette = c("#E7B800", "#2E9FDF"),
  title = "Kaplan-Meier Survival Curve by Cell Cycle Pathway Score",
  xlab = "Time (days)",
  legend.labs = c("High Score", "Low Score"),
  risk.table.height = 0.25
)
dev.off()

# 11. 构建模型
cat("构建预测模型...\n")

# 计算风险分数
merged_data$risk_score <- predict(cox_model, type = "risk", newdata = merged_data)

# 根据风险分数的中位数将样本分为高低风险组
median_risk <- median(merged_data$risk_score, na.rm = TRUE)
merged_data$risk_group <- ifelse(merged_data$risk_score > median_risk, "High Risk", "Low Risk")

# 绘制基于风险分数的生存曲线
risk_fit <- survfit(Surv(OS_time, OS_event) ~ risk_group, data = merged_data)

# 保存风险分组生存曲线图
pdf("D:/R_Project/cellcycle/code/tool/GSVA/risk_group_survival.pdf", width = 8, height = 6)
ggsurvplot(
  risk_fit,
  data = merged_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  linetype = "strata",
  surv.median.line = "hv",
  ggtheme = theme_bw(),
  palette = c("#E7B800", "#2E9FDF"),
  title = "Kaplan-Meier Survival Curve by Risk Group",
  xlab = "Time (days)",
  legend.labs = c("High Risk", "Low Risk"),
  risk.table.height = 0.25
)
dev.off()

# 12. 保存风险评分结果
cat("保存风险评分结果...\n")

# 检查必要列是否存在
required_cols <- c("sample", "OS_time", "OS_event", "CellCycle_Genes", "risk_score", "risk_group")
optional_cols <- c("Age", "Stage")

# 检查必要列
for (col in required_cols) {
  if (!(col %in% colnames(merged_data))) {
    stop(paste("缺失必要列:", col))
  }
}

# 准备选择的列
select_cols <- required_cols

# 添加可选列（如果存在）
for (col in optional_cols) {
  if (col %in% colnames(merged_data)) {
    select_cols <- c(select_cols, col)
  } else {
    cat(paste("警告: 列", col, "不存在，将被忽略\n"))
  }
}

# 输出最终选择的列
cat("选择的列:", paste(select_cols, collapse=", "), "\n")

# 选择列并保存
risk_results <- merged_data %>%
  select(all_of(select_cols))

write.table(risk_results, 
            file = "D:/R_Project/cellcycle/code/tool/GSVA/risk_prediction_results.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 13. 模型验证 - C-index计算
cat("计算模型C-index...\n")
c_index <- concordance(cox_model)$concordance
cat(paste("模型C-index:", round(c_index, 3), "\n"))

# 14. 创建一个简单的预测函数
cat("创建预测函数...\n")
predict_risk <- function(new_data, model = cox_model, median_risk = median_risk) {
  # 计算风险分数
  risk_scores <- predict(model, type = "risk", newdata = new_data)
  
  # 分配风险组
  risk_groups <- ifelse(risk_scores > median_risk, "High Risk", "Low Risk")
  
  # 返回结果
  return(data.frame(
    sample = new_data$sample,
    risk_score = risk_scores,
    risk_group = risk_groups
  ))
}

cat("预后模型构建完成！\n")
