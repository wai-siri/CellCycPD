######20250506 99 zhifangyinzi fenxing 

setwd("E:/01 papers/2025 paper/19 cell cycle-database by Wenyuan Wang/")
rm(list = ls())
options(stringsAsFactors = F)
library(GSVA)
library(ggplot2)
library(limma)
library(pheatmap)
library(reshape2)
library(mclust)
library(stringr)
# 1. 根据cell cycle sig DEgenes 无监督一致性聚类
gene <- read.table("./gene123.txt",header = F,sep = "\t")
head(gene)

#### 
#BiocManager::install('GSVAdata')
#BiocManager::install('GSVA')
#BiocManager::install("limma")
#BiocManager::install("estimate")

cancer <- read.table("D:/TCGA_20220305/GDC_TCGA/sample.number.txt",header = T,sep = "\t")
cancer.1 <- cancer[cancer$tumor > 0,]

cancer.2 <- sort(cancer.1$cancer)
#cancer.2 <- c("ACC","KIRC","KIRP","UVM")
table <- c()
for (m in cancer.2) {
  
  #m <- cancer[i]
  print(m)
  #m="LIHC"
  n= paste("TCGA-",m,sep = "")
  TCGA_gset <- read.table(paste("D:/TCGA_20220305/GDC_TCGA/mRNA_exp_TPM_only_TCGA/",n,".gene_expression_TPM.tsv",sep = ""),header = T,sep = "\t")
  TCGA_gset <- data.frame(TCGA_gset)
  TCGA_gset = as.data.frame(avereps(TCGA_gset[,-1],ID = TCGA_gset$ID) ) ##相同ID的数据取均值
  colnames(TCGA_gset) <- substring(colnames(TCGA_gset),1,15) %>% gsub("-",".",.)
  TCGA_gset[1:4,1:4]
  TCGA_gset <- log2(TCGA_gset+1)
  
  ##筛选癌症样本
  exp <- TCGA_gset[,str_sub(colnames(TCGA_gset),14,15)<10]
  
  df <- exp[rownames(exp) %in% gene$V1,]
  
  # 2.筛选基因(通过中位数绝对偏差度量，MAD)
  mads <- apply(df,1,mad) # MAD测度
  #df <- df[rev(order(mads))[1:5000],] #提取前5000个基因
  df <- df[rev(order(mads)),]
  
  # 3.标准化
  df <-  sweep(df,1, apply(df,1,median,na.rm=T)) # 在行的方向上减去最小值，默认是减法
  df <- as.matrix(df)
  # 4.运行ConsensusClusterPlus
  library(ConsensusClusterPlus)
  maxK <-  6 # 选一个K值进行尝试
  results <-  ConsensusClusterPlus(df,
                                   maxK = maxK,
                                   reps = 1000,              # 抽样次数(一般1000或更多)
                                   pItem = 0.8,              # 抽样比例
                                   pFeature = 1,
                                   clusterAlg = "pam",       # 聚类方法
                                   distance="pearson",       # 距离计算方法
                                   title=paste("./CLUSTER/",m,"-0909",sep = ""), # 结果保存路径
                                   innerLinkage="complete",  # 这里不建议使用默认的方法"average"
                                   plot="pdf")               # 结果保存形式
  
  
  # 5.用PAC的方法确定最佳聚类数
  #   面积最小值对应K为最佳K
  Kvec = 2:maxK
  x1 = 0.1; x2 = 0.9        # threshold defining the intermediate sub-interval
  PAC = rep(NA,length(Kvec)) 
  names(PAC) = paste("K=",Kvec,sep="")  # from 2 to maxK
  for(i in Kvec){
    M = results[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])          # M 为计算出共识矩阵
    PAC[i-1] = Fn(x2) - Fn(x1)
  } 
  
  
  optK = Kvec[which.max(PAC)]  # 理想的K值
  
  icl = calcICL(results,
                title=paste("./CLUSTER/",m,"-0909",sep = ""),
                plot="pdf")
  #optK <- 2
  group <- results[[optK]][["consensusClass"]]
  #group <- results[[2]][["consensusClass"]]
  group <- data.frame(group)
  head(group)
  group$sampleID <- rownames(group)
  group$sampleID <- gsub(".03","",group$sampleID)
  group$sampleID <- gsub(".01","",group$sampleID)
  #group$sampleID <- gsub("\\.","-",group$sampleID)
  #cli <- read.table(paste("G:/TCGA_20220305/clinical-data/clinical-data/TCGA-",m,"_clinical.csv",sep = ""),header = T,sep = "\t")  
  cli <- read.table(paste("D:/TCGA_20220305/clinical/TCGA.",m,".sampleMap/",m,"_clinicalMatrix",sep = ""),header = T,sep = "\t")  
  
  
  cli[1:3,]
  colnames(cli)
  cli.1 <- cli[,c("X_PATIENT","vital_status","days_to_death","days_to_last_followup","age_at_initial_pathologic_diagnosis",
                  "gender") ] 
  
  #cli.1 <- cli[,c("sampleID","vital_status","days_to_death","days_to_last_follow_up","age_at_diagnosis",
  #               "gender","tumor_grade") ] 
  
  
  cli.1$days_to_death = with(cli.1,ifelse(is.na(days_to_death),'',days_to_death))
  cli.1$vital_status <- ifelse(cli.1$vital_status == "LIVING",0,1) #重新赋值alive??? 0，dead??? 1
  
  cli.1$time <- ifelse(cli.1$vital_status == 1, cli.1$days_to_death, cli.1$days_to_last_followup)
  
  #cli.1$time <- with(cli.1,ifelse(vital_status == 1, cli.1$days_to_death, cli.1$days_to_last_follow_up))
  cli.1$gender <- ifelse(cli.1$gender=="MALE",1,2)
  head(cli.1)
  cli.1 <- cli.1[,c("X_PATIENT","vital_status","time","gender","age_at_initial_pathologic_diagnosis")]   
  
  colnames(cli.1) <- c( "sampleID","status","time","gender","age")
  cli.1$time <- as.numeric(cli.1$time)
  head(cli.1)
  
  
  cli.1 <- na.omit(cli.1)
  #df[complete.cases(df[, c(1,3)]), ]
  cli.1$sampleID <- gsub("-",".",cli.1$sampleID)
  
  head(group)
  cli.3 <- merge(group,cli.1, by="sampleID")
  cli.3$time <- as.numeric(cli.3$time)
  
  ###survival analysis
  #### survival analysis based on  clusters
  library(survminer)
  library(survival) 
  cli.3 <- cli.3[,1:4]
  head(cli.3)
  cli.3$group <- as.character(cli.3$group)
  cli.3 <- unique(cli.3)
  fit <- survfit(Surv((time), status) ~ group, data = cli.3)
  summary(fit)
  
  pdf(file = paste("./CLUSTER/",m,"-0909/","survival.pdf",sep = ""),onefile = TRUE,width=8,height=8)
  
  p=ggsurvplot(fit,
               data = cli.3,
               pval = TRUE, 
               conf.int = TRUE,
               palette = "hue",
               # linetype = "strata", # Change line type by groups
               surv.median.line = "hv",
               risk.table = TRUE)
  
  
  print(p)
  
  dev.off()
  aa <- surv_pvalue(fit, method = "survdiff")
  
  test <- cbind(cancer=m,group=optK, pval=aa$pval)
  print(test)
  table <- rbind(table,test)
  
}


table1 <- table
table1 <- as.data.frame(table1)
table1 <- table1[order(table1$pval,decreasing = F),]

write.table(table1,"./CLUSTER/123-gene-clustered.pval.txt",quote = F,sep = "\t",row.names = F)


######### 20250909 123 gene pancancer exp heatmap

library(ggplot2)
library(ggthemes)
library(dplyr)

gene <- read.table("./gene123.txt",header = F,sep = "\t")
head(gene)

#### 
#BiocManager::install('GSVAdata')
#BiocManager::install('GSVA')
#BiocManager::install("limma")
#BiocManager::install("estimate")

cancer <- read.table("D:/TCGA_20220305/GDC_TCGA/sample.number.txt",header = T,sep = "\t")
cancer.1 <- cancer[cancer$normal > 4,]

cancer.2 <- sort(cancer.1$cancer)
#cancer.2 <- c("ACC","KIRC","KIRP","UVM")
table <- c()
for (m in cancer.2) {
  
  #m <- cancer[i]
  print(m)
  #m="LIHC"
  n= paste("TCGA-",m,sep = "")
  nrDEG <- read.table(paste("D:/TCGA_20220305/GDC_TCGA/mRNA_DEGs.exp.TPM_only_TCGA/",m,"_TPM_limma_mRNA_DEG.txt",sep = ""),header = T,sep = "\t")
  nrDEG$Gene <- row.names(nrDEG)

  nrDEG$type<-as.factor(ifelse(nrDEG$adj.P.Val<0.05&abs(nrDEG$logFC)>1.5,ifelse(nrDEG$logFC>1.5,'up','down'),'not'))
  nrDEG$cancer <- m
  head(nrDEG)
  nrDEG <- nrDEG[nrDEG$Gene %in% gene$V1,]
  table <- rbind(nrDEG,table)

}
# ??ɽͼ
DEGs <- table
head(DEGs)

DEGs <- na.omit(DEGs)
DEGs$logFC <- ifelse(DEGs$logFC >3,3,DEGs$logFC)
DEGs$logFC <- ifelse(DEGs$logFC < -3,-3,DEGs$logFC)

p <- ggplot(data=DEGs, aes(x=reorder(Gene,logFC,FUN=median), 
                           y=reorder(cancer,logFC,FUN=median)))+
  geom_tile(aes(fill=logFC))+
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "white", "red"))(100),
                       limits = c(-3, 3),  # 可以根据你的数据范围调整 
                       name = "logFC") +
  
 theme(axis.text.x = element_text(angle=90,size = 6))+
  geom_text(aes(label=ifelse(adj.P.Val <= 0.05 & abs(logFC) >= 1.5 & AveExpr > 2,"*",""),size=8.0))+
# 标题和坐标轴设置
labs(x = "Gene", y ="Cancer Type" ) +
  theme(
    plot.title  = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x  = element_text(angle = 90, hjust = 1, size = 7),
    axis.text.y  = element_text(size = 8),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  )
p


ggsave(p, filename="./123-gene-pancaner-heatmap.pdf", width=14,
       height=6)

ggsave(p, filename="./123-gene-pancaner-heatmap.png", width=14,
       height=6)


##20250917
###PANCER 基因表达小提琴图和箱式图
drug1 <- read.delim("e:/01 papers/2025 paper/19 cell cycle-database by Wenyuan Wang/genelist for evolution present.txt",header =T,sep = "\t")
drug1 <- data.frame(drug1)
head(drug1)
gene <- drug1$gene
load(file="D:/TCGA_20220305/TcgaTargetGtex/TcgaTargetGtex_TPM_log2_exp_matrix/BRCA_TPM+0.001_log2_exp_matrix.Rdata")
load("D:/TCGA_20220305/TcgaTargetGtex/TCGA_GTEx_pancancer_mrna_pheno.rdata")
exp_nc[1:6,1:6]
head(colnames(tcga_gtex_mrna_pheno))
table(tcga_gtex_mrna_pheno$sample_type)

table(tcga_gtex_mrna_pheno$project)

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(latex2exp)
library(rstatix)

gene <- gsub(" ","",gene)

for (i in gene) {
  print(i)
  plot_df <- tcga_gtex_mrna_pheno %>%
    select(1:4,all_of(i)) %>%
    filter(sample_type %in% c("GTEx_normal","TCGA_normal","TCGA_tumor"))
  
  colnames(plot_df)[5] <- "gene"
  head(plot_df)
  # 按中位数由高到低排列???
  data_new <- plot_df %>% 
    group_by(project) %>% 
    mutate(median = median(gene), group_max = max(gene)) %>% 
    arrange(desc(median))
  data_new$sample_type <- ifelse(data_new$sample_type != "TCGA_tumor","Normal",data_new$sample_type )
  data_new$sample_type <- ifelse(data_new$sample_type == "TCGA_tumor","Tumor",data_new$sample_type)
  data_new$gene <- ifelse(data_new$gene < -5,5, data_new$gene)
  # 调整因子顺序???
  data_new$project <- factor(data_new$project, levels = unique(data_new$project))
  
  
  data_new <- data.frame(data_new)
  data_new$sample_type <- as.factor(data_new$sample_type)
  
  head(data_new)
  colnames(data_new)[2] <- "Group"
  
  data_new[1:6,]
  
  library(ggpubr) 
  
  data_new<- as.data.frame(data_new)
  
  if (!is.factor(data_new$Group))  { 
    data_new$Group <- as.factor(data_new$Group)  
  } 
  
  p <- ggplot(data = data_new, aes(x = project, y = gene)) + 
    geom_violin(aes(fill = Group)) + 
    geom_boxplot(aes(colour = Group),fill="white",outlier.shape  = NA,width = 0.3, outlier.color = "black",
                 position = position_dodge(0.9)) + 
    
    scale_color_manual(values = rep("black", 2)) + 
    
    theme_bw() + 
    theme(axis.text.x  = element_blank(), 
          legend.position  = 'top', 
          plot.title  = element_text(size = 12, face = "bold.italic"),  
          axis.title.y  = element_text(size = 12, face = "bold")) + 
    ylab(paste0(i, " expression: log2(TPM+0.01)")) + 
    xlab("Cancer") + 
    facet_wrap(. ~ project, scales = "free_x",  nrow = 3) + 
    scale_x_discrete(breaks = unique(data_new$project)) 
  
  p1=p + stat_compare_means(aes(x = project, y = gene, group = Group), 
                            method = "wilcox.test",  
                            label = "p.signif",  
                            size = 5, 
                            hide.ns  = TRUE, 
                            label.x = 1.2,
                            label.y = max(data_new$gene) -1) 
  
  p1
  ggsave(paste0("E:/01 papers/2025 paper/19 cell cycle-database by Wenyuan Wang/out1/",i,".png"),  plot = last_plot(), width = 10, height = 7, dpi = 300)
  ggsave(paste0("E:/01 papers/2025 paper/19 cell cycle-database by Wenyuan Wang/out1/",i,".pdf"),  plot = last_plot(), width = 10, height = 7, dpi = 300)
}



