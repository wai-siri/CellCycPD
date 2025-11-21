

library(dtplyr)
library(magrittr)
library(ggplot2)
#library(Rserve)
library(reshape)
#library(RMySQL)
library(DBI)
#Rserve()


#lowcol="#0B70B7",highcol="#C1282D"


load('/shujupan/ryj/rdata/pancer_exp_more.Rdata')
exp <-exp

cancer <- read.table("/shujupan/ryj/rdata/TCGA.33cancer.sample.number.TPM.txt",header = T,sep = "\t")
cancers <- sort(cancer$cancer)


rownames(exp) <- sub("^.*?\\.", "", rownames(exp))
load("/shujupan/ryj/rdata/immuscore.Rdata")
cancer <- read.table("/shujupan/ryj/rdata/TCGA.33cancer.sample.number.TPM.txt",header = T,sep = "\t")
rownames(immuscore) <- sub("^.*?\\.", "", rownames(immuscore))
cancers <- sort(cancer$cancer)


rownames(exp) <- sub("^.*?\\.", "", rownames(exp))
cancer <- read.table("/shujupan/ryj/rdata/TCGA.33cancer.sample.number.TPM.txt",header = T,sep = "\t")
cancers <- sort(cancer$cancer)
load("/shujupan/ryj/rdata/immucells.Rdata")

gene_checkpoint_heatmap=function(gene,method="pearson",lowcol="#0B70B7",highcol="#C1282D"){
  clist=list()
  plist=list()
  checkpoint=c('CD44','TNFRSF9','LAG3','CD200','NRP1','CD40','CD40LG', 'CD276', 
               'CD86', 'HHLA2', 'CD48', 'CD160', 'TNFSF4', 'CD274', 'TNFSF18', 'TNFRSF8', 
               'CD80', 'CD244', 'TNFSF9', 'CD70', 'TNFSF14', 'ADORA2A', 'IDO1', 'VTCN1', 
               'HAVCR2', 'CD27', 'TNFRSF14', 'ICOSLG', 'CTLA4', 'ICOS', 'CD200R1','LAIR1', 
               'KIR3DL1', 'TMIGD2', 'LGALS9', 'CD28', 'TNFSF15', 'TIGIT', 'BTLA', 'TNFRSF4', 
               'TNFRSF18', 'PDCD1', 'IDO2', 'PDCD1LG2', 'BTNL2', 'TNFRSF25',"SIGLEC15")
  checkpoint=checkpoint[checkpoint %in% colnames(exp)]
  cancers <- sort(cancer$cancer)
  for (cancer in cancers) {
    tg=subset(exp,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(gene))
    
    ic=subset(exp,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(checkpoint))
    identical(rownames(tg),rownames(ic))
    
    # Correlation
    cor <-psych::corr.test(tg, ic, method = method,adjust="none")
    # Cor, p value
    cmt <-cor$r
    pmt <- cor$p
    clist[[cancer]]=as.data.frame(cmt)
    plist[[cancer]]=as.data.frame(pmt)
  }
  clist=do.call(rbind,clist)
  plist=do.call(rbind,plist)
  clist[is.na(clist)]=0
  plist[is.na(plist)]=0.06
  plist=as.data.frame(ifelse(plist<0.01,'**',ifelse(plist<0.05,'*','')))
  clist$cancer=rownames(clist)
  plist$cancer=rownames(plist)
  clist=reshape2::melt(clist,id.vars = "cancer")
  plist=reshape2::melt(plist,id.vars = "cancer")
  
  df=cbind(clist,plist)
  df=df[,c(1:3,6)]
  colnames(df)=c("cancer","checkpoint","cor","sig")
  df$checkpoint=factor(df$checkpoint,levels=rev(unique(df$checkpoint)))
  
  p <- ggplot(df, aes(cancer, checkpoint, fill = cor)) +
    geom_tile(aes(width = 1, height = 1), linewidth = 10) +
    geom_text(aes(cancer, checkpoint, label = sig), color = "black") +
    scale_fill_gradient2(low = lowcol, mid = "white", high = highcol) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank())+
    theme(axis.text.x = element_text(angle = 45,hjust =1,colour = "black",size=10),
          axis.text.y = element_text(colour = "black",size=10),
          axis.title.x = element_text(hjust = 0.5,size = 18 ,face = "bold",colour = "#333333"),
          plot.title = element_text(hjust = 0.5,vjust = 7,size = 20,face = "bold",colour = "#333333"))+
    labs(x='Cancer Type',y='')+
    ggtitle("Immune CheckPoint HeatMap")+
    theme(axis.ticks = element_blank())+coord_equal()
  
  pic_name <- toString(gene)
  pic_name <- paste0("/shujupan/PicStorage/ryj/RMethyMD/data/result/immune",pic_name,"_checkpoint.png",collapse = "")
  ggsave(pic_name,p,width = 11,height = 15)
  
}
gene_checkpoint_heatmap('ZKSCAN3')


gene_chemokine_heatmap=function(gene,method = "pearson",lowcol="#0B70B7",highcol="#C1282D"){
  clist=list()
  plist=list()
  chemokine=c("CCL1","CCL2","CCL3","CCL4","CCL5","CCL7","CCL8","CCL11","CCL13","CCL14","CCL15","CCL16","CCL17","CCL18","CCL19","CCL20","CCL21","CCL22","CCL23","CCL24","CCL25","CCL26","CCL28","CX3CL1","CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8","CXCL9","CXCL10","CXCL11","CXCL12","CXCL13","CXCL14","CXCL16","CXCL17")
  chemokine=chemokine[chemokine %in% colnames(exp)]
  cancers <- sort(cancer$cancer)
  for (cancer in cancers) {
    tg=subset(exp,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(gene))
    
    ic=subset(exp,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(chemokine))
    identical(rownames(tg),rownames(ic))
    
    # Correlation
    cor <-psych::corr.test(tg, ic, method = method,adjust="none")
    # Cor, p value
    cmt <-cor$r
    pmt <- cor$p
    clist[[cancer]]=as.data.frame(cmt)
    plist[[cancer]]=as.data.frame(pmt)
  }
  clist=do.call(rbind,clist)
  plist=do.call(rbind,plist)
  clist[is.na(clist)]=0
  plist[is.na(plist)]=0.06
  plist=as.data.frame(ifelse(plist<0.01,'**',ifelse(plist<0.05,'*','')))
  
  clist$cancer=rownames(clist)
  plist$cancer=rownames(plist)
  clist=melt(clist,id.vars = "cancer")
  plist=melt(plist,id.vars = "cancer")
  
  df=cbind(clist,plist)
  df=df[,c(1:3,6)]
  colnames(df)=c("cancer","chemokine","cor","sig")
  df$chemokine=factor(df$chemokine,levels=rev(unique(df$chemokine)))
  
  ggplot(df, aes(cancer, chemokine, fill = cor)) +
    geom_tile(aes(width = 1, height = 1), size = 10) +
    geom_text(aes(cancer, chemokine, label = sig), color = "black") +
    scale_fill_gradient2(low = lowcol, mid = "white", high = highcol) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank())+
    theme(axis.text.x = element_text(angle = 45,hjust =1,colour = "black",size=8),
          axis.text.y = element_text(colour = "black",size=8),
          axis.title.x = element_text(hjust = 0.5,size = 10 ,face = "bold",colour = "#333333"),
          plot.title = element_text(hjust = 0.5,vjust = 1,size = 12,face = "bold",colour = "#333333"))+
    labs(x='Cancer Type',y='')+
    ggtitle("Immune Chemokine HeatMap")+
    theme(axis.ticks = element_blank())
  pic_name <- toString(gene)
  pic_name <- paste0("/shujupan/PicStorage/ryj/RMethyMD/data/result/immune",pic_name,"_chemokine.png",collapse = "")
  ggsave(pic_name,width = 11,height = 6,dpi = 300)
  
  #df_chemokine <- cbind(gene,df)
  #write.csv(df_chemokine, file = "D:/data/immune/chemokine.csv", fileEncoding = "UTF-8")
  #print("save chemokine success")
}
gene_chemokine_heatmap('METTL3')

gene_receptor_heatmap=function(gene,method = "pearson",lowcol="#0B70B7",highcol="#C1282D"){
  clist=list()
  plist=list()
  receptor=c("CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CCR8","CCR9","CCR10",
             "CXCR1","CXCR2","CXCR3","CXCR4","CXCR5","CXCR6","XCR1","CX3R1")
  receptor=receptor[receptor %in% colnames(exp)]
  cancers <- sort(cancer$cancer)
  for (cancer in cancers) {
    tg=subset(exp,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(gene))
    
    ic=subset(exp,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(receptor))
    identical(rownames(tg),rownames(ic))
    
    # Correlation
    cor <-psych::corr.test(tg, ic, method = method,adjust="none")
    # Cor, p value
    cmt <-cor$r
    pmt <- cor$p
    clist[[cancer]]=as.data.frame(cmt)
    plist[[cancer]]=as.data.frame(pmt)
  }
  clist=do.call(rbind,clist)
  plist=do.call(rbind,plist)
  clist[is.na(clist)]=0
  plist[is.na(plist)]=0.06
  plist=as.data.frame(ifelse(plist<0.01,'**',ifelse(plist<0.05,'*','')))
  
  clist$cancer=rownames(clist)
  plist$cancer=rownames(plist)
  clist=reshape2::melt(clist,id.vars = "cancer")
  plist=reshape2::melt(plist,id.vars = "cancer")
  
  df=cbind(clist,plist)
  df=df[,c(1:3,6)]
  colnames(df)=c("cancer","receptor","cor","sig")
  df$receptor=factor(df$receptor,levels=rev(unique(df$receptor)))
  
  ggplot(df, aes(cancer, receptor, fill = cor)) +
    geom_tile(aes(width = 1, height = 1), size = 10) +
    geom_text(aes(cancer, receptor, label = sig), color = "black") +
    scale_fill_gradient2(low = lowcol, mid = "white", high = highcol) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank())+
    theme(axis.text.x = element_text(angle = 45,hjust =1,colour = "black",size=12),
          axis.text.y = element_text(colour = "black",size=10),
          axis.title.x = element_text(hjust = 0.5,size = 15 ,face = "bold",colour = "#333333"),
          plot.title = element_text(hjust = 0.5,vjust = 1,size = 16,face = "bold",colour = "#333333"))+
    labs(x='Cancer Type',y='')+
    ggtitle("Immune Receptor HeatMap")+
    theme(axis.ticks = element_blank())
  
  pic_name <- toString(gene)
  pic_name <- paste0("/shujupan/PicStorage/ryj/RMethyMD/data/result/immune",pic_name,"_receptor.png",collapse = "")
  ggsave(pic_name,width = 11,height = 6,dpi = 300)
  
  #df_receptor <- cbind(gene,df)
  #write.csv(df_receptor, file = "D:/data/immune/receptor.csv", fileEncoding = "UTF-8")
  #print("save receptor success")
}
#gene_receptor_heatmap("METTL3")



gene_immustimulator_heatmap=function(gene,method = "pearson",lowcol="#0B70B7",highcol="#C1282D"){
  clist=list()
  plist=list()
  immustimulator=c("CD27","CD276","CD28","CD40","CD40LG","CD48","CD70","CD80","CD86","CXCL12","CXCR4","ENTPD1","HHLA2","ICOS","ICOSLG","IL2RA","IL6","IL6R","KLRC1","KLRK1","LTA","MICB","NT5E","PVR","RAET1E","TMIGD2","TNFRSF13B","TNFRSF13C","TNFRSF14","TNFRSF17","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF13","TNFSF13B","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF9","ULBP1")
  immustimulator=immustimulator[immustimulator %in% colnames(exp)]
  cancers <- sort(cancer$cancer)
  for (cancer in cancers) {
    tg=subset(exp,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(gene))
    
    ic=subset(exp,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(immustimulator))
    identical(rownames(tg),rownames(ic))
    
    # Correlation
    cor <-psych::corr.test(tg, ic, method = method,adjust="none")
    # Cor, p value
    cmt <-cor$r
    pmt <- cor$p
    clist[[cancer]]=as.data.frame(cmt)
    plist[[cancer]]=as.data.frame(pmt)
  }
  clist=do.call(rbind,clist)
  plist=do.call(rbind,plist)
  clist[is.na(clist)]=0
  plist[is.na(plist)]=0.06
  plist=as.data.frame(ifelse(plist<0.01,'**',ifelse(plist<0.05,'*','')))
  
  clist$cancer=rownames(clist)
  plist$cancer=rownames(plist)
  clist=reshape2::melt(clist,id.vars = "cancer")
  plist=reshape2::melt(plist,id.vars = "cancer")
  
  df=cbind(clist,plist)
  df=df[,c(1:3,6)]
  colnames(df)=c("cancer","immustimulator","cor","sig")
  df$immustimulator=factor(df$immustimulator,levels=rev(unique(df$immustimulator)))
  
  ggplot(df, aes(cancer, immustimulator, fill = cor)) +
    geom_tile(aes(width = 1, height = 1), size = 10) +
    geom_text(aes(cancer, immustimulator, label = sig), color = "black") +
    scale_fill_gradient2(low = lowcol, mid = "white", high = highcol) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank())+
    theme(axis.text.x = element_text(angle = 45,hjust =1,colour = "black",size=6),
          axis.text.y = element_text(colour = "black",size=7),
          axis.title.x = element_text(hjust = 0.5,size = 8 ,face = "bold",colour = "#333333"),
          plot.title = element_text(hjust = 0.5,vjust = 1,size = 9,face = "bold",colour = "#333333"))+
    labs(x='Cancer Type',y='')+
    ggtitle("Immustimulator HeatMap")+
    theme(axis.ticks = element_blank())
  pic_name <- toString(gene)
  pic_name <- paste0("/shujupan/PicStorage/ryj/RMethyMD/data/result/immune",pic_name,"_immustimulator.png",collapse = "")
  ggsave(pic_name,dpi = 300)
  
  #df_immustimulator <- cbind(gene,df)
  #write.csv(df_immustimulator, file = "D:/data/immune/immustimulator.csv", fileEncoding = "UTF-8")
  #print("save immustimulator success")
}
#gene_immustimulator_heatmap("METTL3")




gene_immuinhibitor_heatmap=function(gene,method = "pearson",lowcol="#0B70B7",highcol="#C1282D"){
  clist=list()
  plist=list()
  immuinhibitor=c("ADORA2A","BTLA","CD160","CD244","CD274","CD96","CSF1R","CTLA4","HAVCR2","IDO1","IL10","IL10RB","KDR","KIR2DL1","KIR2DL3","LAG3","LGALS9","PDCD1","PDCD1LG2","TGFB1","TGFBR1","TIGIT","VTCN1")
  immuinhibitor=immuinhibitor[immuinhibitor %in% colnames(exp)]
  cancers <- sort(cancer$cancer)
  for (cancer in cancers) {
    tg=subset(exp,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(gene))
    
    ic=subset(exp,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(immuinhibitor))
    identical(rownames(tg),rownames(ic))
    
    # Correlation
    cor <-psych::corr.test(tg, ic, method = method,adjust="none")
    # Cor, p value
    cmt <-cor$r
    pmt <- cor$p
    clist[[cancer]]=as.data.frame(cmt)
    plist[[cancer]]=as.data.frame(pmt)
  }
  clist=do.call(rbind,clist)
  plist=do.call(rbind,plist)
  clist[is.na(clist)]=0
  plist[is.na(plist)]=0.06
  plist=as.data.frame(ifelse(plist<0.01,'**',ifelse(plist<0.05,'*','')))
  
  clist$cancer=rownames(clist)
  plist$cancer=rownames(plist)
  clist=reshape2::melt(clist,id.vars = "cancer")
  plist=reshape2::melt(plist,id.vars = "cancer")
  
  df=cbind(clist,plist)
  df=df[,c(1:3,6)]
  colnames(df)=c("cancer","immuinhibitor","cor","sig")
  df$immuinhibitor=factor(df$immuinhibitor,levels=rev(unique(df$immuinhibitor)))
  
  ggplot(df, aes(cancer, immuinhibitor, fill = cor)) +
    geom_tile(aes(width = 1, height = 1), size = 10) +
    geom_text(aes(cancer, immuinhibitor, label = sig), color = "black") +
    scale_fill_gradient2(low = lowcol, mid = "white", high = highcol) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank())+
    theme(axis.text.x = element_text(angle = 45,hjust =1,colour = "black",size=6),
          axis.text.y = element_text(colour = "black",size=8),
          axis.title.x = element_text(hjust = 0.5,size = 10 ,face = "bold",colour = "#333333"),
          plot.title = element_text(hjust = 0.5,vjust = 1,size = 10,face = "bold",colour = "#333333"))+
    labs(x='Cancer Type',y='')+
    ggtitle("Immuinhibitor HeatMap")+
    theme(axis.ticks = element_blank())
  pic_name <- toString(gene)
  pic_name <- paste0("/shujupan/PicStorage/ryj/RMethyMD/data/result/immune",pic_name,"_immuinhibitor.png",collapse = "")
  ggsave(pic_name,dpi = 300)
  
  #df_immuinhibitor<- cbind(gene,df)
  #write.csv(df_immuinhibitor, file = "D:/data/immune/immuinhibitor.csv", fileEncoding = "UTF-8")
  #print("save immuinhibitor success")
}
#gene_immuinhibitor_heatmap("METTL3")


gene_immunescore_heatmap=function(gene,method="pearson",lowcol="#0B70B7",highcol="#C1282D"){
  clist=list()
  plist=list()
  for (cancer in cancers) {
    tg=subset(exp,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(gene))%>%
      dplyr::filter(rownames(.) %in% rownames(immuscore))
    
    ic=immuscore[rownames(tg),]
    identical(rownames(tg),rownames(ic))
    
    # Correlation
    cor <-psych::corr.test(tg, ic, method = method,adjust="none")
    # Cor, p value
    cmt <-cor$r
    pmt <- cor$p
    clist[[cancer]]=as.data.frame(cmt)
    plist[[cancer]]=as.data.frame(pmt)
    
  }
  
  clist=do.call(rbind,clist)
  plist=do.call(rbind,plist)
  clist[is.na(clist)]=0
  plist[is.na(plist)]=0.06
  plist=as.data.frame(ifelse(plist<0.01,'**',ifelse(plist<0.05,'*','')))
  
  clist$cancer=rownames(clist)
  plist$cancer=rownames(plist)
  clist=reshape2::melt(clist,id.vars = "cancer")
  plist=reshape2::melt(plist,id.vars = "cancer")
  
  df=cbind(clist,plist)
  df=df[,c(1:3,6)]
  colnames(df)=c("cancer","score","cor","sig")
  df$score=factor(df$score,levels=rev(unique(df$score)))
  
  ggplot(df, aes(cancer, score, fill = cor)) +
    geom_tile(aes(width = 2, height = 2), linewidth = 10) +
    geom_text(aes(cancer, score, label = sig), color = "black") +
    scale_fill_gradient2(low = lowcol, mid = "white", high = highcol) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank())+
    theme(axis.text.x = element_text(angle = 45,hjust =1,colour = "black",size=6),
          axis.text.y = element_text(colour = "black",size=8),
          axis.title.x = element_text(hjust = 0.5,size = 8 ,face = "bold",colour = "#333333"),
          plot.title = element_text(hjust = 0.5,vjust = 1,size = 8,face = "bold",colour = "#333333"))+
    labs(x='Cancer Type',y='')+
    ggtitle("ImmuneScore HeatMap")+
    theme(axis.ticks = element_blank())+coord_equal()
  
  pic_name <- toString(gene)
  pic_name <- paste0("/shujupan/PicStorage/ryj/RMethyMD/data/result/immune",pic_name,"_immunescore.png",collapse = "")
  ggsave(pic_name,dpi = 300)
  
  #df_score <-cbind(gene,df)
  #write.csv(df_score,file = "D:/data/immune/score.csv", fileEncoding = "UTF-8")
  #print("save score success")
}
#gene_immunescore_heatmap('METTL3')

gene_immucell_heatmap=function(gene,method="pearson",lowcol="#0B70B7",highcol="#C1282D"){
  clist=list()
  plist=list()
  for (cancer in cancers) {
    tg=subset(exp,Group=="Tumor" & Cancer==cancer)%>%
      dplyr::select(all_of(gene))%>%
      dplyr::filter(rownames(.) %in% rownames(immucells))
    
    ic=immucells[rownames(tg),]
    identical(rownames(tg),rownames(ic))
    
    # Correlation
    cor <-psych::corr.test(tg, ic, method = method,adjust="none")
    # Cor, p value
    cmt <-cor$r
    pmt <- cor$p
    clist[[cancer]]=as.data.frame(cmt)
    plist[[cancer]]=as.data.frame(pmt)
  }
  clist=do.call(rbind,clist)
  plist=do.call(rbind,plist)
  clist[is.na(clist)]=0
  plist[is.na(plist)]=0.06
  plist=as.data.frame(ifelse(plist<0.01,'**',ifelse(plist<0.05,'*','')))
  
  clist$cancer=rownames(clist)
  plist$cancer=rownames(plist)
  clist=reshape2::melt(clist,id.vars = "cancer")
  plist=reshape2::melt(plist,id.vars = "cancer")
  
  df=cbind(clist,plist)
  df=df[,c(1:3,6)]
  colnames(df)=c("cancer","cell","cor","sig")
  df$cell=factor(df$cell,levels=rev(unique(df$cell)))
  
  p <- ggplot(df, aes(cancer, cell, fill = cor)) +
    geom_tile(aes(width = 1, height =1), linewidth = 10) +
    geom_text(aes(cancer, cell, label = sig), color = "black") +
    scale_fill_gradient2(low = lowcol, mid = "white", high = highcol) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank())+
    theme(axis.text.x = element_text(angle = 45,hjust =1,colour = "black",size=14),
          axis.text.y = element_text(colour = "black",size=16),
          axis.title.x = element_text(hjust = 0.5,size = 18 ,face = "bold",colour = "#333333"),
          plot.title = element_text(hjust = 0.5,vjust = 5,size = 18,face = "bold",colour = "#333333"))+
    labs(x='Cancer Type',y='')+
    ggtitle("ImmuCell HeatMap")+
    theme(axis.ticks = element_blank())+coord_equal()
  
  pic_name <- toString(gene)
  pic_name <- paste0("/shujupan/PicStorage/ryj/RMethyMD/data/result/immune",pic_name,"_immucell.png",collapse = "")
  ggsave(pic_name,p,height = 10,width = 15,dpi = 300)
  
  #df_cell <-cbind(gene,df)
  #write.csv(df_cell,file = "D:/data/immune/cell.csv", fileEncoding = "UTF-8")
  #print("save cell success")
}
gene_immucell_heatmap('METTL3')



