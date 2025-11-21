# 载入必要的数据和库  
library(RIdeogram)
library(Cairo)  # 加载 Cairo 包

data(human_karyotype, package="RIdeogram")  
data(gene_density, package="RIdeogram")  

# 使用绝对路径读取文件  
label <- read.table("D:/R_Project/cellcycle/code/tool/gene_on_chr/code_use_position.txt",   
                    sep = "\t", header = TRUE, stringsAsFactors = FALSE)  

# 定义调整后的配色方案  
colorset1 <- c(  
  "#f0f0f0",  
  "#bdbdbd", 
  "#67a9cf",  
  "#ef8a62", 
  "#d73027"   # 深红色，用于非常高密度  
)  

CairoSVG(file = 'chromosome_6.svg', width = 11, height = 12)  # 使用 CairoSVG
ideogram(karyotype = human_karyotype, overlaid = gene_density, label, label_type = "marker",   
         colorset1 = colorset1, output='chromosome_6.svg')  

# 将 SVG 图片转换为 PNG  
convertSVG("chromosome_6.svg", device = "png", file = 'chromosome_6.png', dpi = 400)