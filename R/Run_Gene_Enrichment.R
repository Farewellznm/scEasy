#' Title  Gene enrichment analysis ORA
#' @description
#' This function is used to human or mouse gene ORA analysis
#'
#' @param geneList geneList is which gene list you want to ORA analysis
#' @param org The spesis names. defalt = "hsa"
#'
#' @return null
#' @export
#'
#' @examples  Run_Gene_Enrichment(geneList = genelist,org = "hsa")
Run_Gene_Enrichment <- function(geneList = geneList,org = "hsa") {
  #geneList = gl
  #org = "hsa"

  if (org == "hsa") {
    anno = "org.Hs.eg.db"
  }
  if (org == "mmu") {
    anno = "org.Mm.eg.db"
  }
  if (org == "ath") {
    anno = "org.At.tair.db"
  }
  if (org == "rno") {
    anno = "org.Rn.eg.db"
  }

  if (org == "ssc") {
    anno = "org.Ss.eg.db"
  }
  print(paste("The annotation database:",anno))


  # organism = mouse
  # 小鼠 使用的注释包 ： 	org.Mm.eg.db
  # annotationDB = "org.Mm.eg.db"
  # org = "mmu"

  # organism = Rat
  # 大鼠 使用的注释包 ： 	org.Rn.eg.db
  # annotationDB = "org.Rn.eg.db"
  # org = "rno"

  # organism = pig
  # 猪 使用的注释包 ： org.Ss.eg.db
  # annotationDB = "org.Ss.eg.db"
  # org = "ssc"

  # 可选择模式物种
  # 按蚊（Anopheles）	    org.Ag.eg.db         organm = aga
  # 拟南芥（Arabidopsis）	org.At.tair.db       organm = ath **
  # 牛（Brovine)	        org.Bt.eg.db         organm = bta
  # 犬（Canine）	        org.Cf.eg.db         organm = cfa
  # 黑腹果蝇（Drosophila melanogaster) org.Dm.eg.db         organm = dme
  # 斑马鱼（Zebrafish） 	org.Dr.eg.db         organm = dre
  # 鸡（Chicken）	        org.Gg.eg.db         organm = gga
  # 人（Humanm)	          org.Hs.eg.db         organm = hsa
  # 小鼠（Mouse）	        org.Mm.eg.db         organm = mmu
  # 黑猩猩（Chimp）	      org.Pt.eg.db         organm = ptr
  # 大鼠，褐家鼠（Rat)	  org.Rn.eg.db         organm = rno **
  # 猪（Pig)	            org.Ss.eg.db         organm = ssc **

  # 对于非模式物种 需要手动构建orgDB注释包

  # ------------------------ 环境准备 ----------------------------
  ##################################
  # 选择模式物种需要安装对应的注释包
  ##################################
  suppressPackageStartupMessages({
    if(org == "ssc") {
      if (!require("org.Ss.eg.db")) {BiocManager::install("org.Ss.eg.db")}
    }

    if(org == "rno") {
      if (!require("org.Rn.eg.db")) {BiocManager::install("org.Rn.eg.db")}
    }

    if(org == "ath") {
      if (!require("org.At.tair.db")) {BiocManager::install("org.At.tair.db")}
    }

    if(org == "mmu") {
      if(!require("org.Mm.eg.db")) { BiocManager::install("org.Mm.eg.db")}
    }

    if(org == "hsa") {
      if(!require("org.Hs.eg.db")) { BiocManager::install("org.Hs.eg.db")}
    }

    if(!require(BiocManager)) {install.packages("BiocManager")}
    #---------------------------------
    if(!require(stringr)){install.packages("stringr")}
    if(!require(dplyr)){install.packages("dplyr")}
    if(!require(ggplot2)){install.packages("ggplot2")}
    if(!require(clusterProfiler)) { BiocManager::install("clusterProfiler")}
    if(!require(enrichplot)) { BiocManager::install("enrichplot")
      require(enrichplot)}
    if(!require(cowplot)) { BiocManager::install("cowplot")}
    if(!require("ggtree")){devtools::install_github("YuLab-SMU/ggtree")}
    if(!require("meme")){install.packages("meme")}
    if(!require("GOplot")){install.packages("GOplot")}
    if(!require("yyplot")){devtools::install_github("GuangchuangYu/yyplot")}
    if (!require(R.utils)) {
      install.packages("R.utils")
    }
    R.utils::setOption("clusterProfiler.download.method", 'auto')
  })

  gene = bitr(geneList,
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = anno)

  ################ GO BP ##############################
  go <-
    enrichGO(
      gene          = gene$ENTREZID,
      #universe = all_entrezid$ENTREZID,
      OrgDb         = anno,
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 1,
      readable      = TRUE
    )
  # 去冗余
  go <- clusterProfiler::simplify(go, cutoff = 0.7, by = "p.adjust", select_fun = min)
  go_res = go@result


  barplot(go, showCategory = 10,color = "p.adjust")
  ggsave("GO_pathways_barplot.pdf", width = 8,height = 6)

  dotplot(go, showCategory = 10)
  ggsave("GO_pathways_dotplot.pdf",width = 8,height = 6)

  write.csv(go_res,"go_res.csv")
  #### goplot  zscore (up - down) / sqrt(count-term)
  # 绘制GO富集分析条形图，结果默认按qvalue升序，分别选出前十的term进行绘图即可

  # 如果没有10个会报错
  goBP <- subset(go_res,subset = (ONTOLOGY == "BP"))[1:10,]
  goCC <- subset(go_res,subset = (ONTOLOGY == "CC"))[1:10,]
  goMF <- subset(go_res,subset = (ONTOLOGY == "MF"))[1:10,]
  go.df <- rbind(goBP,goCC,goMF)

  go.df$ONTOLOGY <- factor(go.df$ONTOLOGY,levels = rev(c("BP","CC","MF")))
  go.df <- go.df[order(go.df$ONTOLOGY,go.df$Count,decreasing = T),]
  # 使画出的GO term的顺序与输入一致
  go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
  # 绘图
  go_bar <- ggplot(data = go.df, # 绘图使用的数据
                   aes(x = Description, y = Count,fill = ONTOLOGY))+ # 横轴坐标及颜色分类填充
    geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
    coord_flip()+theme_classic()+ # 横纵坐标反转及去除背景色
    scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
    labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ # 设置坐标轴标题及标题
    scale_fill_manual(values = rev(c("#26547C","#EF476F","#FFD166"))) +
    #scale_fill_manual(values  = c("#A13425","#E0C9A7","#283870")) +  # 配色方案 mumuxi推荐
    theme(axis.title = element_text(size = 13), # 坐标轴标题大小
          axis.text = element_text(size = 11), # 坐标轴标签大小
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
          legend.title = element_text(size = 13), # 图例标题大小
          legend.text = element_text(size = 11)) # 图边距
  go_bar
  ggsave(go_bar,filename = "GO_Barplot.pdf",width = 9,height = 7)
  ggsave(go_bar,filename = "GO_Barplot.png",width = 9,height = 7)
  ###################### KEGG ######################
  kk = enrichKEGG(gene = gene$ENTREZID,organism = org,pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,qvalueCutoff = 1)
  library(DOSE)
  #如果原始的ID号为entrez gene id那么这里keyType设置为ENTREZID
  kk <-setReadable(kk, OrgDb = anno , keyType="ENTREZID")
  reskk = kk@result
  reskk = reskk[reskk$p.adjust < 0.05, ]
  write.csv(reskk, file = "KEGG_result.csv", row.names = F)
  barplot(kk, showCategory = 10)
  ggsave("KEGG_Pathways_barplot.pdf",width = 8,height = 6)
  ggsave("KEGG_Pathways_barplot.png",width = 8,height = 6)
  options(digits  = 2)
  dotplot(kk, showCategory = 10)
  ggsave("KEGG_Pathways_dotplot.pdf",width = 8,height = 6)
  ggsave("KEGG_Pathways_dotplot.png",width = 8,height = 6)

  print("The jobs was finished!")
}

################################################################################
re_plot <- function(gofile,keggfile) {

  require(ggplot2)
  require(tidyverse)
  #require(enrichplot)
  gofile <- readxl::read_excel(gofile)

  keggfile <- readxl::read_excel(keggfile)


  gofile <- as.data.frame(gofile)
  keggfile <- as.data.frame(keggfile)

  ##############################################################################
  go_res <- gofile
  # 绘制GO富集分析条形图，结果默认按qvalue升序，分别选出前十的term进行绘图即可
  goBP <- subset(go_res,subset = (ONTOLOGY == "BP"))
  goCC <- subset(go_res,subset = (ONTOLOGY == "CC"))
  goMF <- subset(go_res,subset = (ONTOLOGY == "MF"))
  go.df <- rbind(goBP,goCC,goMF)

  go.df$ONTOLOGY <- factor(go.df$ONTOLOGY,levels = rev(c("BP","CC","MF")))
  go.df <- go.df[order(go.df$ONTOLOGY,go.df$Count,decreasing = T),]

  # 使画出的GO term的顺序与输入一致
  go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))

  # 绘图
  go_bar <- ggplot(data = go.df, # 绘图使用的数据
                   aes(x = Description, y = Count,fill = ONTOLOGY))+ # 横轴坐标及颜色分类填充
    geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
    coord_flip()+theme_classic()+ # 横纵坐标反转及去除背景色
    scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
    labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ # 设置坐标轴标题及标题
    scale_fill_manual(values = rev(c("#26547C","#EF476F","#FFD166"))) +
    #scale_fill_manual(values  = c("#A13425","#E0C9A7","#283870")) +  # 配色方案 mumuxi推荐
    theme(axis.title = element_text(size = 13), # 坐标轴标题大小
          axis.text = element_text(size = 10), # 坐标轴标签大小
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
          legend.title = element_text(size = 13), # 图例标题大小
          legend.text = element_text(size = 11))# 图边距
  print(go_bar)
  ggsave(go_bar,filename = "GO_Barplot_replot.pdf",width = 9,height = 7)
  ggsave(go_bar,filename = "GO_Barplot_replot.png",width = 9,height = 7)


  kegg_res <- keggfile


  kegg_res <- kegg_res[,colnames(kegg_res) != "ID"]
  #kegg_res <- kegg_res[order(kegg_res$qvalue,kegg_res$GeneRatio),]

  library(plyr)
  library(stringr)
  library(grid)
  library(ggplot2)
  e_data <- kegg_res
  e_data <- e_data[,-1]
  # 分数转小数
  mixedToFloat <- function(x){
    x <- sapply(x, as.character)
    is.integer  <- grepl("^-?\\d+$", x)
    is.fraction <- grepl("^-?\\d+\\/\\d+$", x)
    is.float <- grepl("^-?\\d+\\.\\d+$", x)
    is.mixed    <- grepl("^-?\\d+ \\d+\\/\\d+$", x)
    stopifnot(all(is.integer | is.fraction | is.float | is.mixed))

    numbers <- strsplit(x, "[ /]")

    ifelse(is.integer,  as.numeric(sapply(numbers, `[`, 1)),
           ifelse(is.float,    as.numeric(sapply(numbers, `[`, 1)),
                  ifelse(is.fraction, as.numeric(sapply(numbers, `[`, 1)) /
                           as.numeric(sapply(numbers, `[`, 2)),
                         as.numeric(sapply(numbers, `[`, 1)) +
                           as.numeric(sapply(numbers, `[`, 2)) /
                           as.numeric(sapply(numbers, `[`, 3)))))

  }

  e_data_1 <- e_data
  e_data_1$GeneRatio = mixedToFloat(e_data_1$GeneRatio)
  e_data_1$Count = mixedToFloat(e_data_1$Count)
  log_name <- "-Log10(qvalue)"
  col_name_e_1 <- colnames(e_data_1)
  col_name_e_1 <- c(col_name_e_1,log_name)
  e_data_1$log_name <- log10(e_data_1$qvalue) * (-1)
  colnames(e_data_1) <- col_name_e_1

  e_data_1_freq <- as.data.frame(table(e_data_1$Description))
  colnames(e_data_1_freq) <- c("Description","ID")
  head(e_data_1_freq)


  e_data_2 <- merge(e_data_1,e_data_1_freq,by="Description")
  e_data_3 <- e_data_2[order(e_data_2$ID,
                             e_data_2$GeneRatio,
                             e_data_2$`-Log10(qvalue)`),]


  t_order <- unique(e_data_3$Description)
  e_data_1$Description <- factor(e_data_1$Description,
                                 levels = t_order,ordered = T)

  color_1 <- c("green","red")
  p <- ggplot(e_data_1,aes(x=GeneRatio,y=Description)) +
    labs(x="GeneRatio",y="GO description") + labs(title="")

  p
  p <- p + geom_point(aes(size=Count,color = `-Log10(qvalue)`)) +
    scale_color_gradient(low = color_1[1],high=color_1[2],name="-Log10(qvalue)")
  p
  p <- p + scale_y_discrete(labels=function(x) str_wrap(x,width = 60))
  p
  ggsave(filename = "kegg_dotplot_replot.pdf",width = 9,height = 7)
  ggsave(filename = "kegg_dotplot_replot.png",width = 9,height = 7)
  message("The replot was finished!")
}

