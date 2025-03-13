# 设置工作目录为当前脚本放置目录
# 差异基因文件 xxx.txt放置于当前目录的data目录下


Run_Gene_Enrichment <- function(geneList = geneList,org = "hsa") {

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
    
    if(!require(ggh4x)){install.packages("ggh4x")}
    if(!require(ggfun)){install.packages("ggfun")}
    if(!require(ggnewscale)){install.packages("ggnewscale")}
    if(!require(grid)){install.packages("grid")}
    
    
    if(!require(stringr)){install.packages("stringr")}
    if(!require(DOSE)){BiocManager::install("DOSE")}
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

  
  gene = suppressWarnings(bitr(geneList,
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = anno))

  ################ GO BP ##############################
  gene = na.omit(gene)
  go <- enrichGO(gene = gene$ENTREZID,
      OrgDb         = anno,
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
  # 去冗余,如果不调用clusterProfilter中的simplify,会报错
  go <- clusterProfiler::simplify(go, cutoff = 0.7, by = "p.adjust", select_fun = min)
  go_res <- setReadable(go, OrgDb = anno, keyType = "ENTREZID") %>% as.data.frame()
  

  barplot(go, showCategory = 10,color = "p.adjust")
  ggsave("GO_pathways_barplot.pdf", width = 8,height = 6)
  ggsave("GO_pathways_barplot.png", width = 8,height = 6)

  clusterProfiler::dotplot(go, showCategory = 10)
  ggsave("GO_pathways_dotplot.pdf",width = 8,height = 6)
  ggsave("GO_pathways_dotplot.png",width = 8,height = 6)

  write.csv(go_res,"go_res.csv")
  
  #### goplot  zscore (up - down) / sqrt(count-term)
  # 绘制GO富集分析条形图，结果默认按qvalue升序，分别选出前十的term进行绘图即可
  
  go_res = go_res[order(go_res$p.adjust),]
  # 如果没有10个会报错
  goBP <- subset(go_res,subset = (ONTOLOGY == "BP"))[1:10,]
  goCC <- subset(go_res,subset = (ONTOLOGY == "CC"))[1:10,]
  goMF <- subset(go_res,subset = (ONTOLOGY == "MF"))[1:10,]
  go.df <- rbind(goBP,goCC,goMF)
  go.df <- na.omit(go.df)

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
  #如果原始的ID号为entrez gene id那么这里keyType设置为ENTREZID
  kk <-setReadable(kk, OrgDb = anno , keyType="ENTREZID")
  reskk = kk@result
  
  reskk = reskk[reskk$p.adjust < 0.05, ]
  write.csv(reskk, file = "KEGG_result.csv", row.names = F)
  barplot(kk, showCategory = 10)
  ggsave("KEGG_Pathways_barplot.pdf",width = 6,height = 6)
  ggsave("KEGG_Pathways_barplot.png",width = 6,height = 6)
  options(digits  = 2)
  dotplot(kk, showCategory = 10)
  ggsave("KEGG_Pathways_dotplot.pdf",width = 6,height = 6)
  ggsave("KEGG_Pathways_dotplot.png",width = 6,height = 6)
  
  if (org == "mmu") {
    reskk1 = kk@result
    reskk1 = na.omit(reskk1)
    reskk1$Description = str_split(reskk1$Description,pattern = " - Mus musculus ",simplify = T)[,1]
    reskk1 = reskk1[order(reskk1$p.adjust),]
    write.csv(reskk1, file = "KEGG_result.csv", row.names = F)
    
    kk@result = reskk1
    barplot(kk, showCategory = 10)
    ggsave("KEGG_Pathways_barplot.pdf",width = 8,height = 6)
    ggsave("KEGG_Pathways_barplot.png",width = 8,height = 6)
    options(digits  = 2)
    dotplot(kk, showCategory = 10)
    ggsave("KEGG_Pathways_dotplot.pdf",width = 8,height = 6)
    ggsave("KEGG_Pathways_dotplot.png",width = 8,height = 6)
    reskk = reskk1
  }
  
  
  # 另外一个风格的
  suppressWarnings(replot_2(deg_go_2 = go_res,deg_kegg_2 = reskk))
  
  suppressWarnings(GO_KEGG_plot_bar(go = go_res,kegg = reskk,color = NULL))
  
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
  ggsave(filename = "kegg_dotplot_replot.pdf",width = 7,height = 7)
  ggsave(filename = "kegg_dotplot_replot.png",width = 7,height = 7)
  print("The replot was finished!")
}

replot_2 <- function(deg_go_2 = go_res,deg_kegg_2 = reskk){
  # 挑选top10
  GO_top10 <- deg_go_2 %>%
    dplyr::group_by(ONTOLOGY) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice(1:5) %>%
    dplyr::ungroup()
  
  KEGG_top10 <- deg_kegg_2 %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice(1:5) %>%
    dplyr::select(ID:Count) %>%
    dplyr::mutate(ONTOLOGY = "KEGG") %>%
    dplyr::select(ONTOLOGY, everything())
  
  plot_df <- rbind(GO_top10, KEGG_top10) %>%
    dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")), ordered = T)) %>%
    dplyr::arrange(ONTOLOGY, desc(Count)) %>%
    dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>%
    dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered = T))
  
  ####----Plot----####
  plot <- plot_df %>%
    ggplot() + 
    geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
               aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) + 
    scale_fill_gradient(low = "#a1d99b", high = "#238b45", name = "KEGG p.adjust") + 
    ggnewscale::new_scale_fill() + 
    geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "MF"),
               aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) + 
    scale_fill_gradient(low = "#a6bddb", high = "#0570b0", name = "MF p.adjust") + 
    ggnewscale::new_scale_fill() + 
    geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "CC"),
               aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) + 
    scale_fill_gradient(low = "#fdd49e", high = "#d7301f", name = "CC p.adjust") + 
    ggnewscale::new_scale_fill() +
    geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
               aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) + 
    scale_fill_gradient(low = "#8c96c6", high = "#8c6bb1",  name = "BP p.adjust") + 
    guides(y = "axis_nested",
           y.sec = guide_axis_manual(breaks = 1:nrow(plot_df),
                                     labels = plot_df$Description)) + 
    ggtitle(label = "GO and KEGG annotation") + 
    labs(x = "Count", y = "Description") + 
    scale_size(range = c(3, 7),
               guide = guide_legend(override.aes = list(fill = "#000000"))) + 
    theme_bw() + 
    theme(
      ggh4x.axis.nestline.y = element_line(size = 3, color = c("#74c476", "#41b6c4", "#f46d43", "#9e9ac8")),
      ggh4x.axis.nesttext.y = element_text(colour = c("#74c476", "#41b6c4", "#f46d43", "#9e9ac8")),
      legend.background = element_roundrect(color = "#969696"),
      panel.border = element_rect(size = 0.5),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
      axis.text = element_text(color = "#000000", size = 11),
      axis.text.y = element_text(color = rep(c("#41ae76", "#225ea8", "#fc4e2a", "#88419d"),as.vector(table(plot_df$ONTOLOGY)))),
      axis.text.y.left = element_blank(),
      axis.ticks.length.y.left = unit(10, "pt"),
      axis.ticks.y.left = element_line(color = NA),
      axis.title = element_text(color = "#000000", size = 15),
      plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
    ) + 
    coord_cartesian(clip = "off") + 
    annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                           gp = gpar(col = "#969696", lwd = 1.5)),
                      xmin = unit(10, "native"),
                      xmax = unit(42.25, "native"),
                      ymin = unit(40.85, "native"),
                      ymax = unit(42.25, "native"))
  
  
  
  plot
  
  ggsave(filename = "replot2_GO_KEGG.pdf",
         plot = plot,
         height = 11,
         width = 12.5)
  
  ggsave(filename = "replot2_GO_KEGG.png",
         plot = plot,
         height = 11,
         width = 12.5)
  
  # save enrichment table results
  write.csv(deg_go_2,"replot2_sl_go.csv",row.names = F)
  write.csv(deg_kegg_2,"replot2_sl_kegg.csv",row.names = F)
}

# 2024年12月08日 更新
GO_KEGG_plot_bar <- function(go,kegg,color=NULL){
  library(tidyverse)
  library(ggh4x)
  #假设go和kegg是GO和KEGG富集分析的结果
  #示例：GO 分析
  go_top<-as.data.frame(go)%>%
    group_by(ONTOLOGY)%>%
    slice_head(n=5)%>%#获取每个GO分类的前四个
    arrange(p.adjust)%>%
    ungroup()%>%#取消分组
    dplyr::select(ONTOLOGY,everything())
  
  #示例：KEGG 分析
  kegg_top<-as.data.frame(kegg)%>%
    dplyr::arrange(p.adjust)%>%
    dplyr::slice(1:5)%>%
    dplyr::select(ID:Count)%>%
    dplyr::mutate(ONTOLOGY="KEGG")%>%
    dplyr::select(ONTOLOGY,everything())
  
  #合并GO和KEGG数据
  data<-rbind(go_top,kegg_top)%>%
    dplyr::mutate(ONTOLOGY=factor(ONTOLOGY,levels=rev(c("BP","CC","MF","KEGG")),ordered=TRUE))%>%
    dplyr::arrange(desc(ONTOLOGY),-log10(p.adjust))
  
  #清洗数据
  data$geneID<-gsub("/",", ",data$geneID)
  data$CountNumber<-data$Count/1e3
  
  #选择颜色
  if(is.null(color)){
    color<-rev(c("#EFA39F","#F7CB65","#A0D8EA","#66C2A5"))
  }
  
  #绘制图形
  plot<-ggplot(data)+
    geom_bar(aes(x=-log10(p.adjust),y=interaction(Description,ONTOLOGY),
                 fill=ONTOLOGY),stat="identity")+
    scale_fill_manual(values=color,name="ONTOLOGY")+
    geom_text(aes(x=0.1,y=interaction(Description,ONTOLOGY),
                  label=Description),size=3,hjust=0,color="black")+
    geom_text(aes(x=0.1,y=interaction(Description,ONTOLOGY),
                  label=geneID),size=2,hjust=0,vjust=2.5,color="black")+
    geom_point(aes(x=-0.5,y=interaction(Description,ONTOLOGY),
                   size=Count,fill=ONTOLOGY),shape=21)+
    geom_text(aes(x=-0.5,y=interaction(Description,ONTOLOGY),label=Count),size=3)+
    scale_size(range=c(4,8),guide=guide_legend(override.aes=list(fill="black")))+
    #scale_x_continuous(expand = expansion(mult=c(0,0.2)),limits=c(-1,10))+
    guides(y = "axis_nested",
           y.sec = guide_axis_manual(breaks = 1:nrow(data),
                                     labels = data$Description)) + 
    #guides(fill=guide_legend(reverse=TRUE))+
    labs(x="-log10(FDR)",y="Description")+
    theme(
      legend.title=element_text(color="#000080",size=12),
      legend.text=element_text(color="#000000",size=10),
      axis.text.x=element_text(color="#000000",size=12),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title=element_text(color="#000080",size=14),
      panel.background=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      legend.background=element_blank(),
      ggh4x.axis.nestline.y = element_line(size = 2,color  = color),
      ggh4x.axis.nesttext.y = element_text(colour = color,hjust = 0.5,size = 12,angle = 90)
    )
  print(plot)
  ggsave(plot=plot,filename="Newbarplot.pdf",height=7,width=9)
  ggsave(plot=plot,filename="Newbarplot.png",height=7,width=9)
}


# enrich_GO = read.csv("./go_res.csv",row.names = 1)
# enrich_KEGG = read.csv("./KEGG_result.csv")
# library(RColorBrewer)
# #需要指定四种颜色：分别对应GO的三种类型和KEGG：
# col_value<-brewer.pal(4,"Accent")
# plot<-GO_KEGG_plot(enrich_GO,enrich_KEGG,color=col_value)
# ggsave(plot=plot,filename="plot.pdf",height=7,width=7)





Col_list = function(n = col_num){
  
  colorSpace <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                  "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", 
                  "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99", 
                  "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
                  "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#DCF0B9", 
                  "#8DD3C7", "#999999")
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}

