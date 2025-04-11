#####涉及Extended Data Fig.9o-9q
############该分析主要是为探究H4K77ac 与H3K27ac 增强子的相关性  因此找了公共数据库中的H3k27ac HSC的 CHIP-seq数据进行了联合分析
###参考文献提供了H3K27ac相关的增强子的list 因此读取list直接使用
##list的位置：/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/enchri_peaks_K27.csv
###参考文献网址：https://aacrjournals.org/cancerdiscovery/article/9/8/1080/42233/Aging-Human-Hematopoietic-Stem-Cells-Manifest

library("ChIPseeker")
library("org.Mm.eg.db")
library("clusterProfiler")
library("GenomicFeatures")
library(corrplot)
library(enrichplot)
library(ReactomePA)
library(dplyr)
library("TxDb.Mmusculus.UCSC.mm10.knownGene")

get_gsea_file <- function(
    object, 
    idents_vars = NULL,
    cls_name = 'cls.cls', 
    exp_name = 'matrix.gct',
    gseapy = T, 
    my_sep = '\t', 
    add_prefix=T, 
    prefix = 'C_'
){  
    write_cls <- function(x, object, cls_name){
      sink(file = cls_name)
      cat(ncol(object), my_sep, length(unique(x)), my_sep, 1, '\n', sep = '')
      cat('#', my_sep, paste(unique(x), collapse = my_sep), '\n', sep = '')
      cat(paste(x, collapse = my_sep), '\n', sep = '')
      sink()
    }
    if(!is.null(idents_vars)){
      Idents(object) <- object[[idents_vars]]
    }    
    clu_num <- length(levels(object)) 
    list <- NULL
    clu_name <- NULL
    lbl <- NULL
    for( i in 1:clu_num) {
      list2 <- WhichCells(object, idents = levels(object)[i])
      list <- if(is.null(list)) {list2} else {c(list, list2)}
      lbl2 <- if(add_prefix){
        rep(paste0('C_', levels(object)[i]), length(list2))
      } else {
        rep(levels(object)[i], length(list2))
      }
      lbl <- if(is.null(lbl)) {lbl2} else {c(lbl, lbl2)}
    }
    matrix <- as.data.frame(object$RNA@data[, list])
    matrix <- cbind(data.frame(Name = rownames(matrix), Description = rep("", nrow(matrix))), matrix)
    if(gseapy){
      cls <- matrix(lbl, nrow=1, ncol=length(lbl))
      write.table(cls, file = cls_name, quote = F, col.names = F, row.names = F, sep = my_sep)
    }else {
       write_cls(x = lbl, object = object, cls_name = cls_name)
    }    
    write.table(matrix, quote = F, sep = my_sep, row.names = F, file = exp_name)
}

Young<-readPeakFile("/data4/cut_TAG/H4K77ac/bedgraph_callpeakSEACR/newcall/SEACR_005/annobed/only_Y_peaks_005.bed")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
only_Y_peakAnno <- annotatePeak(Young, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoPie(only_Y_peakAnno)
only_Y_adata<-as.data.frame(only_Y_peakAnno)
only_Y_adata$anno<-""
for (i in 1:length(rownames(only_Y_adata))){
    a<-only_Y_adata$annotation[i]
    b<-unlist(strsplit(a,''))
    n<-length(which(b=="("))
    if (n==0){
        only_Y_adata$anno[i]<-a
    }else{
        n1<-which(b=="(")
        c<-b[1:n1-1]
        only_Y_adata$anno[i]<-c(paste(c,collapse=''))
    }
}
only_Y_gene_TSS<-unique(only_Y_adata$SYMBOL[only_Y_adata$anno=='Promoter '])
only_Y_gene<-unique(only_Y_adata$SYMBOL)
length(only_Y_gene)
length(only_Y_gene_TSS)

df_K27<-read.csv("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/enchri_peaks_K27.csv")

K27_Y_up_enhancer<-df_K27[(df_K27$H3K27ac_FC)<(-2)&
                         df_K27$EnhancerType=="Active","gene_symbol"]

library(Hmisc)
K27_Y_up_enhancer_mm<-capitalize(tolower(K27_Y_up_enhancer))

library(Seurat)
HSC<-readRDS("/data1/02.private/other/11_5/object2.rds")

K27_only_gene<-setdiff(K27_Y_up_enhancer_mm,only_Y_gene_TSS)

K77_only_gene<-setdiff(only_Y_gene_TSS,K27_Y_up_enhancer_mm)

length(K27_only_gene)
length(K77_only_gene)

genelist_GSEA<-unique(K27_only_gene)
ff_HSC<-HSC[rownames(HSC)%in%genelist_GSEA,]
ff_HSC

genelist_GSEA<-unique(K77_only_gene)
ff_HSC<-HSC[rownames(HSC)%in%genelist_GSEA,]
ff_HSC

length(overlap_gene)

genelist_GSEA<-unique(K27_only_gene)
ff_HSC<-HSC[rownames(HSC)%in%genelist_GSEA,]
ff_HSC

get_gsea_file(ff_HSC, idents_vars = 'time', 
              exp_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/input/K27_only_gene_matrix.txt', 
              cls_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/input/K27_only_gene_cls.cls',gseapy = F, add_prefix = F)

genelist_GSEA<-unique(c(K27_only_gene,overlap_gene))
ff_HSC<-HSC[rownames(HSC)%in%genelist_GSEA,]
ff_HSC

get_gsea_file(ff_HSC, idents_vars = 'time', 
              exp_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/input/K27_only_and_overlap_gene_matrix.txt', 
              cls_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/input/K27_only_and_overlap_gene_cls.cls',gseapy = F, add_prefix = F)

genelist_GSEA<-unique(K77_only_gene)
ff_HSC<-HSC[rownames(HSC)%in%genelist_GSEA,]
ff_HSC

get_gsea_file(ff_HSC, idents_vars = 'time', 
              exp_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/input/K77_only_gene_matrix.txt', 
              cls_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/input/K77_only_gene_cls.cls',gseapy = F, add_prefix = F)

genelist_GSEA<-unique(c(K77_only_gene,overlap_gene))
ff_HSC<-HSC[rownames(HSC)%in%genelist_GSEA,]
ff_HSC

get_gsea_file(ff_HSC, idents_vars = 'time', 
              exp_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/input/K77_only_and_overlap_gene_matrix.txt', 
              cls_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/input/K77_only_and_overlap_gene_cls.cls',gseapy = F, add_prefix = F)

library(Seurat)

overlap_gene<-intersect(K27_Y_up_enhancer_mm,only_Y_gene_TSS)

Idents(HSC)<-"time"

genelist_GSEA<-unique(overlap_gene)
ff_HSC<-HSC[rownames(HSC)%in%genelist_GSEA,]

#####304是Extended Data Fig.9o  三部分overlap的基因数量  因此知道这个数量的获取办法后，就知道了三部分分别是什么，通过计算就可得到其余部分数量
ff_HSC

plan("multisession", workers = 16)
options(future.globals.maxSize= 4*1024^4)
markers_features<-FindMarkers(HSC,`ident.1` = "2M",`ident.2` = "28M",logfc.threshold =log2(1),min.pct = 0,pseudocount.use = 0,features =rownames(ff_HSC) )

length(intersect(rownames(HSC),only_Y_gene_TSS))

# 加载包、读取数据
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(tidyverse)


markers_features$adj_p_val<-markers_features$p_val_adj

markers_features$adj_p_val[markers_features$p_val_adj==0]<-min(markers_features$p_val_adj[markers_features$p_val_adj!=0])

degdata<-markers_features


degdata$status <- 'stable'
degdata$status[degdata$avg_log2FC>log2(1.5) & degdata$p_val_adj<0.05] <- 'down'
degdata$status[degdata$avg_log2FC<(-log2(1.5)) & degdata$p_val_adj<0.05] <- 'up'

labeldata <- subset(degdata, abs(avg_log2FC)>log2(1.5) & p_val_adj<0.05)

id <- order(labeldata$avg_log2FC,decreasing = T)

labeldata <- labeldata[id,]


degdata$status<-factor(degdata$status,levels=c("up","stable","down"))

labeldata$gene<-rownames(labeldata)

options(repr.plot.width=12, repr.plot.height=8)
# (1)基础版-------------------------------------------
ggplot(degdata,aes(-avg_log2FC,-log10(adj_p_val)))+
  # 绘制散点
  geom_point(aes(color=status),
             size=5)+
  # 绘制垂直线
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)), linetype = "dashed")+
  # 绘制水平线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  theme_bw()+
geom_text_repel(data = labeldata[c(1:10,159:168),], 
                  aes(label = gene),
                  size=8)+xlab("Log2FC")+
  # 设置点的颜色
  scale_color_manual(values = c("down"='#84E0FF',"stable"='grey',"up"='#FFC284'))+
NULL
ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/plot/overlap_gene_dotplot.svg",height =8 ,width = 12,dpi=600)

library(openxlsx)
library(ggplot2)

df_pathway<-read.xlsx("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/GSEA/out/GSEA_GO.xlsx",sheet="Sheet1")

name<-c('Positive regulation of transcription by RNA polymerase II',
'Positive regulation of RNA metabolic process',
'Positive regulation of nucleobase containing compound metabolic process',
'Negative regulation of biosynthetic process',
'Nucleolus',
'Negative regulation of gene expression',
'Positive regulation of cellular biosynthetic process',
'Positive regulation of cellular component organization',
'Regulation of transferase activity',
'Phosphorylation',
'Nucleoside triphosphatase regulator activity',
'Response to growth factor',
'Chromosome',
'Posttranscriptional regulation of gene expression',
'Cellular response to endogenous stimulus',
'Negative regulation of nucleobase containing compound metabolic process',
'Organic acid metabolic process',
'DNA binding transcription factor activity',
'Nuclear body',
'Transcription regulator activity')

df_pathway$name<-name
df_pathway<-df_pathway[1:15,]
df_pathway<-df_pathway[order(df_pathway$SIZE),]

df_pathway$NOM_p_val

options(repr.plot.width=18, repr.plot.height=10)
ggplot(df_pathway,aes(x=NES,y=NAME))+
    #geom_segment( aes(x=0, xend=NES, y=NAME,yend=NAME)) +
    geom_point( size=8,aes(color=NOM_p_val,fill=NOM_p_val), shape=21, stroke=2)+
scale_y_discrete(limits = rev(df_pathway$NAME), labels = rev(name))+
#scale_colour_gradientn(colours = c("#D26B44","white","blue"),breaks=c(0,0.01,0.05))+
theme(axis.text.x = element_text(size=20),
          axis.text.y = element_text(size = 26),
          axis.title.y=element_blank(),
          axis.title.x=element_text(size=20),
          panel.background = element_rect(color = "black",fill="white",linetype="solid"),
          axis.line=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          #panel.grid.major=element_blank(),
          #panel.grid.minor=element_blank(),
          legend.text = element_text(size = 20), 
          legend.title = element_text(size =20),plot.title=element_text(hjust = 0.5,size=25),
          legend.position="right",
          #panel.grid.minor = element_line(colour = "red"), 
     )+
    guides(size=guide_legend(reverse = T,ncol=1))+
    xlab("NES")+
scale_fill_gradientn("NOM_p_val",colours = c("#008EBD","#00ACE6","#AEEBFF","white"),limits=c(0,0.05),breaks=c(0,0.05))+
scale_color_gradientn("NOM_p_val",colours = c("#008EBD","#00ACE6","#AEEBFF","white"),limits=c(0,0.05),breaks=c(0,0.05))+
#scale_size_continuous(range = c(8,4),guide = FALSE,limits=c(0,0.1),breaks=c(0,0.05,0.1))+
    #scale_x_continuous(limits =c(0,0.06) , breaks=c(0,75,150,225))+
NULL
ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/plot/K77ac_K27_ovelap.svg",height =10 ,width = 20,dpi=600)
