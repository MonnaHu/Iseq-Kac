##运行环境：/home/hulab/.conda/envs/Gpb_4.2.0
######Fig2g 艳秋师姐进行前期的分析 之后把结果给我 我进行可视化 下列为可视化代码 
library(ggplot2)
library(openxlsx)

HSC_pathway<-read.xlsx("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/choose_pathway_20241231.xlsx",sheet = "HSC")

HSC_pathway_ft<-HSC_pathway

HSC_pathway_ft$type<-factor(HSC_pathway_ft$type,levels=c("ref_HSC","our_HSC"))


options(repr.plot.width=12, repr.plot.height=6)
ggplot(HSC_pathway_ft,aes_string(x="NES",y="NAME",fill="type"))+
    geom_bar(stat='identity',width = 0.65,position = position_dodge(0.65))+
scale_fill_manual(values = c("our_HSC"="#94271F","ref_HSC"="#DC885C"))+
theme(
          axis.text.y = element_text(size=15),
          axis.title=element_blank(),
          panel.background = element_rect(color = "black",fill="white",linetype="solid"),
          axis.line=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.text = element_text(size = 10), 
          legend.title = element_text(size =0),plot.title=element_text(hjust = 0.5,size=25),legend.position="NA")+
geom_text(aes(label=SIZE), vjust = 0.5,hjust=0.8, position = position_dodge(.5),size=5,angle=360)+

NULL

HSC_pathway_ft$type<-factor(HSC_pathway_ft$type,levels=c("ref_HSC","our_HSC"))


options(repr.plot.width=3, repr.plot.height=5)
ggplot(HSC_pathway_ft,aes_string(x="NES",y="NAME",fill="type"))+
    geom_bar(stat='identity',width = 0.8,position = position_dodge(0.8))+
scale_fill_manual(values = c("our_HSC"="#94271F","ref_HSC"="#DC885C"))+
theme(
          axis.text = element_blank(),
          axis.title=element_blank(),
          panel.background = element_rect(color = "black",fill="white",linetype="solid"),
          axis.line=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.text = element_text(size = 10), 
          legend.title = element_text(size =0),plot.title=element_text(hjust = 0.5,size=25),legend.position="NA")+
NULL
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/plot/pathway_Plot/HSC_pathway_protein_cor.svg",width=3,height = 5,dpi=300)

MPP1_pathway<-read.xlsx("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/choose_pathway_20241231.xlsx",sheet = "MPP1")

MPP1_pathway_ft<-MPP1_pathway

MPP1_pathway_ft$type<-factor(MPP1_pathway_ft$type,levels=c("ref_MPP1","our_MPP1"))

options(repr.plot.width=10, repr.plot.height=4)
ggplot(MPP1_pathway_ft[order(MPP1_pathway_ft$type,decreasing = F),],aes(x=NES,y=NAME,fill=type))+
    geom_bar(stat='identity',width = 0.65,position = position_dodge(0.65))+
scale_fill_manual(values = c("our_MPP1"="#67A158","ref_MPP1"="#AEC9BA"))+
theme(
          axis.text.y = element_text(size=15),
          axis.title=element_blank(),
          panel.background = element_rect(color = "black",fill="white",linetype="solid"),
          axis.line=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.text = element_text(size = 10), 
          legend.title = element_text(size =0),plot.title=element_text(hjust = 0.5,size=25),legend.position="NA")+
geom_text(aes(label=SIZE), vjust = 0.5,hjust=0.8, position = position_dodge(.5),size=5,angle=360)+

NULL



options(repr.plot.width=5, repr.plot.height=5)
ggplot(MPP1_pathway_ft[order(MPP1_pathway_ft$type,decreasing = F),],aes(x=NES,y=NAME,fill=type))+
    geom_bar(stat='identity',width = 1,position = position_dodge(0.8))+
scale_fill_manual(values = c("our_MPP1"="#67A158","ref_MPP1"="#AEC9BA"))+
theme(
          axis.text = element_blank(),
          axis.title=element_blank(),
          panel.background = element_rect(color = "black",fill="white",linetype="solid"),
          axis.line=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.text = element_text(size = 10), 
          legend.title = element_text(size =0),plot.title=element_text(hjust = 0.5,size=25),legend.position="NA")+
 scale_x_continuous(breaks = c(-2,-1,0,1,2),limits=c(-2,2))+
NULL
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/plot/pathway_Plot/MPP1_pathway_protein_cor.svg",width=5,height = 5,dpi=300)

MPP2_pathway<-read.xlsx("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/choose_pathway_20241231.xlsx",sheet = "MPP2")

MPP2_pathway_ft<-MPP2_pathway

MPP2_pathway_ft$type<-factor(MPP2_pathway_ft$type,levels=c("ref_MPP2","our_MPP2"))

options(repr.plot.width=10, repr.plot.height=4)
ggplot(MPP2_pathway_ft[order(MPP2_pathway_ft$type,decreasing = F),],aes(x=NES,y=NAME,fill=type))+
    geom_bar(stat='identity',width = 0.65,position = position_dodge(0.65))+
scale_fill_manual(values = c("our_MPP2"="#7CA7B6","ref_MPP2"="#B6CDDC"))+
theme(
          axis.text.y = element_text(size=15),
          axis.title=element_blank(),
          panel.background = element_rect(color = "black",fill="white",linetype="solid"),
          axis.line=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.text = element_text(size = 10), 
          legend.title = element_text(size =0),plot.title=element_text(hjust = 0.5,size=25),legend.position="NA")+
geom_text(aes(label=SIZE), vjust = 0.5,hjust=0.8, position = position_dodge(.5),size=5,angle=360)+

NULL



options(repr.plot.width=3, repr.plot.height=5)
ggplot(MPP2_pathway_ft[order(MPP2_pathway_ft$type,decreasing = F),],aes(x=NES,y=NAME,fill=type))+
    geom_bar(stat='identity',width = 1,position = position_dodge(0.8))+
scale_fill_manual(values = c("our_MPP2"="#7CA7B6","ref_MPP2"="#B6CDDC"))+
theme(
          axis.text = element_blank(),
          axis.title=element_blank(),
          panel.background = element_rect(color = "black",fill="white",linetype="solid"),
          axis.line=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.text = element_text(size = 10), 
          legend.title = element_text(size =0),plot.title=element_text(hjust = 0.5,size=25),legend.position="NA")+
NULL
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/plot/pathway_Plot/MPP2_pathway_protein_cor.svg",width=3,height = 5,dpi=300)

MPP3_pathway<-read.xlsx("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/choose_pathway_20241231.xlsx",sheet = "MPP3")

MPP3_pathway_ft<-MPP3_pathway

MPP3_pathway_ft$type<-factor(MPP3_pathway_ft$type,levels=c("ref_MPP3","our_MPP3"))

options(repr.plot.width=11, repr.plot.height=4)
ggplot(MPP3_pathway_ft[order(MPP3_pathway_ft$type,decreasing = F),],aes(x=NES,y=NAME,fill=type))+
    geom_bar(stat='identity',width = 0.65,position = position_dodge(0.65))+
scale_fill_manual(values = c("our_MPP3"="#7F5A8D","ref_MPP3"="#C5C5E7"))+
theme(
          axis.text.y = element_text(size=15),
          axis.title=element_blank(),
          panel.background = element_rect(color = "black",fill="white",linetype="solid"),
          axis.line=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.text = element_text(size = 10), 
          legend.title = element_text(size =0),plot.title=element_text(hjust = 0.5,size=25),legend.position="NA")+
geom_text(aes(label=SIZE), vjust = 0.5,hjust=0.8, position = position_dodge(.5),size=5,angle=360)+

NULL



options(repr.plot.width=3, repr.plot.height=5)
ggplot(MPP3_pathway_ft[order(MPP3_pathway_ft$type,decreasing = F),],aes(x=NES,y=NAME,fill=type))+
    geom_bar(stat='identity',width = 1,position = position_dodge(0.8))+
scale_fill_manual(values = c("our_MPP3"="#7F5A8D","ref_MPP3"="#C5C5E7"))+
theme(
          axis.text = element_blank(),
          axis.title=element_blank(),
          panel.background = element_rect(color = "black",fill="white",linetype="solid"),
          axis.line=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.text = element_text(size = 10), 
          legend.title = element_text(size =0),plot.title=element_text(hjust = 0.5,size=25),legend.position="NA")+
NULL
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/plot/pathway_Plot/MPP3_pathway_protein_cor.svg",width=3,height = 5,dpi=300)

##运行环境：/home/hulab/.conda/envs/Gpb_4.2.0
######Fig 3h  第一次投稿（日期：20240311）是我画的 之后这次投稿（日期：20250204）是我把代码发给艳秋师姐 她自己画的 
###下列整理的第一次投稿的可视化代码 
library(ggplot2)
library(enrichplot)
library(data.table)
library(DOSE)
library(RColorBrewer)
library(ggrepel)
gsea_plot <- function(
    gene_set_dir = 'edb/', 
    result_folder = '/data1/02.private/dengyj/analysis/thymus/clustering_10_20/ETP_DN/Enrichment/nov16/12_vs_60.Gsea.1605493817058', 
    set_name = "KEGG_RIBOSOME",
    title = set_name,
    color="#87c444",
    base_size = 11,
    rel_heights=c(1.5, .5, 1),
    subplots = 1:3,
    pvalue_table = FALSE,
    ES_geom="line",
    color.vline="NA",
    color.hline="gray67",
    line_size = 2, 
    label_a = 'A',
    label_b = 'B',
    title_size = 15,
    val_x = c(0.25, 0.25), 
    val_y = c(0.725, 0.675),  
    label_x = c(0.19, 0.88),
    label_y = c(0.315, 0.315),
    label_size = 5, 
    val_size = 5
){
   if(!stringr::str_detect(result_folder, '/$')){
       result_folder <- paste0(result_folder, '/')
   }
    require(ggplot2)
    gene_set_dir <- paste0(result_folder, gene_set_dir)
    gene_set_dir_files <- dir(gene_set_dir)
    gene_set <- paste0(gene_set_dir, gene_set_dir_files[grep('gmt$', gene_set_dir_files)])
    set_name <- toupper(set_name)
    file <- dir(result_folder)
    report_file <- file[grep('(^gsea_report){1}(.*)xls$', file)]
 
    result_1 <- as.data.frame(data.table::fread(paste0(result_folder, "/",report_file[1]),header = T))
    rownames(result_1)<-result_1$NAME
    result_1<-result_1[,-1]
    result_2 <- as.data.frame(data.table::fread(paste0(result_folder, "/",report_file[2]),header = T))
    rownames(result_2)<-result_2$NAME
    result_2<-result_2[,-1]
    result <- rbind(result_1, result_2)
    line <- strsplit(readLines(gene_set), '\t')
    location <- which(sapply(line, function(x) x[1]) == set_name)
    gene <- line[[location]][3:length(line[[location]])]
    rank_dir <- file[grep('^ranked_gene_list', file)]
    rank_list <- data.table::fread(paste0(result_folder, rank_dir))
    score <- rank_list$SCORE
    names(score) <- as.character(as.matrix(rank_list$NAME))

    ES_geom <- match.arg(ES_geom, c("line", "dot"))

    geneList <- position <- NULL ## to satisfy codetool

    #gsinfo函数处理
    gseaScores <- getFromNamespace("gseaScores", "DOSE")
    gsdata <- gseaScores(geneList=score, geneSet=gene, exponent=1, fortify=TRUE)
    gsdata$ymin=0
    gsdata$ymax=0
    pos <- gsdata$position == 1
    h <- diff(range(gsdata$runningScore))/20
    gsdata$ymin[pos] <- -h
    gsdata$ymax[pos] <- h
    gsdata$geneList <- score
    gsdata$Description <- "extracellular matrix"

    p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
            theme_classic(base_size) +
            theme(panel.grid.major = element_line(colour = "grey92"),
                panel.grid.minor = element_line(colour = "grey92"),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank()) +
            scale_x_continuous(expand=c(0,0))

        if (ES_geom == "line") {
            es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description), size=line_size)
        } else {
            es_layer <- geom_point(aes_(y = ~runningScore, color= ~Description), size=1, data = subset(gsdata, position == 1))
        }

        p.res <- p + es_layer +
            theme(legend.position = c(.8, .8), legend.title = element_blank(),
                legend.background = element_rect(fill = "transparent"))

        p.res <- p.res + ylab("Running Enrichment Score") +
            theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.line.x=element_blank(),
                plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))

        i <- 0
        for (term in unique(gsdata$Description)) {
            idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
            gsdata[idx, "ymin"] <- i
            gsdata[idx, "ymax"] <- i + 1
            i <- i + 1
        }
        p2 <- ggplot(gsdata, aes_(x = ~x)) +
            geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
            xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
            theme(legend.position = "none",
                plot.margin = margin(t=-.1, b=0,unit="cm"),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                axis.line.x = element_blank()) +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0))


            v <- seq(1, sum(gsdata$position), length.out=9)
            inv <- findInterval(rev(cumsum(gsdata$position)), v)
            if (min(inv) == 0) inv <- inv + 1

            col=c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))

            ymin <- min(p2$data$ymin)
            yy <- max(p2$data$ymax - p2$data$ymin) * .3
            xmin <- which(!duplicated(inv))
            xmax <- xmin + as.numeric(table(inv)[unique(inv)])
            d <- data.frame(ymin = ymin, ymax = yy,
                            xmin = xmin,
                            xmax = xmax,
                            col = col[unique(inv)])
            p2 <- p2 + geom_rect(
                        aes_(xmin=~xmin,
                                xmax=~xmax,
                                ymin=~ymin,
                                ymax=~ymax,
                                fill=~I(col)),
                        data=d,
                        alpha=.9,
                        inherit.aes=FALSE)


        df2 <- p$data #data.frame(x = which(p$data$position == 1))
        df2$y <- p$data$geneList[df2$x]
        p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0), color="grey")
        p.pos <- p.pos + ylab("Ranked list metric") +
            xlab("Rank in Ordered Dataset") +
            theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))

        if (!is.null(title) && !is.na(title) && title != "")
            p.res <- p.res + ggtitle(title) + theme(plot.title=element_text(size=title_size, hjust = 0.5))

            p.res <- p.res + scale_color_manual(values=color)
            if (length(color) == 1) {
                p.res <- p.res + theme(legend.position = "none")
                p2 <- p2 + scale_color_manual(values = "black")
            } else {
                p2 <- p2 + scale_color_manual(values = color)
            }

    enrichmentScore <- if(abs(max(gsdata$runningScore)) >= abs(min(gsdata$runningScore))) {max(gsdata$runningScore)} else {min(gsdata$runningScore)}
    es.df <- data.frame(es = which.min(abs(p$data$runningScore - enrichmentScore)))
    p.res <- p.res + geom_vline(data = es.df, aes_(xintercept = ~es), colour = color.vline, linetype = "dashed") + geom_hline(aes(yintercept=0), colour=color.hline, linetype="dashed")

        plotlist <- list(p.res, p2, p.pos)[subplots]
        n <- length(plotlist)
        plotlist[[n]] <- plotlist[[n]] +
            theme(axis.line.x = element_line(),
                axis.ticks.x=element_line(),
                axis.text.x = element_text())

        if (length(subplots) == 1)
            return(plotlist[[1]] + theme(plot.margin=margin(t=.2, r = .2, b=.2, l=.2, unit="cm")))


        if (length(rel_heights) > length(subplots))
            rel_heights <- rel_heights[subplots]
    pic <- plot_grid(plotlist = plotlist, ncol=1, align="v", rel_heights=rel_heights)
    val_to_plot <- data.frame(x=val_x,y=val_y,
                               col=c("black","black"),
                               text=c(paste0("NES:", round(result[set_name, "NES"],digits=3)), 
                               paste0("italic('P')~val:", round(result[set_name, "NOM p-val"],digits=3)) 
                               ))
    label_to_plot <- data.frame(x = label_x, y = label_y, col = c('red', 'blue'), text = c(label_a, label_b))                           
    pic <- pic + geom_text(data=val_to_plot,mapping = aes(x, y, color = col, label = text),parse = T, 
                           size = val_size)+
           theme(legend.position = "none") + 
           geom_text_repel(parse = TRUE) +
           geom_text(data = label_to_plot, mapping = aes(x,y, color = col, label = text), size = label_size) + 
           scale_colour_manual(values = c('red' = 'red', 'blue' = 'blue', 'black' = "black"))
    return(pic)
}
                             
                             

library(cowplot)

library(openxlsx)
HSC_pathway<-read.xlsx("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/PreRankedGSEA/Preranked-gsea-LSK/HSC_choose_pathway.xlsx")

MPP1_pathway<-read.xlsx("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/PreRankedGSEA/Preranked-gsea-LSK/MPP1_choose_pathway.xlsx")
MPP2_pathway<-read.xlsx("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/PreRankedGSEA/Preranked-gsea-LSK/MPP2_choose_pathway.xlsx")
MPP3_pathway<-read.xlsx("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/PreRankedGSEA/Preranked-gsea-LSK/MPP3_choose_pathway.xlsx")

for (i in (1:length(HSC_pathway$NAME))){
    set_name<-HSC_pathway$NAME[i]
    result_folder<-HSC_pathway$file[i]
    p<-gsea_plot(    
    result_folder = result_folder, 
    set_name = set_name, label_a = 'pos', label_b = 'neg', 
    label_y = c(0.315, 0.1))
    print(p)
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/GSEA_Plot/ATAC_GSEA/GO_POSITIVE_REGULATION_OF_NOTCH_SIGNALING_PATHWAY.svg",height = 6,width=6,dpi=300)
    ggsave(paste0("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/GSEA_Plot/HSC/",set_name,".svg"),height = 6,width=6,dpi=300)

}

for (i in (1:length(MPP1_pathway$NAME))){
    set_name<-MPP1_pathway$NAME[i]
    result_folder<-MPP1_pathway$file[i]
    p<-gsea_plot(    
    result_folder = result_folder, 
    set_name = set_name, label_a = 'pos', label_b = 'neg', 
    label_y = c(0.315, 0.1))
    print(p)
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/GSEA_Plot/ATAC_GSEA/GO_POSITIVE_REGULATION_OF_NOTCH_SIGNALING_PATHWAY.svg",height = 6,width=6,dpi=300)
    ggsave(paste0("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/GSEA_Plot/MPP1/",set_name,".svg"),height = 6,width=6,dpi=300)

}

for (i in (1:length(MPP2_pathway$NAME))){
    set_name<-MPP2_pathway$NAME[i]
    result_folder<-MPP2_pathway$file[i]
    p<-gsea_plot(    
    result_folder = result_folder, 
    set_name = set_name, label_a = 'pos', label_b = 'neg', 
    label_y = c(0.315, 0.1))
    print(p)
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/GSEA_Plot/ATAC_GSEA/GO_POSITIVE_REGULATION_OF_NOTCH_SIGNALING_PATHWAY.svg",height = 6,width=6,dpi=300)
    ggsave(paste0("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/GSEA_Plot/MPP2/",set_name,".svg"),height = 6,width=6,dpi=300)

}

for (i in (1:length(MPP3_pathway$NAME))){
    set_name<-MPP3_pathway$NAME[i]
    result_folder<-MPP3_pathway$file[i]
    p<-gsea_plot(    
    result_folder = result_folder, 
    set_name = set_name, label_a = 'pos', label_b = 'neg', 
    label_y = c(0.315, 0.1))
    print(p)
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/GSEA_Plot/ATAC_GSEA/GO_POSITIVE_REGULATION_OF_NOTCH_SIGNALING_PATHWAY.svg",height = 6,width=6,dpi=300)
    ggsave(paste0("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/GSEA_Plot/MPP3/",set_name,".svg"),height = 6,width=6,dpi=300)

}

###############Fig.4a-4b
###Supplementary Table 7b 7c peaks注释列表
##callpeaks 方法：SEACR Top5%
#############运行环境：R_base 3.6
library("ChIPseeker")
library("org.Mm.eg.db")
library("clusterProfiler")
library("GenomicFeatures")
library(corrplot)
library(enrichplot)
library(ReactomePA)
library(dplyr)
library("TxDb.Mmusculus.UCSC.mm10.knownGene") 

mypie <- function (data, area_var, fill_var, pie_colors = "black", pie_size = 0.5, 
    plot_label = T, labels = NULL, digits = 2, min_cutoff = 0, 
    lbl_adjust = 0, lbl_size = 5, lbl_face = "plain", lbl_colors = "black", 
    lbl_family = "") 
{
    if (class(data[, fill_var]) != "factor") {
        data[, fill_var] <- factor(data[, fill_var])
    }
    pie_colors <- match.arg(pie_colors)
    data$x_pos <- 1
    .y <- data[, area_var][match(levels(data[, fill_var]), as.character(data[, 
        fill_var]))]
    .y <- rev(.y)
    data$lbl_ypos <- cumsum(.y) - .y/2
    p <- ggplot(data, aes(x_pos)) + geom_bar(stat = "identity", 
        mapping = aes_string(y = area_var, fill = fill_var), 
        color = pie_colors, size = pie_size) + coord_polar(theta = "y", 
        start = 0) + theme(axis.title = element_blank(), axis.text = element_blank(), 
        axis.line = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = NA))
    if (plot_label) {
        if (is.null(labels)) {
            labels <- paste0(sprintf(paste0("%0.", digits, "f"), 
                .y/sum(.y) * 100), "%")
            labels[.y/sum(.y) <= min_cutoff] <- ""
        }
        else {
            if (class(labels) != "character" | length(labels) != 
                ncol(data)) {
                stop("Invalid input: labels!")
            }
        }
        data <- cbind(data, labels = labels)
        p <- p + geom_text(data = data, mapping = aes(x = x_pos + 
            lbl_adjust, y = lbl_ypos, label = labels), size = lbl_size, 
            fontface = lbl_face, colour = lbl_colors, family = lbl_family)
    }
    return(p)
}



all_Old<-readPeakFile("/data4/cut_TAG/H4K77ac/bedgraph_callpeakSEACR/newcall/SEACR_005/Old_005.rmbl.bed")
all_Young<-readPeakFile("/data4/cut_TAG/H4K77ac/bedgraph_callpeakSEACR/newcall/SEACR_005/You_005.rmbl.bed")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
all_OldpeakAnno <- annotatePeak(all_Old, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
all_YoungpeakAnno <- annotatePeak(all_Young, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
all_Y_adata<-as.data.frame(all_YoungpeakAnno)
all_O_adata<-as.data.frame(all_OldpeakAnno)

all_Y_adata$anno<-""
for (i in 1:length(rownames(all_Y_adata))){
    a<-all_Y_adata$annotation[i]
    b<-unlist(strsplit(a,''))
    n<-length(which(b=="("))
    if (n==0){
        all_Y_adata$anno[i]<-a
    }else{
        n1<-which(b=="(")
        c<-b[1:n1-1]
        all_Y_adata$anno[i]<-c(paste(c,collapse=''))
    }
}

all_O_adata$anno<-""
for (i in 1:length(rownames(all_O_adata))){
    a<-all_O_adata$annotation[i]
    b<-unlist(strsplit(a,''))
    n<-length(which(b=="("))
    if (n==0){
        all_O_adata$anno[i]<-a
    }else{
        n1<-which(b=="(")
        c<-b[1:n1-1]
        all_O_adata$anno[i]<-c(paste(c,collapse=''))
    }
}

df_num_of_sample  <- table(all_O_adata$anno) 

df <- df_num_of_sample/length(rownames(all_O_adata))

data <- data.frame(number = df, Identity = rownames(df))

data$Identity <- 
    factor(data$Identity, levels=c("Promoter ","3\' UTR","5\' UTR",'Exon ','Intron ','Downstream ','Distal Intergenic'))

data<-data[order(data$Identity),]

data$Llab<-paste0(data$Identity," (",round(data$number.Freq,4)*100,"%)")

library(ggplot2)

p <- mypie( data=data, area_var = 'number.Freq', min_cutoff = 0.01,
         fill_var = 'Identity',
         plot_label = F, lbl_adjust = 0.6,lbl_size = 10,digits = 1,
         pie_colors = 'black', lbl_face = 'bold', pie_size = 1)+
        
        scale_fill_manual(values = c("#95D2CA","#6867A9","#F18E89","#E86086","#0084BC","#D03366","#DAD379"),
                         labels=data$Llab)+
        
NULL

options(repr.plot.width=12, repr.plot.height=6)

p+guides(fill=guide_legend(override.aes = list(color = "white"),keywidth =2,keyheight =1))+
theme(
   legend.key.size = unit(35, "pt"),
    legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)
)

#ggsave("/data1/02.private/guanpb/data/Trim68/cut_tag/H4K77_1st/new analysis/plot/SEACR005_pie/K77_1st_Old_anno_SEACR005_pie_1231.svg",dpi=600,width=12,height=6)

df_num_of_sample  <- table(all_Y_adata$anno) 

df <- df_num_of_sample/length(rownames(all_Y_adata))

data <- data.frame(number = df, Identity = rownames(df))

data$Identity <- 
    factor(data$Identity, levels=c("Promoter ","3\' UTR","5\' UTR",'Exon ','Intron ','Downstream ','Distal Intergenic'))

data<-data[order(data$Identity),]

data$Llab<-paste0(data$Identity," (",round(data$number.Freq,4)*100,"%)")

library(ggplot2)

p <- mypie( data=data, area_var = 'number.Freq', min_cutoff = 0.01,
         fill_var = 'Identity',
         plot_label = F, lbl_adjust = 0.6,lbl_size = 10,digits = 1,
         pie_colors = 'black', lbl_face = 'bold', pie_size = 1)+
        
        scale_fill_manual(values = c("#95D2CA","#6867A9","#F18E89","#E86086","#0084BC","#D03366","#DAD379"),
                         labels=data$Llab)+
        
NULL

options(repr.plot.width=12, repr.plot.height=6)

p+guides(fill=guide_legend(override.aes = list(color = "white"),keywidth =2,keyheight =1))+
theme(
   legend.key.size = unit(35, "pt"),
    legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)
)

#ggsave("/data1/02.private/guanpb/data/Trim68/cut_tag/H4K77_1st/new analysis/plot/SEACR005_pie/K77_1st_Young_anno_SEACR005_pie_1231.svg",dpi=600,width=12,height=6)


library("ChIPseeker")
library("org.Mm.eg.db")
library("clusterProfiler")
library("GenomicFeatures")
library(corrplot)
library(enrichplot)
library(ReactomePA)
library(dplyr)
library("TxDb.Mmusculus.UCSC.mm10.knownGene")

sampleList<-c("Y1","Y2","Y3","O1","O2","O3")
alignResult = c()
for(i in sampleList){
  alignRes = read.table(paste0("/data4/cut_TAG/H4K77ac/new_analysis/alignment/sam/bowtie2_summary/", i, "_bowtie2.txt"), header = FALSE, fill = TRUE)
  alignRate = substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  histInfo = strsplit(i, "_")[[1]]
  alignResult = data.frame(Name = i,
                           SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric, 
                           MappedFragNum_mm10 = alignRes$V1[4] %>% as.character %>% as.numeric + alignRes$V1[5] %>% as.character %>% as.numeric, 
                           AlignmentRate_mm10 = alignRate %>% as.numeric)  %>% rbind(alignResult, .)
}
#alignResult$Histone = factor(alignResult$Name, levels = histList)
alignResult %>% mutate(AlignmentRate_mm10 = paste0(AlignmentRate_mm10, "%"))

peakN<-c()
peakWidth<-c()

for (i in sampleList){
    peakInfo<-read.table(paste0("/data4/cut_TAG/H4K77ac/new_analysis/SEACR_005/",i,"_005.rmbl.bed"),
                         header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
    
    peakN = data.frame(peakN = nrow(peakInfo) ,Name=i) %>% rbind(peakN, .)
    peakWidth = data.frame(width = peakInfo$width, Name = i )  %>% rbind(peakWidth, .)
}
peakN

peakOverlap = c() 
for(i in sampleList){
     overlap.gr = GRanges() 
     peakInfo = read.table(paste0("/data4/cut_TAG/H4K77ac/new_analysis/SEACR_005/",i,"_005.rmbl.bed"), header = FALSE, fill = TRUE)
      peakInfo.gr = GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
      if(length(overlap.gr) >0){
        overlap.gr = overlap.gr[findOverlaps(overlap.gr, peakInfo.gr)@from]
      }else{
        overlap.gr = peakInfo.gr
        
            }
    peakOverlap = data.frame(peakReprod = length(overlap.gr), Name =i) %>% rbind(peakOverlap, .)
}
peakReprod = left_join(peakN, peakOverlap, by = c("Name")) %>% mutate(peakReprodRate = peakReprod/peakN * 100)

library(chromVAR)

bamDir = c("/data4/cut_TAG/H4K77ac/new_analysis/rmdup/rmdupbam")
inPeakData = c()
for(i in sampleList){
    peakRes = read.table(paste0("/data4/cut_TAG/H4K77ac/new_analysis/SEACR_005/",i,"_005.rmbl.bed"), header = FALSE, fill = TRUE)
    peak.gr = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
    bamFile = paste0(bamDir, "/",i, ".sorted.markdup.bam")
    fragment_counts <- getCounts(bamFile, peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
    inPeakN = counts(fragment_counts)[,1] %>% sum
    inPeakData = rbind(inPeakData, data.frame(inPeakN = inPeakN, Name = i))
}

frip = left_join(inPeakData, alignResult, by = c("Name")) %>% 
    mutate(frip = inPeakN/MappedFragNum_mm10 * 100)
#frip = left_join(inPeakData, alignResult, by = c("Name")) %>% mutate(Frip = paste0(round(inPeakN/MappedFragNum_mm10*100 ,2),"%"))

frip

#############计算的逻辑  上述代码数据的inPeakN的例   除以 去重后的该样本对应的reads剩余数量
###去重后的该样本对应的reads数量 获取方式：
##/data4/cut_TAG/H4K77ac/txt 该文件夹下对应的样本名字 打开后ESTIMATED_LIBRARY_SIZE 对应的数字

276602/1109968
329634/1324735
151403/602282

240929/1037078
272555/1200048
189453/760603

library(ggplot2)
mydata_frip<-data.frame(
  name=c("Young","Young","Young","old","old","old"),
  value=c("0.2492","0.2488","0.2514","0.2323","0.2271","0.2491" )
)

t.test(value~name,mydata_frip)

##Fig.4c 
##callpeaks 方法：SEACR Top5%
##运行环境：/home/hulab/.conda/envs/Gpb_4.2.0
library(ggplot2)
mydata_frip<-data.frame(
  name=c("Young","Young","Young","old","old","old"),
  value=c("0.2491","0.2488","0.2514","0.2323","0.2271","0.2491" )
)
mydata_frip$value<-as.numeric(mydata_frip$value)
mydata_frip$name<-factor(mydata_frip$name,levels=c("Young","old"))
mydata_frip$value<-as.numeric(mydata_frip$value)
options(repr.plot.width=5, repr.plot.height=5)
ggplot(mydata_frip, aes(x=name,y=value,fill=name)) + 
scale_y_continuous(limits =c(0.22,0.26) , breaks=c(0.22,0.24,0.26))+
  geom_boxplot()+
scale_fill_manual(values = c("#33CCFF","#FF9933"))
NULL
#ggsave("H4K77ac_FRIP_box.svg",width=5,height=5) 

###############Fig.4d
#######基于Linux运行！！！！！！！！！！！！！！！！！！！！ 
####运行环境 base
##callpeaks 方法：MACS2 默认参数
####运行位置：/data4/cut_TAG/H4K77ac

#!/bin/bash

mkdir -p heatmap/mat
mkdir -p heatmap/plot

bedfile='/data4/cut_TAG/H4K77ac/use/bed'
bwfile='/data4/cut_TAG/H4K77ac/bw'

computeMatrix reference-point \
--referencePoint center \
-b 3000 -a 3000 \
-R ${bedfile}/merge.common.bed ${bedfile}/merge.Y.bed ${bedfile}/merge.O.bed \
-S  ${bwfile}/Y.sorted.markdup.bw ${bwfile}/O.sorted.markdup.bw \
-o heatmap/mat/COY_center.mat.gz -p 8 \
--outFileSortedRegions heatmap/mat/regions_COY_genes.bed

plotHeatmap -m heatmap/mat/COY_center.mat.gz \
--whatToShow 'heatmap and colorbar' \
--heatmapWidth 8 \
-out heatmap/plot/COY_heatmap.png \
--colorMap Blues --missingDataColor 1 \
--dpi 600 --regionsLabel 'common_peaks' 'Y_peaks' 'O_peaks'

plotProfile -m heatmap/mat/COY_center.mat.gz  \
-out heatmap/plot/COY_plot.png \
--colors blue "#f29438" "#4bc1ee" \
--dpi 600 --plotWidth 8 \
--regionsLabel "Common" "Only_Young" "Only_Old" \
--samplesLabel "Young.bw" "Old.bw"


####输出文件位置
##mat:/data4/cut_TAG/H4K77ac/use/old/bed/heatmap/mat
##plot:/data4/cut_TAG/H4K77ac/use/old/bed/heatmap/plot

#################Fig.4e 
#####Supplementary Table 7e. List of genes expressed 转录组基因表达列表
##callpeaks 方法：SEACR Top5%  only_peaks TSS±3Kb
#######上部分Venv的获取代码:<1.2.注释并联合转录组提取运算GSEA输入文件>
###整合转录组后 2Mvs28M 的GSEA
##运行环境：R_base

library(ggplot2)
df<-read.csv("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/GSEA.csv")
df<-df[df$choose==1,]
df$class<-factor(df$class,levels=c('Lineage_differentiation','Immune','Signaling_transduction','Protein_activities',
                                  'Metabolism','TRANSLATION&Translation'))
table(df$class)
df<-df[order(df$class),]
colourlist<-c(rep("#295734",4),
              rep("#C60000",5),
              rep("#00C8C8",3),
              rep("#601694",4),
              rep("#805C20",4),
              rep("#C46138",4)
             )
options(repr.plot.width=14, repr.plot.height=15)
ggplot(df,aes(x=SIZE,y=NAME))+
    geom_point(aes(size=FDR.q.val,colour=NES))+
    scale_y_discrete(limits = rev(df$NAME), labels = rev(df$NAME))+
    scale_x_continuous(limits =c(0,350) , breaks=c(0,70,140,210,280,350))+
    #geom_vline(xintercept=" ",linetype="solid")+
    scale_colour_gradientn(colours = c("#f29438","white","#4bc1ee"))+
    scale_size(range=c(8,4),limits=c(0,0.2),breaks=c(0,0.01,0.05,0.1,0.2))+
    theme(axis.text.x = element_text(size=20,angle = 60,hjust = 1),
          axis.text.y = element_text(size = 20,colour =rev(colourlist )),
          axis.title.y=element_blank(),
          axis.title.x=element_text(size=20),
          panel.background = element_rect(color = "white",fill="white",linetype="solid"),
          axis.line.y=element_blank(),
          axis.line.x=element_line(color = "black"),
          #panel.grid.major=element_blank(),
          #panel.grid.minor=element_blank(),
          legend.text = element_text(size = 20), 
          legend.title = element_text(size =20),plot.title=element_text(hjust = 0.5,size=25),
          legend.position="bottom",
          #panel.grid.minor = element_line(colour = "red"), 
          panel.grid.major.y = element_line(colour = "black",linetype =3),
          panel.grid.major.x = element_blank())+
    guides(size=guide_legend(reverse = T,ncol=1))+
    xlab("Gene Size")+

NULL
#ggsave("gsea_plot_0210.svg",dpi=600,height = 15,width=14)

#########Fig.4f-4g

##运行环境：/home/hulab/.conda/envs/Gpb_4.2.0

library(ggplot2)
library(enrichplot)
library(data.table)
library(DOSE)
library(RColorBrewer)
library(ggrepel)
gsea_plot <- function(
    gene_set_dir = 'edb/', 
    result_folder = '/data1/02.private/dengyj/analysis/thymus/clustering_10_20/ETP_DN/Enrichment/nov16/12_vs_60.Gsea.1605493817058', 
    set_name = "KEGG_RIBOSOME",
    title = set_name,
    color="#87c444",
    base_size = 11,
    rel_heights=c(1.5, .5, 1),
    subplots = 1:3,
    pvalue_table = FALSE,
    ES_geom="line",
    color.vline="NA",
    color.hline="gray67",
    line_size = 2, 
    label_a = 'A',
    label_b = 'B',
    title_size = 15,
    val_x = c(0.25, 0.25), 
    val_y = c(0.725, 0.675),  
    label_x = c(0.19, 0.88),
    label_y = c(0.315, 0.315),
    label_size = 5, 
    val_size = 5
){
   if(!stringr::str_detect(result_folder, '/$')){
       result_folder <- paste0(result_folder, '/')
   }
    require(ggplot2)
    gene_set_dir <- paste0(result_folder, gene_set_dir)
    gene_set_dir_files <- dir(gene_set_dir)
    gene_set <- paste0(gene_set_dir, gene_set_dir_files[grep('gmt$', gene_set_dir_files)])
    set_name <- toupper(set_name)
    file <- dir(result_folder)
    report_file <- file[grep('(^gsea_report){1}(.*)tsv$', file)]
    result_1 <- data.frame(data.table::fread(paste0(result_folder, report_file[1])), row.names = 'NAME')
    result_2 <- data.frame(data.table::fread(paste0(result_folder, report_file[2])), row.names = 'NAME')
    result <- rbind(result_1, result_2)
    line <- strsplit(readLines(gene_set), '\t')
    location <- which(sapply(line, function(x) x[1]) == set_name)
    gene <- line[[location]][3:length(line[[location]])]
    rank_dir <- file[grep('^ranked_gene_list', file)]
    rank_list <- data.table::fread(paste0(result_folder, rank_dir))
    score <- rank_list$SCORE
    names(score) <- as.character(as.matrix(rank_list$NAME))

    ES_geom <- match.arg(ES_geom, c("line", "dot"))

    geneList <- position <- NULL ## to satisfy codetool

    #gsinfo函数处理
    gseaScores <- getFromNamespace("gseaScores", "DOSE")
    gsdata <- gseaScores(geneList=score, geneSet=gene, exponent=1, fortify=TRUE)
    gsdata$ymin=0
    gsdata$ymax=0
    pos <- gsdata$position == 1
    h <- diff(range(gsdata$runningScore))/20
    gsdata$ymin[pos] <- -h
    gsdata$ymax[pos] <- h
    gsdata$geneList <- score
    gsdata$Description <- "extracellular matrix"

    p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
            theme_classic(base_size) +
            theme(panel.grid.major = element_line(colour = "grey92"),
                panel.grid.minor = element_line(colour = "grey92"),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank()) +
            scale_x_continuous(expand=c(0,0))

        if (ES_geom == "line") {
            es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description), size=line_size)
        } else {
            es_layer <- geom_point(aes_(y = ~runningScore, color= ~Description), size=1, data = subset(gsdata, position == 1))
        }

        p.res <- p + es_layer +
            theme(legend.position = c(.8, .8), legend.title = element_blank(),
                legend.background = element_rect(fill = "transparent"))

        p.res <- p.res + ylab("Running Enrichment Score") +
            theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.line.x=element_blank(),
                plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))

        i <- 0
        for (term in unique(gsdata$Description)) {
            idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
            gsdata[idx, "ymin"] <- i
            gsdata[idx, "ymax"] <- i + 1
            i <- i + 1
        }
        p2 <- ggplot(gsdata, aes_(x = ~x)) +
            geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
            xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
            theme(legend.position = "none",
                plot.margin = margin(t=-.1, b=0,unit="cm"),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                axis.line.x = element_blank()) +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0))


            v <- seq(1, sum(gsdata$position), length.out=9)
            inv <- findInterval(rev(cumsum(gsdata$position)), v)
            if (min(inv) == 0) inv <- inv + 1

            col=c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))

            ymin <- min(p2$data$ymin)
            yy <- max(p2$data$ymax - p2$data$ymin) * .3
            xmin <- which(!duplicated(inv))
            xmax <- xmin + as.numeric(table(inv)[unique(inv)])
            d <- data.frame(ymin = ymin, ymax = yy,
                            xmin = xmin,
                            xmax = xmax,
                            col = col[unique(inv)])
            p2 <- p2 + geom_rect(
                        aes_(xmin=~xmin,
                                xmax=~xmax,
                                ymin=~ymin,
                                ymax=~ymax,
                                fill=~I(col)),
                        data=d,
                        alpha=.9,
                        inherit.aes=FALSE)


        df2 <- p$data #data.frame(x = which(p$data$position == 1))
        df2$y <- p$data$geneList[df2$x]
        p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0), color="grey")
        p.pos <- p.pos + ylab("Ranked list metric") +
            xlab("Rank in Ordered Dataset") +
            theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))

        if (!is.null(title) && !is.na(title) && title != "")
            p.res <- p.res + ggtitle(title) + theme(plot.title=element_text(size=title_size, hjust = 0.5))

            p.res <- p.res + scale_color_manual(values=color)
            if (length(color) == 1) {
                p.res <- p.res + theme(legend.position = "none")
                p2 <- p2 + scale_color_manual(values = "black")
            } else {
                p2 <- p2 + scale_color_manual(values = color)
            }

    enrichmentScore <- if(abs(max(gsdata$runningScore)) >= abs(min(gsdata$runningScore))) {max(gsdata$runningScore)} else {min(gsdata$runningScore)}
    es.df <- data.frame(es = which.min(abs(p$data$runningScore - enrichmentScore)))
    p.res <- p.res + geom_vline(data = es.df, aes_(xintercept = ~es), colour = color.vline, linetype = "dashed") + geom_hline(aes(yintercept=0), colour=color.hline, linetype="dashed")

        plotlist <- list(p.res, p2, p.pos)[subplots]
        n <- length(plotlist)
        plotlist[[n]] <- plotlist[[n]] +
            theme(axis.line.x = element_line(),
                axis.ticks.x=element_line(),
                axis.text.x = element_text())

        if (length(subplots) == 1)
            return(plotlist[[1]] + theme(plot.margin=margin(t=.2, r = .2, b=.2, l=.2, unit="cm")))


        if (length(rel_heights) > length(subplots))
            rel_heights <- rel_heights[subplots]
    pic <- plot_grid(plotlist = plotlist, ncol=1, align="v", rel_heights=rel_heights)
    val_to_plot <- data.frame(x=val_x,y=val_y,
                               col=c("black","black"),
                               text=c(paste0("NES:", round(result[set_name, "NES"],digits=3)), 
                               paste0("italic('P')~val:", round(result[set_name, "NOM.p.val"],digits=3)) 
                               ))
    label_to_plot <- data.frame(x = label_x, y = label_y, col = c('red', 'blue'), text = c(label_a, label_b))                           
    pic <- pic + geom_text(data=val_to_plot,mapping = aes(x, y, color = col, label = text),parse = T, 
                           size = val_size)+
           theme(legend.position = "none") + 
           geom_text_repel(parse = TRUE) +
           geom_text(data = label_to_plot, mapping = aes(x,y, color = col, label = text), size = label_size) + 
           scale_colour_manual(values = c('red' = 'red', 'blue' = 'blue', 'black' = "black"))
    return(pic)
}
library(cowplot)
                         


options(repr.plot.width=6, repr.plot.height=6)
result_folder = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/SEACR005/H4K77ac_only_Y/GO.Gsea.1636628869166'
gsea_plot(    
    result_folder = result_folder, 
    set_name = "GOBP_T_CELL_DIFFERENTIATION", label_a = 'Young', label_b = 'Aged', 
label_y = c(0.315, 0.1),val_x = c(0.3, 0.3))
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/20241211-GSEA-As/plot/GOMF_ATP_HYDROLYSIS_ACTIVITY.pdf",height = 8,width=8,dpi=300)


options(repr.plot.width=6, repr.plot.height=6)
result_folder = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/SEACR005/H4K77ac_only_Y/ReactomePathways.Gsea.1636642463098'
gsea_plot(    
    result_folder = result_folder, 
    set_name = "SIGNALING BY NOTCH", label_a = 'Young', label_b = 'Aged', 
label_y = c(0.315, 0.1),val_x = c(0.3, 0.3))
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/20241211-GSEA-As/plot/GOMF_ATP_HYDROLYSIS_ACTIVITY.pdf",height = 8,width=8,dpi=300)


options(repr.plot.width=6, repr.plot.height=6)
result_folder = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/SEACR005/H4K77ac_only_Y/ReactomePathways.Gsea.1636642463098'
gsea_plot(    
    result_folder = result_folder, 
    set_name = "SIGNALING BY THE B CELL RECEPTOR (BCR)", label_a = 'Young', label_b = 'Aged', 
label_y = c(0.315, 0.1),val_x = c(0.3, 0.3))
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/20241211-GSEA-As/plot/GOMF_ATP_HYDROLYSIS_ACTIVITY.pdf",height = 8,width=8,dpi=300)


options(repr.plot.width=6, repr.plot.height=6)
result_folder = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/SEACR005/H4K77ac_only_O/KEGG.Gsea.1636691512294'
gsea_plot(    
    result_folder = result_folder, 
    set_name = "KEGG_RIBOSOME", label_a = 'Young', label_b = 'Aged', 
label_y = c(0.315, 0.1),val_x = c(0.3, 0.3))
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/20241211-GSEA-As/plot/GOMF_ATP_HYDROLYSIS_ACTIVITY.pdf",height = 8,width=8,dpi=300)


options(repr.plot.width=6, repr.plot.height=6)
result_folder = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/SEACR005/H4K77ac_only_O/ReactomePathways.Gsea.1636698437895'
gsea_plot(    
    result_folder = result_folder, 
    set_name = "CAP-DEPENDENT TRANSLATION INITIATION", label_a = 'Young', label_b = 'Aged', 
label_y = c(0.315, 0.1),val_x = c(0.3, 0.3))
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/20241211-GSEA-As/plot/GOMF_ATP_HYDROLYSIS_ACTIVITY.pdf",height = 8,width=8,dpi=300)


options(repr.plot.width=6, repr.plot.height=6)
result_folder = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/SEACR005/H4K77ac_only_O/ReactomePathways.Gsea.1636698437895'
gsea_plot(    
    result_folder = result_folder, 
    set_name = "EUKARYOTIC TRANSLATION TERMINATION", label_a = 'Young', label_b = 'Aged', 
label_y = c(0.315, 0.1),val_x = c(0.3, 0.3))
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/20241211-GSEA-As/plot/GOMF_ATP_HYDROLYSIS_ACTIVITY.pdf",height = 8,width=8,dpi=300)


#######Fig.4h
###需要在本地使用IGV软件进行可视化 基于mm10
##IGV输入文件位置：
#young:/data4/cut_TAG/H4K77ac/bedgraph_callpeakSEACR/newcall/bgfile/You.bedGraph
###/data4/cut_TAG/H4K77ac/bedgraph_callpeakSEACR/newcall/SEACR_005/You_005.rmbl.bed
#old:/data4/cut_TAG/H4K77ac/bedgraph_callpeakSEACR/newcall/bgfile/Old.bedGraph
###/data4/cut_TAG/H4K77ac/bedgraph_callpeakSEACR/newcall/SEACR_005/Old_005.rmbl.bed

###Fig.4i
######callpeaks 方法：MACS2 默认参数 only peaks
library(ggseqlogo)
library(data.table)
library(ggplot2)
options(repr.plot.width=6, repr.plot.height=2)
for (i in 1:5){
    motif_file<-paste0("/data4/cut_TAG/H4K77ac/use/old/bed/motif/O/knownResults/known",i,".motif")
    data<-fread(motif_file,fill=TRUE,col.names = c("A","C","G","T","V","M"))
    x<-data[,1:4]
    d<-t(as.matrix(x))
    p<-ggseqlogo(d,method="prob")+theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y=element_blank())
    print(p)
    #savefile=paste0("/data1/02.private/guanpb/data/Trim68/cut_tag/H4K77_1st/new analysis/Motif_plot/only_O",i,".svg")
    #ggsave(savefile,width=6,height=2,dpi=600)
}
options(repr.plot.width=6, repr.plot.height=2)
for (i in 1:5){
    motif_file<-paste0("/data4/cut_TAG/H4K77ac/use/old/bed/motif/Y/knownResults/known",i,".motif")
    data<-fread(motif_file,fill=TRUE,col.names = c("A","C","G","T","V","M"))
    x<-data[,1:4]
    d<-t(as.matrix(x))
    p<-ggseqlogo(d,method="prob")+theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y=element_blank())
    print(p)
    #savefile=paste0("/data1/02.private/guanpb/data/Trim68/cut_tag/H4K77_1st/new analysis/Motif_plot/only_Y",i,".svg")
    #ggsave(savefile,width=6,height=2,dpi=600)
}

############Fig.4j是魏霓师姐用AI直接绘制的

###Fig.4k
###运行环境 R_base 3.6.1 
#####callpeaks 方法：MACS2 ATAC-seq对应参数 代码中体现了
###上游分析位置：《1.4 联合ATAC进行DiffBind及下游分析文件获取》
library(DiffBind)
a<-readRDS("/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_seq/GSE162551_dba.rds")
hmap <- colorRampPalette(c("#0080C0","#FF5353"))(n=10)
p<-dba.plotHeatmap(a, contrast=1,
                                           correlations=FALSE, th=1,
                                           scale="row", colScheme =hmap)

##################Fig.4l ---Lymphoid
###上游分析位置：《1.4 联合ATAC进行DiffBind及下游分析文件获取》
####运行环境 Gpb_R42

library(ggplot2)
library(ggpubr)    
matrix<-read.csv("/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/heatmap/mat/Lymphoid_cell_3Kb.mat_x.tab.csv")
matrix<-matrix[,-which(apply(matrix,2,function(x) all(is.na(x))))]
rownames(matrix)<-matrix[,1]


matrix<-matrix[,-1]
matrix<-matrix[,-1]

T_matrix<-as.data.frame(t(matrix))
T_matrix$xlab<-rownames(T_matrix)

mean_Y<-apply(T_matrix[,c("sorted_SRR13186609_uM","sorted_SRR13186610_uM","sorted_SRR13186611_uM","sorted_SRR13186612_uM")],1,sum)

T_matrix$Y.sorted.markdup<-mean_Y

mean_O<-apply(T_matrix[,c("sorted_SRR13186627_uM","sorted_SRR13186628_uM","sorted_SRR13186629_uM","sorted_SRR13186630_uM")],1,sum)

T_matrix$O.sorted.markdup<-mean_O




df_Y<-T_matrix[,c("Y.sorted.markdup","xlab"),drop=F]
df_O<-T_matrix[,c("O.sorted.markdup","xlab"),drop=F]

df_Y$genotype<-"Young"
df_O$genotype<-"Old"

rownames(df_Y)<-NULL
rownames(df_O)<-NULL
colnames(df_Y)<-c("Values","pos","genotype")
colnames(df_O)<-c("Values","pos","genotype")
                             
                             
boxplot(df_Y[,1])
m0<-boxplot.stats(df_Y[,1])$out
df_Y_FT<-df_Y[!df_Y$Values%in%m0,]
boxplot(df_O[,1])
m<-boxplot.stats(df_O[,1])$out
df_O_FT<-df_O[!df_O$Values%in%m,]
#df_c<-rbind(df_Y_FT,df_O_FT)                             
                        
df_c<-rbind(df_Y_FT,df_O_FT)

df_c$genotype<-factor(df_c$genotype,levels=c("Young","Old"))


kruskal.test(Values~genotype, df_c)
                         
df_lym<-df_c
df_lym$Values_nor_3<-(df_lym$Values-min(df_lym$Values))/(max(df_lym$Values)-min(df_lym$Values))
                             
options(repr.plot.width=8, repr.plot.height=7)
ggplot(df_lym, aes(x=genotype,y = Values_nor_3, color=genotype))+
        geom_boxplot(outlier.colour = NA)+
        scale_color_manual(values = c("#33CCFF","#FF9933"))+
        theme(
        # 移除网格
        panel.background = element_rect(color = "white",fill="white",linetype="solid"),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size = 20),
            axis.title=element_blank(),
        #panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
            
          axis.line=element_line(color = "black"),
            legend.text = element_text(size = 10), 
          legend.title = element_text(size =10)
        )+
#geom_point(pch=17)+
geom_jitter(pch=17,size=3,width = 0.2)  +
geom_signif(comparisons = list(c("Young", "Old")),
              test="wilcox.test",map_signif_level=TRUE,textsize=8)
NULL
#ggsave("H4K77ac_ATAC_common_peaks_in_ATAC_Lymphoid_cell_TSS_3Kb_reads_boxplot.svg", height = 7,width=8,dpi=300)

##################Fig.4l ---Myeloid
###上游分析位置：《1.4 联合ATAC进行DiffBind及下游分析文件获取》
####运行环境 Gpb_R42

library(ggplot2)
library(ggpubr)    
matrix<-read.csv("/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/heatmap/mat/Myeloid_cell_3Kb.mat_x.tab.csv")
matrix<-matrix[,-which(apply(matrix,2,function(x) all(is.na(x))))]
rownames(matrix)<-matrix[,1]


matrix<-matrix[,-1]
matrix<-matrix[,-1]

T_matrix<-as.data.frame(t(matrix))
T_matrix$xlab<-rownames(T_matrix)

mean_Y<-apply(T_matrix[,c("sorted_SRR13186609_uM","sorted_SRR13186610_uM","sorted_SRR13186611_uM","sorted_SRR13186612_uM")],1,sum)

T_matrix$Y.sorted.markdup<-mean_Y

mean_O<-apply(T_matrix[,c("sorted_SRR13186627_uM","sorted_SRR13186628_uM","sorted_SRR13186629_uM","sorted_SRR13186630_uM")],1,sum)

T_matrix$O.sorted.markdup<-mean_O




df_Y<-T_matrix[,c("Y.sorted.markdup","xlab"),drop=F]
df_O<-T_matrix[,c("O.sorted.markdup","xlab"),drop=F]

df_Y$genotype<-"Young"
df_O$genotype<-"Old"

rownames(df_Y)<-NULL
rownames(df_O)<-NULL
colnames(df_Y)<-c("Values","pos","genotype")
colnames(df_O)<-c("Values","pos","genotype")
df_c<-rbind(df_Y,df_O)

df_c$genotype<-factor(df_c$genotype,levels=c("Young","Old"))


kruskal.test(Values~genotype, df_c)
df_Mye<-df_c
df_Mye$Values_nor_3<-(df_Mye$Values-min(df_Mye$Values))/(max(df_Mye$Values)-min(df_Mye$Values))
options(repr.plot.width=8, repr.plot.height=7)
ggplot(df_Mye, aes(x=genotype,y = Values_nor_3, color=genotype))+
        geom_boxplot(outlier.colour = NA)+
        scale_colour_manual(values = c("#33CCFF","#FF9933"))+
        theme(
        # 移除网格
        panel.background = element_rect(color = "white",fill="white",linetype="solid"),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size = 20),
            axis.title=element_blank(),
        #panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
            
          axis.line=element_line(color = "black"),
            legend.text = element_text(size = 10), 
          legend.title = element_text(size =10)
        )+
#geom_point(pch=17)+
geom_jitter(pch=21,size=2,width = 0.2,aes(fill=genotype,color=genotype))  +
scale_fill_manual(values = c("#33CCFF","#FF9933"))+
geom_signif(comparisons = list(c("Young", "Old")),
              test="wilcox.test",map_signif_level=TRUE,textsize=8)+
#scale_y_continuous(limits =c(0,8) , breaks=c(0:8))+
NULL
#ggsave("H4K77ac_ATAC_common_peaks_in_ATAC_Myeloid_cell_TSS_3Kb_reads_boxplot.svg", height = 7,width=8,dpi=300)                             

####Fig.4m 基于禹君师兄计算出的GSEA结果进行GSEA 可视化
###DiffBind之后PreRankGSEA 没限制是否位于promoter
####Supplementary Table 8 CUT&Tag与ATAC-seq的 commom peaks 注释文件list
####Supplementary Table 9a 9b 9c PreRankedGSEA 输出结果
library(ggplot2)
library(enrichplot)
library(data.table)
library(DOSE)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
gsea_plot <- function(
    gene_set_dir = 'edb/', 
    result_folder = '/data1/02.private/dengyj/analysis/thymus/clustering_10_20/ETP_DN/Enrichment/nov16/12_vs_60.Gsea.1605493817058', 
    set_name = "KEGG_RIBOSOME",
    title = set_name,
    color="#87c444",
    base_size = 11,
    rel_heights=c(1.5, .5, 1),
    subplots = 1:3,
    pvalue_table = FALSE,
    ES_geom="line",
    color.vline="NA",
    color.hline="gray67",
    line_size = 2, 
    label_a = 'A',
    label_b = 'B',
    title_size = 15,
    val_x = c(0.25, 0.25), 
    val_y = c(0.725, 0.675),  
    label_x = c(0.19, 0.88),
    label_y = c(0.315, 0.315),
    label_size = 5, 
    val_size = 5
){
   if(!stringr::str_detect(result_folder, '/$')){
       result_folder <- paste0(result_folder, '/')
   }
    require(ggplot2)
    gene_set_dir <- paste0(result_folder, gene_set_dir)
    gene_set_dir_files <- dir(gene_set_dir)
    gene_set <- paste0(gene_set_dir, gene_set_dir_files[grep('gmt$', gene_set_dir_files)])
    set_name <- toupper(set_name)
    file <- dir(result_folder)
    report_file <- file[grep('(^gsea_report){1}(.*)tsv$', file)]
    result_1 <- data.frame(data.table::fread(paste0(result_folder, report_file[1])), row.names = 'NAME')
    result_2 <- data.frame(data.table::fread(paste0(result_folder, report_file[2])), row.names = 'NAME')
    result <- rbind(result_1, result_2)
    line <- strsplit(readLines(gene_set), '\t')
    location <- which(sapply(line, function(x) x[1]) == set_name)
    gene <- line[[location]][3:length(line[[location]])]
    rank_dir <- file[grep('^ranked_gene_list', file)]
    rank_list <- data.table::fread(paste0(result_folder, rank_dir))
    score <- rank_list$SCORE
    names(score) <- as.character(as.matrix(rank_list$NAME))

    ES_geom <- match.arg(ES_geom, c("line", "dot"))

    geneList <- position <- NULL ## to satisfy codetool

    #gsinfo函数处理
    gseaScores <- getFromNamespace("gseaScores", "DOSE")
    gsdata <- gseaScores(geneList=score, geneSet=gene, exponent=1, fortify=TRUE)
    gsdata$ymin=0
    gsdata$ymax=0
    pos <- gsdata$position == 1
    h <- diff(range(gsdata$runningScore))/20
    gsdata$ymin[pos] <- -h
    gsdata$ymax[pos] <- h
    gsdata$geneList <- score
    gsdata$Description <- "extracellular matrix"

    p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
            theme_classic(base_size) +
            theme(panel.grid.major = element_line(colour = "grey92"),
                panel.grid.minor = element_line(colour = "grey92"),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank()) +
            scale_x_continuous(expand=c(0,0))

        if (ES_geom == "line") {
            es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description), size=line_size)
        } else {
            es_layer <- geom_point(aes_(y = ~runningScore, color= ~Description), size=1, data = subset(gsdata, position == 1))
        }

        p.res <- p + es_layer +
            theme(legend.position = c(.8, .8), legend.title = element_blank(),
                legend.background = element_rect(fill = "transparent"))

        p.res <- p.res + ylab("Running Enrichment Score") +
            theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.line.x=element_blank(),
                plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))

        i <- 0
        for (term in unique(gsdata$Description)) {
            idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
            gsdata[idx, "ymin"] <- i
            gsdata[idx, "ymax"] <- i + 1
            i <- i + 1
        }
        p2 <- ggplot(gsdata, aes_(x = ~x)) +
            geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
            xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
            theme(legend.position = "none",
                plot.margin = margin(t=-.1, b=0,unit="cm"),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                axis.line.x = element_blank()) +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0))


            v <- seq(1, sum(gsdata$position), length.out=9)
            inv <- findInterval(rev(cumsum(gsdata$position)), v)
            if (min(inv) == 0) inv <- inv + 1

            col=c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))

            ymin <- min(p2$data$ymin)
            yy <- max(p2$data$ymax - p2$data$ymin) * .3
            xmin <- which(!duplicated(inv))
            xmax <- xmin + as.numeric(table(inv)[unique(inv)])
            d <- data.frame(ymin = ymin, ymax = yy,
                            xmin = xmin,
                            xmax = xmax,
                            col = col[unique(inv)])
            p2 <- p2 + geom_rect(
                        aes_(xmin=~xmin,
                                xmax=~xmax,
                                ymin=~ymin,
                                ymax=~ymax,
                                fill=~I(col)),
                        data=d,
                        alpha=.9,
                        inherit.aes=FALSE)


        df2 <- p$data #data.frame(x = which(p$data$position == 1))
        df2$y <- p$data$geneList[df2$x]
        p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0), color="grey")
        p.pos <- p.pos + ylab("Ranked list metric") +
            xlab("Rank in Ordered Dataset") +
            theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))

        if (!is.null(title) && !is.na(title) && title != "")
            p.res <- p.res + ggtitle(title) + theme(plot.title=element_text(size=title_size, hjust = 0.5))

            p.res <- p.res + scale_color_manual(values=color)
            if (length(color) == 1) {
                p.res <- p.res + theme(legend.position = "none")
                p2 <- p2 + scale_color_manual(values = "black")
            } else {
                p2 <- p2 + scale_color_manual(values = color)
            }

    enrichmentScore <- if(abs(max(gsdata$runningScore)) >= abs(min(gsdata$runningScore))) {max(gsdata$runningScore)} else {min(gsdata$runningScore)}
    es.df <- data.frame(es = which.min(abs(p$data$runningScore - enrichmentScore)))
    p.res <- p.res + geom_vline(data = es.df, aes_(xintercept = ~es), colour = color.vline, linetype = "dashed") + geom_hline(aes(yintercept=0), colour=color.hline, linetype="dashed")

        plotlist <- list(p.res, p2, p.pos)[subplots]
        n <- length(plotlist)
        plotlist[[n]] <- plotlist[[n]] +
            theme(axis.line.x = element_line(),
                axis.ticks.x=element_line(),
                axis.text.x = element_text())

        if (length(subplots) == 1)
            return(plotlist[[1]] + theme(plot.margin=margin(t=.2, r = .2, b=.2, l=.2, unit="cm")))


        if (length(rel_heights) > length(subplots))
            rel_heights <- rel_heights[subplots]
    pic <- plot_grid(plotlist = plotlist, ncol=1, align="v", rel_heights=rel_heights)
    val_to_plot <- data.frame(x=val_x,y=val_y,
                               col=c("black","black"),
                               text=c(paste0("NES:", round(result[set_name, "NES"],digits=3)), 
                               paste0("italic('P')~val:", round(result[set_name, "NOM.p.val"],digits=3)) 
                               ))
    label_to_plot <- data.frame(x = label_x, y = label_y, col = c('red', 'blue'), text = c(label_a, label_b))                           
    pic <- pic + geom_text(data=val_to_plot,mapping = aes(x, y, color = col, label = text),parse = T, 
                           size = val_size)+
           theme(legend.position = "none") + 
           geom_text_repel(parse = TRUE) +
           geom_text(data = label_to_plot, mapping = aes(x,y, color = col, label = text), size = label_size) + 
           scale_colour_manual(values = c('red' = 'red', 'blue' = 'blue', 'black' = "black"))
    return(pic)
}
result_folder = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NewGSEA/combian_ATAC/GO.GseaPreranked.1674904593108'
gsea_plot(    
    result_folder = result_folder, 
    set_name = "GO_REGULATION_OF_NOTCH_SIGNALING_PATHWAY", label_a = 'pos', label_b = 'neg', 
label_y = c(0.315, 0.1))
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/GSEA_Plot/ATAC_GSEA/GO_REGULATION_OF_NOTCH_SIGNALING_PATHWAY.svg",height = 6,width=6,dpi=300)
result_folder = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NewGSEA/combian_ATAC/GO.GseaPreranked.1674904593108'
gsea_plot(    
    result_folder = result_folder, 
    set_name = "GO_T_CELL_DIFFERENTIATION", label_a = 'pos', label_b = 'neg', 
label_y = c(0.315, 0.1))
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/GSEA_Plot/ATAC_GSEA/GO_T_CELL_DIFFERENTIATION.svg",height = 6,width=6,dpi=300)
result_folder = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NewGSEA/combian_ATAC/GO.GseaPreranked.1674904593108'
gsea_plot(    
    result_folder = result_folder, 
    set_name = "GO_T_CELL_ACTIVATION", label_a = 'pos', label_b = 'neg', 
label_y = c(0.315, 0.1))
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/GSEA_Plot/ATAC_GSEA/GO_T_CELL_ACTIVATION.svg",height = 6,width=6,dpi=300)
result_folder = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NewGSEA/combian_ATAC/Reactiome.GseaPreranked.1674904688870'
gsea_plot(    
    result_folder = result_folder, 
    set_name = "SIGNALING BY THE B CELL RECEPTOR (BCR)", label_a = 'pos', label_b = 'neg', 
label_y = c(0.315, 0.1))
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/GSEA_Plot/ATAC_GSEA/SIGNALING BY THE B CELL RECEPTOR (BCR).svg",height = 6,width=6,dpi=300)
                          

###Extended Data Fig.8b  该图分为2部分 一部分为上游统计文件计算 第二部分为可视化结果

####第一部分 上游统计文件计算 
##bw文件的生成以及淋系和髓系基因集的生成在《1.1上游处理过程》 代码中  

#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！ 
##运行环境 base
##运行位置：/data4/cut_TAG/H4K77ac/use/old/heatmap/old

#######Young
computeMatrix reference-point \
	--referencePoint TSS \
	-b 3000 -a 3000  -bs 3 \
	-R /data4/cut_TAG/H4K77ac/use/Lymphoid_Cell_genes.bed /data4/cut_TAG/H4K77ac/use/Myeloid_Cell_genes.bed  \
	-S /data4/cut_TAG/H4K77ac/bw/Y.sorted.markdup.bw \
	-o heatmap/mat/H4K77ac_Young_LvsM_TSS_3Kb.mat.gz -p 8 --missingDataAsZero \
	--outFileSortedRegions heatmap/mat/regions_H4K77ac_Young_LvsM_TSS_3Kb.bed


	plotProfile -m heatmap/mat/H4K77ac_Young_LvsM_TSS_3Kb.mat.gz  \
	-out heatmap/plot/H4K77ac_Young_LvsM_TSS_3Kb.png \
	--dpi 600 --plotWidth 8 \
	--outFileNameData heatmap/mat/H4K77ac_Young_LvsM_TSS_3Kb.mat_x.tab 

####################Aged
computeMatrix reference-point \
	--referencePoint TSS \
	-b 3000 -a 3000  -bs 3 \
	-R /data4/cut_TAG/H4K77ac/use/Lymphoid_Cell_genes.bed /data4/cut_TAG/H4K77ac/use/Myeloid_Cell_genes.bed  \
	-S /data4/cut_TAG/H4K77ac/bw/O.sorted.markdup.bw \
	-o heatmap/mat/H4K77ac_Old_LvsM_TSS_3Kb.mat.gz -p 8 --missingDataAsZero \
	--outFileSortedRegions heatmap/mat/regions_H4K77ac_Old_LvsM_TSS_3Kb.bed


	plotProfile -m heatmap/mat/H4K77ac_Old_LvsM_TSS_3Kb.mat.gz  \
	-out heatmap/plot/H4K77ac_Old_LvsM_TSS_3Kb.png \
	--dpi 600 --plotWidth 8 \
	--outFileNameData heatmap/mat/H4K77ac_OLd_LvsM_TSS_3Kb.mat_x.tab 


###Extended Data Fig.8b  该图分为2部分 一部分为上游统计文件计算 第二部分为可视化结果

####第二部分 可视化结果
###输入文件位置:/data4/cut_TAG/H4K77ac/use/old/heatmap/old/mat


##运行环境 Gpb_R42

library(ggplot2)
matrix<-read.csv("/data4/cut_TAG/H4K77ac/use/old/heatmap/old/mat/H4K77ac_Young_LvsM_TSS_3Kb.mat_x.tab.csv")
#matrix<-matrix[,-which(apply(matrix,2,function(x) all(is.na(x))))]
rownames(matrix)<-matrix[,1]

matrix<-matrix[,-1]
matrix<-matrix[,-1]
T_matrix<-as.data.frame(t(matrix))
T_matrix$xlab<-rownames(T_matrix)
df_L<-T_matrix[,c("Lymphoid_Cell_genes.bed","xlab"),drop=F]
df_M<-T_matrix[,c("Myeloid_Cell_genes.bed","xlab"),drop=F]

df_L$genotype<-"Lymphoid"
df_M$genotype<-"Myeloid"

rownames(df_L)<-NULL
rownames(df_M)<-NULL
colnames(df_L)<-c("Values","pos","genotype")
colnames(df_M)<-c("Values","pos","genotype")

boxplot(df_L[,1])
m0<-boxplot.stats(df_L[,1])$out
df_L_FT<-df_L[!df_L$Values%in%m0,]
boxplot(df_M[,1])
m<-boxplot.stats(df_M[,1])$out
df_M_FT<-df_M[!df_M$Values%in%m,]
df_c<-rbind(df_L_FT,df_M_FT)



df_c$genotype<-factor(df_c$genotype,levels=c("Lymphoid","Myeloid"))

kruskal.test(Values~genotype, df_c)
df_c$x_Values<-df_c$Values-(median(df_c$Values[df_c$genotype=="Lymphoid"])-1)
ggplot(df_c, aes(pos, group=genotype, fill=genotype))+geom_boxplot(aes(y = x_Values))+
                             scale_fill_manual(values = c("Lymphoid"="#8CA9C7","Myeloid"="#E5C4DF"))
#ggsave("H4K77ac_young_LvsM_signal_violinplot_fill.svg",width=6,height=6)   



matrix<-read.csv("/data4/cut_TAG/H4K77ac/use/old/heatmap/old/mat/H4K77ac_OLd_LvsM_TSS_3Kb.mat_x.tab.csv")
#matrix<-matrix[,-which(apply(matrix,2,function(x) all(is.na(x))))]
rownames(matrix)<-matrix[,1]

matrix<-matrix[,-1]
matrix<-matrix[,-1]
T_matrix<-as.data.frame(t(matrix))
T_matrix$xlab<-rownames(T_matrix)
df_L<-T_matrix[,c("Lymphoid_Cell_genes.bed","xlab"),drop=F]
df_M<-T_matrix[,c("Myeloid_Cell_genes.bed","xlab"),drop=F]

df_L$genotype<-"Lymphoid"
df_M$genotype<-"Myeloid"

rownames(df_L)<-NULL
rownames(df_M)<-NULL
colnames(df_L)<-c("Values","pos","genotype")
colnames(df_M)<-c("Values","pos","genotype")

boxplot(df_L[,1])
m0<-boxplot.stats(df_L[,1])$out
df_L_FT<-df_L[!df_L$Values%in%m0,]
boxplot(df_M[,1])
m<-boxplot.stats(df_M[,1])$out
df_M_FT<-df_M[!df_M$Values%in%m,]
df_c<-rbind(df_L_FT,df_M_FT)



df_c$genotype<-factor(df_c$genotype,levels=c("Lymphoid","Myeloid"))

kruskal.test(Values~genotype, df_c)
df_c$x_Values<-df_c$Values-(median(df_c$Values[df_c$genotype=="Lymphoid"])-1)
ggplot(df_c, aes(pos, group=genotype, fill=genotype))+geom_boxplot(aes(y = x_Values))+
                             scale_fill_manual(values = c("Lymphoid"="#8CA9C7","Myeloid"="#E5C4DF"))
#ggsave("H4K77ac_old_LvsM_signal_violinplot_fill.svg",width=6,height=6)    

###########涉及Extended Data Fig.9o-9q
##分析过程以及可视化代码在《1.5 H4K77ac Young联合Pubilc H3K27ac Young》
##############################

#####Fig.5j ------Treated in Young HSC
########Supplementary Table 11a 11b peaks注释list
####Supplementary Table 10 DMSO vs RGFP966 HSC 蛋白组原始结果
####GSEA相关的运行代码在《2.2 联合蛋白质组PreRanked GSEA分析》
##callpeaks 方法：SEACR Top5% only peaks TSS±3Kb
###only_peaks 注释基因与蛋白质组overlap的基因 之后基因蛋白质的FC 生成。rnk 之后PreRankedGSEA 
###GSEA输入文件位置：/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/RGFP966
###GSEA输出结果文件位置：/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/RGFP966/GSEA/Protein_only_RGFP966
###运行环境：Gpb_R42

library(openxlsx)
pathway<-read.csv("/data1/02.private/guanpb/data/Trim68/cut_tag/RGFP966/GSEA/Highlight_pathway_2.csv")
pathway<-pathway[order(pathway$SIZE),]
options(repr.plot.width=14.5, repr.plot.height=10)

ggplot(pathway,aes(x=-NES,y=NAME))+
    #geom_segment( aes(x=0, xend=NES, y=NAME,yend=NAME)) +
    geom_point( size=8,aes(color=NOM.p.val,fill=NOM.p.val), shape=21, stroke=2)+
scale_y_discrete(limits = (pathway$NAME), labels = (pathway$NAME))+
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
scale_fill_gradientn("NOM.p.val",colours = c("#E3728A","white"),limits=c(0,0.1),breaks=c(0,0.05,0.1))+
scale_color_gradientn("NOM.p.val",colours = c("#E3728A","white"),limits=c(0,0.1),breaks=c(0,0.05,0.1))+
#scale_size_continuous(range = c(8,4),guide = FALSE,limits=c(0,0.1),breaks=c(0,0.05,0.1))+
    #scale_x_continuous(limits =c(0,0.06) , breaks=c(0,75,150,225))+
NULL
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/plot/H4K77ac_RGFP966_Aged_CT_DiffBind_pathway_Young_dotplot.svg",width=15,height=10,dpi=600)

#####Extended Data Fig.9n
###运行环境：R_base 3.6.1

library(DiffBind)
H4K77ac<-readRDS("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/H4K77ac_RGFP966_Aged_DiffBind_ht_df.rds")
#svg('/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/plot/H4K77ac_RGFP966_Aged_CT_DiffBind_ht.svg',width=6,height=6)

hmap <- colorRampPalette(c("#363666","#427471", "#FCE525"))(n=10)

p<-dba.plotHeatmap(H4K77ac, contrast=1,
                                           correlations=FALSE, th=1,
                                           scale="row", colScheme =hmap)
print(p)
#dev.off()

#####Fig.7j ------Treated in Aged HSC
####Supplementary Table 11c 11d peaks 注释list
##Supplementary Table 11e 11f PreRankedGSEA 输出结果
####GSEA相关的运行代码在《3.2 DiffBind后通路富集分析》
##callpeaks 方法：SEACR Top5% 
###直接使用DiffBind计算的FC 生成。rnk 之后PreRankedGSEA 
###GSEA输入文件位置：/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/input
###GSEA输出结果文件位置：/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/out/H4K77ac_RGFP966_DiffBind_prerankedGSEA
###运行环境：Gpb_R42

library(openxlsx)
library(Hmisc)
pathway<-read.xlsx("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/pathway_01.xlsx",
                  sheet=4)
name<-c('REGULATION OF T CELL MIGRATION',
'REGULATION OF LYMPHOCYTE MIGRATION',
'POSITIVE REGULATION OF B CELL PROLIFERATION',
'REGULATION OF T CELL MEDIATED CYTOTOXICITY',
'MATURE B CELL DIFFERENTIATION',
'REGULATION OF B CELL RECEPTOR SIGNALING PATHWAY',
'T CELL DIFFERENTIATION INVOLVED IN IMMUNE RESPONSE',
'POSITIVE REGULATION OF CHEMOKINE PRODUCTION',
'REGULATION OF PROTEIN KINASE A SIGNALING',
'HEMATOPOIETIC CELL LINEAGE')
pathway$name<-name
df_pathway<-pathway[order(pathway$SIZE,decreasing = T),]
options(repr.plot.width=17, repr.plot.height=10)

ggplot(df_pathway,aes(x=-NES,y=name))+
    #geom_segment( aes(x=0, xend=NES, y=NAME,yend=NAME)) +
    geom_point( size=8,aes(color=NOM_p_val,fill=NOM_p_val), shape=21, stroke=2)+
scale_y_discrete(limits = rev(df_pathway$name), labels = rev(df_pathway$name))+
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
scale_fill_gradientn("NOM_p_val",colours = c("#D94811","white"),limits=c(0,0.1),breaks=c(0,0.05,0.1))+
scale_color_gradientn("NOM_p_val",colours = c("#D94811","white"),limits=c(0,0.1),breaks=c(0,0.05,0.1))+
#scale_size_continuous(range = c(8,4),guide = FALSE,limits=c(0,0.1),breaks=c(0,0.05,0.1))+
    #scale_x_continuous(limits =c(0,0.06) , breaks=c(0,75,150,225))+
NULL
##ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/plot/H4K77ac_RGFP966_Aged_CT_DiffBind_pathway_Aged_dotplot.svg",width=18,height=10,dpi=600)

#########Fig.5k
##IGV可视化输入文件位置 基于mm10
###DMSO:/data/02.private/guanpb/H4K77ac/CUTTag/bgfile/K77_Ct_0um.bedGraph
##RGFP966:/data/02.private/guanpb/H4K77ac/CUTTag/bgfile/K77_Ct_1um.bedGraph
##文件获取代码在《3.1 上游分析过程》

###Fig.5m 
##魏霓师姐通过AI绘制完成

##Fig.5n
##运行环境：R_base 3.6.3
#####上游分析位置：《4.2 联合Aged CUT&tag 数据 DiffBind分析后GSEA分析》
library(DiffBind)
RGFP966<-readRDS("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/CT_ATAC_combind_DiffBind_ht.rds")
#svg('/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/plot/CT_ATAC_combind_DiffBind_ht.svg',width=6,height=6)

hmap <- colorRampPalette(c("#8DCED6","#DAE2F0","#E44832"))(n=10)

p<-dba.plotHeatmap(RGFP966, contrast=1,
                                           correlations=FALSE, th=1,
                                           scale="row", colScheme =hmap)
print(p)
#dev.off()

##Fig.5o
####Supplementary Table 9d 9e 9f PreRankedGSEA 输出结果
####GSEA相关的运行代码在《4.2 联合Aged CUT&tag 数据 DiffBind分析后GSEA分析》
##callpeaks 方法：CUT&Tag:SEACR Top5% only peaks TSS±3Kb ATAC-seq:MASC2
###基于ATAC-seq的值进行DiffBind计算的FC 生成。rnk 之后PreRankedGSEA 
###GSEA输入文件位置：/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/input
###GSEA输出结果文件位置：/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/out/ATAC_CT_combined_DiffBind_By_CT_prerankedGSEA
###运行环境：Gpb_R42
##运行环境：Gpb_R42
library(openxlsx)
library(ggplot2)
pathway<-read.xlsx("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/pathway_combine_by_CT.xlsx",
                  sheet=4)
pathway<-pathway[order(pathway$SIZE),]
options(repr.plot.width=18, repr.plot.height=10)
ggplot(pathway,aes(x=-NES,y=NAME,color=NOM_p_val,size=SIZE))+
    geom_point(stat='identity')+
scale_y_discrete(limits = (pathway$NAME), labels = (pathway$NAME))+
theme(axis.text.x = element_text(size=20),
          axis.text.y = element_text(size = 36, face = 'bold'),
          axis.title.x = element_text(size=20),
          axis.title.y=element_blank(),
          panel.background = element_rect(color = "black",fill="white",linetype="solid"),
          axis.line=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.text = element_text(size = 20), 
          legend.title = element_text(size =20),plot.title=element_text(hjust = 0.5,size=25))+
scale_color_gradientn(limits=c(0,0.1),breaks=c(0,0.05,0.1),colours = c("#F9D5DD","#E3728A","#703745"))+
scale_size(range = c(6,12))+
#geom_text(aes(label=SIZE),size=10,angle=360,hjust=1)
NULL
#ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/plot/pathway_combine_by_CT_Dotplot.svg",width=18,height=10,dpi=600)



######Extended Data Fig.9o 9n 9q 
###分析过程及可视化具体代码：《1.5 H4K77ac Young联合Pubilc H3K27ac Young.ipynb》

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

######Extended Data Fig.9q 
##Supplementary Table 12a  CUT&Tag 和 ATAC-seq common peaks 注释list
##Supplementary Table 12b PreRankedGSEA 输出结果
###分析过程及可视化具体代码：《1.5 H4K77ac Young联合Pubilc H3K27ac Young.ipynb》
###基于overlap的基因后在转录组中做的GSEA
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


