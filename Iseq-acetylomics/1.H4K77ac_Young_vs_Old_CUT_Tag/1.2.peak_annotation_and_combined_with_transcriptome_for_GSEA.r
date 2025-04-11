##############SEACR callpeak 之后only_Young Only_Old 文件获取 

###callpeaks为上游分析 因此位于1.1上游处理过程代码文件中

#######基于Linux base 环境运行！！！！！！！！！！！！！！！！！！！！ 
#运行位置：/data4/cut_TAG/H4K77ac/bedgraph_callpeakSEACR/newcall/SEACR_005
mkdir annobed
bedtools intersect -a You_005.rmbl.bed -b Old_005.rmbl.bed -wa -wb >m.bed
bedtools merge -i m.bed > annobed/common_peaks_005.bed
##
bedtools intersect -a You_005.rmbl.bed -b annobed/common_peaks_005.bed -v > annobed/only_Y_peaks_005.bed
bedtools intersect -a Old_005.rmbl.bed -b annobed/common_peaks_005.bed -v > annobed/only_O_peaks_005.bed

###########以下代码基于 Gpb_R_base (3.6.1)

library("ChIPseeker")
library("org.Mm.eg.db")
library("clusterProfiler")
library("GenomicFeatures")
library(corrplot)
library(enrichplot)
library(ReactomePA)
library(dplyr)
library("TxDb.Mmusculus.UCSC.mm10.knownGene") 
library(Seurat)

#########是用的only的
Old<-readPeakFile("/data4/cut_TAG/H4K77ac/bedgraph_callpeakSEACR/newcall/SEACR_005/annobed/only_O_peaks_005.bed")
Young<-readPeakFile("/data4/cut_TAG/H4K77ac/bedgraph_callpeakSEACR/newcall/SEACR_005/annobed/only_Y_peaks_005.bed")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

OldpeakAnno <- annotatePeak(Old, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")

YoungpeakAnno <- annotatePeak(Young, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")

plotAnnoPie(OldpeakAnno)
plotAnnoPie(YoungpeakAnno)
plotAnnoPie(CommonpeakAnno)

Y_adata<-as.data.frame(YoungpeakAnno)
O_adata<-as.data.frame(OldpeakAnno)

Y_adata$anno<-""
for (i in 1:length(rownames(Y_adata))){
    a<-Y_adata$annotation[i]
    b<-unlist(strsplit(a,''))
    n<-length(which(b=="("))
    if (n==0){
        Y_adata$anno[i]<-a
    }else{
        n1<-which(b=="(")
        c<-b[1:n1-1]
        Y_adata$anno[i]<-c(paste(c,collapse=''))
    }
}

O_adata$anno<-""
for (i in 1:length(rownames(O_adata))){
    a<-O_adata$annotation[i]
    b<-unlist(strsplit(a,''))
    n<-length(which(b=="("))
    if (n==0){
        O_adata$anno[i]<-a
    }else{
        n1<-which(b=="(")
        c<-b[1:n1-1]
        O_adata$anno[i]<-c(paste(c,collapse=''))
    }
}

Y_gene_Tss<-unique(Y_adata$SYMBOL[Y_adata$anno=='Promoter '])
O_gene_Tss<-unique(O_adata$SYMBOL[O_adata$anno=='Promoter '])
Y_gene<-unique(Y_adata$SYMBOL)
O_gene<-unique(O_adata$SYMBOL)

length(unique(Y_adata$SYMBOL[Y_adata$anno=='Promoter ']))
length(unique(O_adata$SYMBOL[O_adata$anno=='Promoter ']))
length(unique(Y_gene))
length(unique(O_gene))

#######转录组数据来源于禹君师兄

HSC<-readRDS("/data1/02.private/other/11_5/object2.rds")
genelist_GSEA<-unique(c(Y_gene_Tss))

ff_HSC<-HSC[rownames(HSC)%in%genelist_GSEA,]

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
get_gsea_file(ff_HSC, idents_vars = 'time', exp_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/ct_SEACR005_H4K77ac_only_Y_GSEA_matrix.txt', 
              cls_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/ct_SEACR005_H4K77ac_only_Y_GSEA_cls.cls',gseapy = F, add_prefix = F)

length(genelist_GSEA)


HSC

ff_HSC

18729-3416
4062-3416

HSC<-readRDS("/data1/02.private/other/11_5/object2.rds")
genelist_GSEA<-unique(c(O_gene_Tss))

ff_HSC<-HSC[rownames(HSC)%in%genelist_GSEA,]

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
get_gsea_file(ff_HSC, idents_vars = 'time', exp_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/ct_SEACR005_H4K77ac_only_O_GSEA_matrix.txt', 
              cls_name = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/ct_SEACR005_H4K77ac_only_O_GSEA_cls.cls',gseapy = F, add_prefix = F)



length(genelist_GSEA)
ff_HSC
HSC

18729-4297
4972-4297

#####进行GSEA运算

#######基于Linux base 环境运行！！！！！！！！！！！！！！！！！！！！ 

###运行环境：base
##运行位置：/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/SEACR005
############

#!/bin/bash

cluster="H4K77ac_only_Y H4K77ac_only_O"

gmtlist='GO KEGG ReactomePathways'

for name in ${cluster}
do
	mkdir ${name}
	for gmt in ${gmtlist}
	do	
		gsea-cli.sh GSEA -res /data1/02.private/guanpb/data/Trim68/cut_tag/H4K77_1st/new\ analysis/ct_SEACR005_${name}_GSEA_matrix.txt \
		-cls /data1/02.private/guanpb/data/Trim68/cut_tag/H4K77_1st/new\ analysis/ct_SEACR005_${name}_GSEA_cls.cls#2M_versus_28M \
		-gmx /data1/02.private/guanpb/GSEA_GMT/${gmt}.gmt \
		-collapse Collapse -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
		-zip_report false -out ${name} -rpt_label ${gmt}
	done
done

###输出文件位置：
#only Young GSEA:/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/SEACR005/H4K77ac_only_Y
#only Old GSEA:/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/SEACR005/H4K77ac_only_O
