############基于DiffBind 计算差异peaks
###运行环境：R_base 3.6.1

library(DiffBind)

#####DiffBind输入文件的获取代码在 《3.1 上游分析过程》
samples <- read.csv("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/RGFP966_Aged_CT_DiffBind.csv")

samples

H4K77ac <- dba(sampleSheet=samples)
H4K77ac <- dba.count(H4K77ac, bUseSummarizeOverlaps=TRUE)
plot(H4K77ac)
H4K77ac <- dba.contrast(H4K77ac, categories = DBA_TREATMENT,minMembers = 2)
H4K77ac <- dba.analyze(H4K77ac)
hmap <- colorRampPalette(c("green", "black", "red"))(n = 13)
readscores <-dba.plotHeatmap(H4K77ac, contrast=1,correlations=FALSE, scale="row", colScheme = hmap)

####涉及Extended Data Fig.9n
saveRDS(H4K77ac,"/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/H4K77ac_RGFP966_Aged_DiffBind_ht_df.rds")


dba.plotMA(H4K77ac, contrast = 1,bUsePval=T)
dba.plotMA(H4K77ac, contrast = 1,bUsePval=F,th=0.25)
dba.plotVolcano(H4K77ac,,bUsePval=T)
H4K77ac <- dba.report(H4K77ac,  method=DBA_DESEQ2, contrast = 1, th=1)

H4K77ac

####保存差异peak计算结果
saveRDS(H4K77ac,"/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/H4K77ac_CT_DiffBind.rds")

###########对差异peaks进行注释
##运行环境：Gpb_R42

library("ChIPseeker")
library("org.Mm.eg.db")
library("clusterProfiler")
library("GenomicFeatures")
library(corrplot)
library(enrichplot)
library(ReactomePA)
library(dplyr)
library("TxDb.Mmusculus.UCSC.mm10.knownGene")


library(plyr)

H4K77ac_RGF966<-readRDS("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/H4K77ac_CT_DiffBind.rds")

H4K77ac_RGF966<-readRDS("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/H4K77ac_CT_DiffBind.rds")

H4K77ac_RGF966

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
H4K77ac_RGF966_peakAnno <- annotatePeak(H4K77ac_RGF966, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoPie(H4K77ac_RGF966_peakAnno)
H4K77ac_RGF966_adata<-as.data.frame(H4K77ac_RGF966_peakAnno)
H4K77ac_RGF966_adata$anno<-""
for (i in 1:length(rownames(H4K77ac_RGF966_adata))){
    a<-H4K77ac_RGF966_adata$annotation[i]
    b<-unlist(strsplit(a,''))
    n<-length(which(b=="("))
    if (n==0){
        H4K77ac_RGF966_adata$anno[i]<-a
    }else{
        n1<-which(b=="(")
        c<-b[1:n1-1]
        H4K77ac_RGF966_adata$anno[i]<-c(paste(c,collapse=''))
    }
}


H4K77ac_RGF966_adata_ft<-H4K77ac_RGF966_adata[H4K77ac_RGF966_adata$p.value<0.1,]

write.csv(H4K77ac_RGF966_adata_ft,"/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/H4K77ac_RGF966_vs_DMSO_DiffBind_Aged_adata.csv")

filtered_H4K77ac_RGFP966 <- 
ddply(H4K77ac_RGF966_adata_ft, 'SYMBOL', summarise, Fold = mean(Fold))

filtered_H4K77ac_RGFP966 <- 
filtered_H4K77ac_RGFP966[!is.na(filtered_H4K77ac_RGFP966$SYMBOL), ]
filtered_H4K77ac_RGFP966 <- filtered_H4K77ac_RGFP966[
    order(filtered_H4K77ac_RGFP966$Fold, decreasing = T),]

write.table(filtered_H4K77ac_RGFP966[, c("SYMBOL", 'Fold')], 
            file = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/input/H4K77ac_RGFP966_DiffBind.rnk',
            col.names = F, quote = F, sep = '\t', row.names = F )

#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！
##运行环境：base
##运行位置：/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/out

#!/bin/bash
gmtlist='GO KEGG ReactomePathways'
mkdir H4K77ac_RGFP966_DiffBind_prerankedGSEA
for gmt in ${gmtlist}
do
        gsea-cli.sh GSEAPreranked -gmx /data1/02.private/guanpb/GSEA_GMT/${gmt}.gmt \
        -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/input/H4K77ac_RGFP966_DiffBind.rnk' \
        -zip_report false -out H4K77ac_RGFP966_DiffBind_prerankedGSEA -rpt_label ${gmt} \
        -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
        -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 10
done
