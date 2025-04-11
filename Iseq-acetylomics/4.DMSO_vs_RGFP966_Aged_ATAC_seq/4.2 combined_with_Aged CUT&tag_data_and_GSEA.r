########具体思路：
#先找CUT&Tag里的 0uM 和1uM内部保守的peaks（至少在2个样本出现）  之后和ATAC-seq 0，1 分别做交集   
###之后把得到交集的结果 合并 输入到CUT&tag中进行DiffBind

####先找CUT&Tag中保守的Pekas

#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！

##运行位置：/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005
##运行环境：base

#!/bin/bash
mkdir common

bedtools intersect -a 0uM-501-X-701_005.rmbl.bed -b 0uM-501-X-702_005.rmbl.bed -wa -wb >common/c12_0uM.bed 
bedtools intersect -a 0uM-501-X-701_005.rmbl.bed -b 0uM-501-X-703_005.rmbl.bed -wa -wb >common/c13_0uM.bed
bedtools intersect -a 0uM-501-X-702_005.rmbl.bed -b 0uM-501-X-703_005.rmbl.bed -wa -wb >common/c23_0uM.bed


bedtools intersect -a 1uM-502-X-701_005.rmbl.bed -b 1uM-502-X-702_005.rmbl.bed -wa -wb >common/c12_1uM.bed
bedtools intersect -a 1uM-502-X-701_005.rmbl.bed -b 1uM-502-X-703_005.rmbl.bed -wa -wb >common/c13_1uM.bed
bedtools intersect -a 1uM-502-X-702_005.rmbl.bed -b 1uM-502-X-703_005.rmbl.bed -wa -wb >common/c23_1uM.bed



############
#######基于Gpb_R42运行！！！！！！！！！！！
############

library(data.table)

common<-as.data.frame(fread("/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/c12_0uM.bed"))
for(i in (1:length(rownames(common)))){
    a<-common$V2[i]
    b<-common$V3[i]
    c<-common$V12[i]
    d<-common$V13[i]
    chr<-common$V1[i]
    common$cp_chr[i]<-chr
    common$cp_s[i]<-min(a,c)
    common$cp_e[i]<-max(b,d)
}

common<-common[order(as.integer(common$cp_s)),]
common<-common[order(common$cp_chr),]
common$cp_s<-as.integer(common$cp_s)
write.table(unique(common[,c("cp_chr","cp_s","cp_e")]),file='/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/common_peaks_c12_0uM.bed',quote = F,sep = "\t",row.names = F,
            col.names = F)

common<-as.data.frame(fread("/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/c13_0uM.bed"))
for(i in (1:length(rownames(common)))){
    a<-common$V2[i]
    b<-common$V3[i]
    c<-common$V12[i]
    d<-common$V13[i]
    chr<-common$V1[i]
    common$cp_chr[i]<-chr
    common$cp_s[i]<-min(a,c)
    common$cp_e[i]<-max(b,d)
}

common<-common[order(as.integer(common$cp_s)),]
common<-common[order(common$cp_chr),]
common$cp_s<-as.integer(common$cp_s)
write.table(unique(common[,c("cp_chr","cp_s","cp_e")]),file='/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/common_peaks_c13_0uM.bed',quote = F,sep = "\t",row.names = F,
            col.names = F)

common<-as.data.frame(fread("/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/c23_0uM.bed"))
for(i in (1:length(rownames(common)))){
    a<-common$V2[i]
    b<-common$V3[i]
    c<-common$V12[i]
    d<-common$V13[i]
    chr<-common$V1[i]
    common$cp_chr[i]<-chr
    common$cp_s[i]<-min(a,c)
    common$cp_e[i]<-max(b,d)
}

common<-common[order(as.integer(common$cp_s)),]
common<-common[order(common$cp_chr),]
common$cp_s<-as.integer(common$cp_s)
write.table(unique(common[,c("cp_chr","cp_s","cp_e")]),file='/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/common_peaks_c23_0uM.bed',quote = F,sep = "\t",row.names = F,
            col.names = F)

common<-as.data.frame(fread("/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/c12_1uM.bed"))
for(i in (1:length(rownames(common)))){
    a<-common$V2[i]
    b<-common$V3[i]
    c<-common$V12[i]
    d<-common$V13[i]
    chr<-common$V1[i]
    common$cp_chr[i]<-chr
    common$cp_s[i]<-min(a,c)
    common$cp_e[i]<-max(b,d)
}

common<-common[order(as.integer(common$cp_s)),]
common<-common[order(common$cp_chr),]
common$cp_s<-as.integer(common$cp_s)
write.table(unique(common[,c("cp_chr","cp_s","cp_e")]),file='/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/common_peaks_c12_1uM.bed',quote = F,sep = "\t",row.names = F,
            col.names = F)

common<-as.data.frame(fread("/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/c13_1uM.bed"))
for(i in (1:length(rownames(common)))){
    a<-common$V2[i]
    b<-common$V3[i]
    c<-common$V12[i]
    d<-common$V13[i]
    chr<-common$V1[i]
    common$cp_chr[i]<-chr
    common$cp_s[i]<-min(a,c)
    common$cp_e[i]<-max(b,d)
}

common<-common[order(as.integer(common$cp_s)),]
common<-common[order(common$cp_chr),]
common$cp_s<-as.integer(common$cp_s)
write.table(unique(common[,c("cp_chr","cp_s","cp_e")]),file='/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/common_peaks_c13_1uM.bed',quote = F,sep = "\t",row.names = F,
            col.names = F)

common<-as.data.frame(fread("/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/c23_1uM.bed"))
for(i in (1:length(rownames(common)))){
    a<-common$V2[i]
    b<-common$V3[i]
    c<-common$V12[i]
    d<-common$V13[i]
    chr<-common$V1[i]
    common$cp_chr[i]<-chr
    common$cp_s[i]<-min(a,c)
    common$cp_e[i]<-max(b,d)
}

common<-common[order(as.integer(common$cp_s)),]
common<-common[order(common$cp_chr),]
common$cp_s<-as.integer(common$cp_s)
write.table(unique(common[,c("cp_chr","cp_s","cp_e")]),file='/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/common_peaks_c23_1uM.bed',quote = F,sep = "\t",row.names = F,
            col.names = F)

##############################

####先找CUT&Tag中保守的Pekas

#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！

##运行位置：/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common
##运行环境：base
#!/bin/bash

bedtools merge -i common_peaks_c12_0uM.bed > common_peaks_c12_0uM_dp.bed
bedtools merge -i common_peaks_c13_0uM.bed > common_peaks_c13_0uM_dp.bed
bedtools merge -i common_peaks_c23_0uM.bed > common_peaks_c23_0uM_dp.bed


bedtools merge -i common_peaks_c12_1uM.bed > common_peaks_c12_1uM_dp.bed
bedtools merge -i common_peaks_c13_1uM.bed > common_peaks_c13_1uM_dp.bed
bedtools merge -i common_peaks_c23_1uM.bed > common_peaks_c23_1uM_dp.bed


cat common_peaks_c12_0uM_dp.bed common_peaks_c13_0uM_dp.bed common_peaks_c23_0uM_dp.bed > CT_0uM_fin_peaks.bed

cat common_peaks_c12_1uM_dp.bed common_peaks_c13_1uM_dp.bed common_peaks_c23_1uM_dp.bed > CT_1uM_fin_peaks.bed

sort -k1,1V -k2,2n -k3,3n ATAC_0uM_fin_peaks.bed > CT_0uM_fin_peaks_sorted.bed
bedtools merge -i ATAC_0uM_fin_peaks_sorted.bed > CT_0uM_fin_peaks_sorted_dp.bed

sort -k1,1V -k2,2n -k3,3n ATAC_1uM_fin_peaks.bed > CT_1uM_fin_peaks_sorted.bed
bedtools merge -i ATAC_1uM_fin_peaks_sorted.bed > CT_1uM_fin_peaks_sorted_dp.bed

######ATAC-seq的callpeaks文件获取代码在《4.1 上游分析过程》
bedtools intersect -a CT_0uM_fin_peaks_sorted_dp.bed -b /data/02.private/guanpb/H4K77ac/ATAC/MACS2_callpeak/rmbl/ATAC_0um_peaks.rmbl.bed -wa -wb >c0uM.bed

bedtools intersect -a CT_1uM_fin_peaks_sorted_dp.bed -b /data/02.private/guanpb/H4K77ac/ATAC/MACS2_callpeak/rmbl/ATAC_1um_peaks.rmbl.bed -wa -wb >c1uM.bed

############
#######基于Gpb_R42运行！！！！！！！！！！！
############

library(data.table)
common<-as.data.frame(fread("/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/c0uM.bed"))
for(i in (1:length(rownames(common)))){
    a<-common$V2[i]
    b<-common$V3[i]
    c<-common$V5[i]
    d<-common$V6[i]
    chr<-common$V1[i]
    common$cp_chr[i]<-chr
    common$cp_s[i]<-min(a,c)
    common$cp_e[i]<-max(b,d)
}

common<-common[order(common$cp_s),]
common<-common[order(common$cp_chr),]
write.table(common[,c("cp_chr","cp_s","cp_e")],file='/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/CT_ATAC_0uM_common_peaks.bed',quote = F,sep = "\t",row.names = F,
            col.names = F)

common<-as.data.frame(fread("/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/c1uM.bed"))
for(i in (1:length(rownames(common)))){
    a<-common$V2[i]
    b<-common$V3[i]
    c<-common$V5[i]
    d<-common$V6[i]
    chr<-common$V1[i]
    common$cp_chr[i]<-chr
    common$cp_s[i]<-min(a,c)
    common$cp_e[i]<-max(b,d)
}

common<-common[order(common$cp_s),]
common<-common[order(common$cp_chr),]
write.table(common[,c("cp_chr","cp_s","cp_e")],file='/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common/CT_ATAC_1uM_common_peaks.bed',quote = F,sep = "\t",row.names = F,
            col.names = F)

####################

#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！

##运行位置：/data/02.private/guanpb/H4K77ac/CUTTag/SEACR_005/common
##运行环境：base
bedtools merge -i CT_ATAC_0uM_common_peaks.bed > CT_ATAC_0uM_common_peaks_dp.bed
bedtools merge -i CT_ATAC_1uM_common_peaks.bed > CT_ATAC_1uM_common_peaks_dp.bed

cat CT_ATAC_0uM_common_peaks_dp.bed CT_ATAC_1uM_common_peaks_dp.bed > combine_DiffBind_Peaks_CT.bed


############
#######基于R_base 3.6.1运行！！！！！！！！！！！
############DiffBind运算

library(DiffBind)

samples <- read.csv("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/RGFP966_Aged_CT_ATAC_DiffBind_by_CT.csv")

RGFP966 <- dba(sampleSheet=samples)
RGFP966 <- dba.count(RGFP966, bUseSummarizeOverlaps=TRUE)
plot(RGFP966)
RGFP966 <- dba.contrast(RGFP966, categories = DBA_TREATMENT,minMembers = 2)
RGFP966 <- dba.analyze(RGFP966)
hmap <- colorRampPalette(c("green", "black", "red"))(n = 13)
readscores <-dba.plotHeatmap(RGFP966, contrast=1,correlations=FALSE, scale="row", colScheme = hmap)

saveRDS(RGFP966,"CT_ATAC_combind_DiffBind_ht.rds")

svg('/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/plot/CT_ATAC_combind_DiffBind_ht.svg',width=6,height=6)

hmap <- colorRampPalette(c("#8DCED6","#DAE2F0","#E44832"))(n=10)

p<-dba.plotHeatmap(RGFP966, contrast=1,
                                           correlations=FALSE, th=1,
                                           scale="row", colScheme =hmap)
dev.off()

dba.plotMA(RGFP966, contrast = 1,bUsePval=T)
dba.plotMA(RGFP966, contrast = 1,bUsePval=F,th=0.25)
dba.plotVolcano(RGFP966,,bUsePval=T)
RGFP966 <- dba.report(RGFP966,  method=DBA_DESEQ2, contrast = 1, th=1)

saveRDS(RGFP966,"CT_ATAC_combind_DiffBind.rds")

############
#######基于Gpb_R42运行！！！！！！！！！！！
############
############对DiffBind计算出的peaks进行注释

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

object<-readRDS("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/CT_ATAC_combind_DiffBind.rds")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
object_peakAnno <- annotatePeak(object, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoPie(object_peakAnno)
object_adata<-as.data.frame(object_peakAnno)
object_adata$anno<-""
for (i in 1:length(rownames(object_adata))){
    a<-object_adata$annotation[i]
    b<-unlist(strsplit(a,''))
    n<-length(which(b=="("))
    if (n==0){
        object_adata$anno[i]<-a
    }else{
        n1<-which(b=="(")
        c<-b[1:n1-1]
        object_adata$anno[i]<-c(paste(c,collapse=''))
    }
}


#object_adata_ft<-object_adata[object_adata$anno=="Promoter "&object_adata$p.value<0.1,]
object_adata_ft<-object_adata[object_adata$p.value<0.1,]

write.csv(object_adata_ft,"H4K77ac_ATAC_RGF966_vs_DMSO_DiffBind_Aged_adata.csv")

filtered_new_dba_result_combined_HSC_anno <- 
ddply(object_adata_ft, 'SYMBOL', summarise, Fold = mean(Fold))

filtered_new_dba_result_combined_HSC_anno <- 
filtered_new_dba_result_combined_HSC_anno[!is.na(filtered_new_dba_result_combined_HSC_anno$SYMBOL), ]
filtered_new_dba_result_combined_HSC_anno <- filtered_new_dba_result_combined_HSC_anno[
    order(filtered_new_dba_result_combined_HSC_anno$Fold, decreasing = T),]

filtered_new_dba_result_combined_HSC_anno

write.table(filtered_new_dba_result_combined_HSC_anno[, c("SYMBOL", 'Fold')], 
            file = '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/input/combined_DiffBind_by_CT.rnk',
            col.names = F, quote = F, sep = '\t', row.names = F )

#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！
##运行环境：base
##运行位置：/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/out
#!/bin/bash

gmtlist='GO KEGG ReactomePathways'
mkdir ATAC_CT_combined_DiffBind_By_CT_prerankedGSEA
for gmt in ${gmtlist}
do
        gsea-cli.sh GSEAPreranked -gmx /data1/02.private/guanpb/GSEA_GMT/${gmt}.gmt \
        -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk '/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/GSEA/input/combined_DiffBind_by_CT.rnk' \
        -zip_report false -out ATAC_CT_combined_DiffBind_By_CT_prerankedGSEA -rpt_label ${gmt} \
        -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
        -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 10
done
