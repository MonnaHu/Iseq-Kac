#########涉及Fig.5j
##############SEACR callpeak 之后Only_RGFP966 文件获取方法

###callpeaks为上游分析 因此位于《2.1上游处理过程》代码文件中

#######基于Linux base 环境运行！！！！！！！！！！！！！！！！！！！！ 
#运行位置：/home/hulab/RGFP966/SEACR_005
mkdir annobed
bedtools intersect -a K77-0_NO_B4_005.rmbl.bed -b K77-1_005.rmbl.bed -wa -wb >m.bed
bedtools merge -i m.bed > annobed/common_peaks_SEACR.bed
##
bedtools intersect -a K77-0_NO_B4_005.rmbl.bed -b annobed/common_peaks_005.bed -v > annobed/Only_DMSO.bed
bedtools intersect -a K77-1_005.rmbl.bed -b annobed/common_peaks_005.bed -v > annobed/Only_RGFP966.bed

library("ChIPseeker")
library("org.Mm.eg.db")
library("clusterProfiler")
library("GenomicFeatures")
library(corrplot)
library(enrichplot)
library(ReactomePA)
library(dplyr)
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(DiffBind)

RGFP966<-readPeakFile("/home/hulab/RGFP966/SEACR_005/annobed/Only_RGFP966.bed")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
DEpeakAnno <- annotatePeak(RGFP966, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoPie(DEpeakAnno)
adata<-as.data.frame(DEpeakAnno)
adata$anno<-""
for (i in 1:length(rownames(adata))){
    a<-adata$annotation[i]
    b<-unlist(strsplit(a,''))
    n<-length(which(b=="("))
    if (n==0){
        adata$anno[i]<-a
    }else{
        n1<-which(b=="(")
        c<-b[1:n1-1]
        adata$anno[i]<-c(paste(c,collapse=''))
    }
}

Pro<-read.csv("/data1/02.private/guanpb/data/Trim68/cut_tag/RGFP966/Proteins_RGPF966HSC.csv")

R_Gene<-unique(adata$SYMBOL[adata$anno=='Promoter '])
overlapgene<-intersect(Pro$Genes,R_Gene)
length(overlapgene)

df_gene<-Pro[Pro$Genes%in%overlapgene,c("Genes","median.ratio.966.NC","p.value")]
df_gene<-df_gene[order(df_gene$median.ratio.966.NC,decreasing = T),]
write.table(df_gene[,c("Genes","median.ratio.966.NC")],"/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/RGFP966/df_only_RGFP966_Proteion_GSEA.rnk",quote = F, sep ="\t", row.names = F,col.names = F)

#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！
##运行环境：base
##运行位置：/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/RGFP966/GSEA
#!/bin/bash
gmtlist='c2.cp.reactome.v7.5.1.symbols c5.go.v7.5.1.symbols Total_kegg'
mkdir Protein_only_RGFP966
for gmt in ${gmtlist}
do	
	gsea-cli.sh GSEAPreranked -gmx /data1/02.private/guanpb/GSEA_GMT/new/${gmt}.gmt \
	-collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk /data1/02.private/guanpb/data/Trim68/cut_tag/RGFP966/GSEA/df_only_RGFP966_Proteion_GSEA.rnk \
	-scoring_scheme weighted \
	-zip_report false -out Protein_only_RGFP966 -rpt_label ${gmt} \
	-chip /data1/02.private/guanpb/GSEA_GMT/new/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.5.chip \
	-create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 10 
done
