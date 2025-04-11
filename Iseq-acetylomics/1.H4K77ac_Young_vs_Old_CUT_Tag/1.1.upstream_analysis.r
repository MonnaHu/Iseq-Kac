###涉及文章20250204 revised版本的 Fig 4a-4m
一共 3例H4K77ac-Young HSC 3例 H4K77ac-Old HSC （物种：小鼠）
包括质控 比对 去重复 整合 peak-calling 获取共定位peaks 注释 联合转录进行GSEA文件获取
#####原始fastq
数据位置：/data4/cuttag_rawdata/HSC/H4k77ac_1st_rawdata




#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！

##直接使用的公司返回来的clean的fastq进行参考基因组比对

###运行位置：/data4/cut_TAG/H4K77ac

#!/bin/bash
SIMPLE_NAME='O1 O2 O3 Y1 Y2 Y3'
mkdir sam
mkdir txt
mkdir -p bam/bwsortbam
cores=8
for name in $SIMPLE_NAME
do

        bowtie2 --very-sensitive-local --no-mixed --no-discordant -X 1000 -x /data4/cut_TAG/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome \
        -1 /data4/12.18datatocopy/data——toshare/201123_X591_0375_AHF225CCX2/002/${name}-*.R1.clean.fastq.gz \
        -2 /data4/12.18datatocopy/data——toshare/201123_X591_0375_AHF225CCX2/002/${name}-*.R2.clean.fastq.gz -p 8 -S sam/${name}.sam > txt/${name}_bowtie2.txt
        samtools view -o bam/${name}.bam sam/${name}.sam
        samtools sort -o bam/bwsortbam/${name}.sorted.bam bam/${name}.bam
done

#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！

##去重复

#workdir='//data4/cut_TAG/H4K77ac/rmdup'

#!/bin/bash
#O1 O2 O3 Y1 Y2 Y3
###alignment refgene#####
mkdir rmdup
mkdir sortbam_name
SIMPLE_NAME='O1 O2 O3 Y1 Y2 Y3'

for name in $SIMPLE_NAME
do
        samtools sort -n bam/${name}.bam -o sortbam_name/${name}.sorted.bam
        samtools fixmate -m sortbam_name/${name}.sorted.bam rmdup/${name}.fixmate.bam
        samtools sort  rmdup/${name}.fixmate.bam -o rmdup/${name}.positionsort.bam
        samtools markdup -r rmdup/${name}.positionsort.bam rmdup/${name}.sorted.markdup.bam -s -f txt/${name}.rmdup.txt
done


#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！ 

#######合并样本 3个Young 3个Old （单个样本已经去重完毕）

#运行位置：/data4/cut_TAG/H4K77ac/rmdup

###合并
#!/bin/bash 
samtools merge Old_sample.sorted.markdup.bam \
O1.sorted.markdup.bam \
O2.sorted.markdup.bam \
O3.sorted.markdup.bam


samtools merge You_sample.sorted.markdup.bam \
Y1.sorted.markdup.bam \
Y2.sorted.markdup.bam \
Y3.sorted.markdup.bam


#######生成bw格式的文件 绘制peak的丰度热图使用
mkdir bw
SAMPLE="Old_sample You_sample"

for name in $SAMPLE
do
        samtools index ${name}.sorted.markdup.bam
        bamCoverage -b ${name}.sorted.markdup.bam -o bw/${name}.sorted.markdup.bw --extendReads --normalizeUsing RPGC --skipNAs --effectiveGenomeSize 2652783500
done

SIMPLE_NAME='O1 O2 O3 Y1 Y2 Y3'

for name in $SAMPLE
do
        samtools index ${name}.sorted.markdup.bam
        bamCoverage -b ${name}.sorted.markdup.bam -o bw/${name}.sorted.markdup.bw --extendReads --normalizeUsing RPGC --skipNAs --effectiveGenomeSize 2652783500
done


#######基于Linux环境运行！！！！！！！！！！！！！！！！！！！！ 

###运行位置：/data4/cut_TAG/H4K77ac/bedgraph_callpeakSEACR/newcall

###peak calling

#!/bin/bash
SAMPLE_NAME="Old You"
mmbt='/data1/00.software/01.database/peak_blacklist/mm10-blacklist.v2.bed'
mkdir bgfile
for name in $SAMPLE_NAME
do
        bedtools genomecov -bg -ibam /data4/cut_TAG/H4K77ac/rmdup/${name}_sample.sorted.markdup.bam > bgfile/${name}.bedGraph
done
mkdir SEACR_005
SEACR="/data1/00.software/00.common/SEACR/SEACR-1.3/SEACR_1.3.sh"
for name in $SAMPLE_NAME
do
        bash $SEACR bgfile/${name}.bedGraph 0.05 non stringent SEACR_005/${name}_005
        bedtools intersect -a SEACR_005/${name}_005.stringent.bed -b ${mmbt} -v >SEACR_005/${name}_005.rmbl.bed
done

##############获取only young peaks and only old peaks的方法
#######基于Linux base 环境运行！！！！！！！！！！！！！！！！！！！！ 
##运行位置：/data4/cut_TAG/H4K77ac/use/old/bed
bedfile="/data4/cut_TAG/H4K77ac/use"
bamfile="/data4/cut_TAG/H4K77ac/rmdup"
mmbt='/data1/00.software/01.database/peak_blacklist/mm10-blacklist.v2.bed'
macs2 callpeak -t ${bamfile}/Old.sorted.markdup.bam -f BAMPE -g mm -p 1e-5 --keep-dup 1 -n ${bedfile}/old_sample
bedtools intersect -a ${bedfile}/old_sample_peaks.narrowPeak -b ${mmbt} -v >${bedfile}/old.dup1_peaks.rmbl.narrowPeak

macs2 callpeak -t ${bamfile}/You.sorted.markdup.bam -f BAMPE -g mm -p 1e-5 --keep-dup 1 -n ${bedfile}/You_sample
bedtools intersect -a ${bedfile}/You_sample_peaks.narrowPeak -b ${mmbt} -v >${bedfile}/You.dup1_peaks.rmbl.narrowPeak
##运行位置：/data4/cut_TAG/H4K77ac/use/old/bed
bedtools intersect -a ${bedfile}/You.dup1_peaks.rmbl.narrowPeak -b ${bedfile}/old.dup1_peaks.rmbl.narrowPeak-wa -wb >c1.bed
bedtools merge -i c1.bed > merge.common.bed
##
bedtools intersect -a ${bedfile}/You.dup1_peaks.rmbl.narrowPeak -b merge.common.bed -v > merge.Y.bed
bedtools intersect -a ${bedfile}/old.dup1_peaks.rmbl.narrowPeak -b merge.common.bed -v > merge.O.bed






###############
#######涉及Extended Data Fig.8 
###############

#########淋系和髓系基因集的生成

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
library(data.table)
####callpeaks  在上述代码块中提供了获取办法
all_Old<-readPeakFile("/data4/cut_TAG/H4K77ac/bedgraph_callpeakSEACR/newcall/SEACR_005/Old_005.rmbl.bed")
all_Young<-readPeakFile("/data4/cut_TAG/H4K77ac/bedgraph_callpeakSEACR/newcall/SEACR_005/You_005.rmbl.bed")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
all_OldpeakAnno <- annotatePeak(all_Old, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
all_YoungpeakAnno <- annotatePeak(all_Young, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
all_Y_adata<-as.data.frame(all_YoungpeakAnno)
all_O_adata<-as.data.frame(all_OldpeakAnno)

Young_anno<-all_Y_adata

whole_gene_pos<-readRDS("/data1/02.private/guanpb/data/Trim68/cut_tag/whole_gene_annotation.rds")

h2m_gene<-as.data.frame(fread("/data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip"))

gmt_GO <- strsplit(readLines("/data1/02.private/guanpb/GSEA_GMT/new/c5.go.v7.5.1.symbols.gmt"), split = '\t')
geneset_names <- sapply(1:length(gmt_GO), function(x) {gmt_GO[[x]][1]})
gmt_GO <- sapply(1:length(gmt_GO), function(x) {gmt_GO[[x]][-c(1,2)]})
names(gmt_GO) <- geneset_names
for (i in (1:length(gmt_GO))){
    gmt_GO[[i]]<-h2m_gene$"Probe Set ID"[h2m_gene$"Gene Symbol"%in%gmt_GO[[i]]]
}

for (i in 1:length(rownames(Young_anno))){
    a<-Young_anno$annotation[i]
    b<-unlist(strsplit(a,''))
    n<-length(which(b=="("))
    if (n==0){
        Young_anno$anno[i]<-a
    }else{
        n1<-which(b=="(")
        c<-b[1:n1-1]
        Young_anno$anno[i]<-c(paste(c,collapse=''))
    }
}

Myeloid_Cell_genes<-gmt_GO$GOBP_MYELOID_CELL_DIFFERENTIATION
Lymphoid_Cell_genes<-c(gmt_GO$GOBP_B_CELL_DIFFERENTIATION,gmt_GO$GOBP_T_CELL_DIFFERENTIATION)

Myeloid_Cell_genes_ft<-intersect(Young_anno$SYMBOL[Young_anno$anno=="Promoter "],Myeloid_Cell_genes)
Lymphoid_Cell_genes_ft<-intersect(Young_anno$SYMBOL[Young_anno$anno=="Promoter "],Lymphoid_Cell_genes)

all_gene_M<-whole_gene_pos[whole_gene_pos$SYMBOL%in%Myeloid_Cell_genes_ft,c("seqnames","start","end","SYMBOL")]

all_gene_L<-whole_gene_pos[whole_gene_pos$SYMBOL%in%Lymphoid_Cell_genes_ft,c("seqnames","start","end","SYMBOL")]

write.table(all_gene_M,"/data4/cut_TAG/H4K77ac/use/Myeloid_Cell_genes.bed",row.names = F,
            col.names = F,quote = F,sep = "\t")
write.table(all_gene_L,"/data4/cut_TAG/H4K77ac/use/Lymphoid_Cell_genes.bed",row.names = F,
            col.names = F,quote = F,sep = "\t")
