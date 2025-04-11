##########ATAC-seq部分的处理由禹君师兄负责 (下面也整理进来了)
###参考文献链接： doi:10.1038/s41467-022-30440-2
##GEO编号：GSE162551
###ATAC-seq下载原始fastq位置：/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data
###处理代码位置：/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_seq/code.ipynb 
####data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_seq/code2.ipynb
###

#/bin/bash
#######基于Linux base 环境运行！！！！！！！！！！！！！！！！！！！！ 
##运行位置：/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data
bowtie2_dir=/data1/02.private/dengyj/software/bowtie2-2.5.1-linux-x86_64/
for name in $(ls | grep fastq)
do
    name_sub=$(echo $name| cut -d. -f1)
    ${bowtie2_dir}bowtie2 -p 10 \
    -x /data1/00.software/01.database/cut_Tag/Mm/UCSC/mm10/Sequence/Bowtie2Index/genome \
    -q ${name_sub}.fastq -S ${name_sub}.sam
    samtools view -bS ${name_sub}.sam > ${name_sub}.bam

    samtools view -hF 4 ${name_sub}.bam | grep -vF chrM | samtools view -bS > ${name_sub}_uM.bam

    samtools sort ${name_sub}_uM.bam > sorted_${name_sub}_uM.bam
    samtools index sorted_${name_sub}_uM.bam
    ratio=".6378"
    samtools view -h -b -s 660${ratio} sorted_${name_sub}_uM.bam > ${name_sub}_d${ratio}.bam
    samtools index ${name_sub}_d${ratio}.bam    
done

library("ChIPseeker")
library("org.Mm.eg.db")
library("clusterProfiler")
library("GenomicFeatures")
library("TxDb.Mmusculus.UCSC.mm10.knownGene") 
library(data.table)

library(openxlsx)

ATAC_seq_dir <- '../GSE162551_Mm_HSC_ATAC_raw_data/'
ATAC_seq_files <- paste0(ATAC_seq_dir, dir(ATAC_seq_dir))
ATAC_seq_files <- ATAC_seq_files[grepl('rmbl', ATAC_seq_files)]

K77_cuttag_dir <- 'K77/'

for(i in dir(K77_cuttag_dir)[!grepl('processed', dir(K77_cuttag_dir))]){
    df <- fread(paste0(K77_cuttag_dir, i))
    df <- as.data.frame(df)
    file_name <- paste0('processed_',i)
    write.table(df[, 1:3], sep='\t',file = paste0(K77_cuttag_dir, file_name),
                quote = F, col.names = F, row.names = F)
}

get_pairs <- function(char){
  char <- unique(char)
  result_x <- lapply(char, function(x){
    if(which(char == x) != length(char)){
      lapply(char[(which(char == x)+ 1):length(char)], function(y){
        return(c(x, y))
      })            
    }
  })
  result_x <- unlist(result_x, recursive = F)
  result_x <- result_x[!sapply(result_x, is.null)]
  return(result_x)
}

#######基本逻辑：
#先找到在各自组内至少出现2次的peaks  之后再找年龄对应的CUT&Tag和ATAC-seq的 common peaks 
##以此来获取 CUT&Tag和ATAC-seq young 和 aged 分别的共定位的peaks

############peaks occure at least in 2 samples

ATAC_seq_young_list <- get_pairs(ATAC_seq_files[grep('609|610|611|612', ATAC_seq_files)])

ATAC_seq_young_combined <- lapply(ATAC_seq_young_list, function(tmp){
    system2('/data1/00.software/00.common/bedtools/bedtools2/bin/bedtools', 
        args = paste0('intersect -a ', tmp[1], ' -b ', 
                     tmp[2], ' -wa -wb > tmp.bed'))
    df <- read.table('tmp.bed')
    row_1 <- df[1, ]
    chr_col <- grep('chr', row_1)
    chr_col1 <- chr_col[1]
    chr_start1 <- chr_col[1]+1
    chr_end1 <- chr_col[1]+2
    
    chr_col2 <- chr_col[2]
    chr_start2 <- chr_col[2]+1
    chr_end2 <- chr_col[2]+2
    
    new_df <- data.frame('seqnames' = df[, chr_col1],'start' = pmin(df[,chr_start1], df[,chr_start2]),
                         'end' = pmax(df[,chr_end1], df[,chr_end2]))
    new_df
})



ATAC_seq_young_combined_df <- do.call(rbind, ATAC_seq_young_combined)

write.table(ATAC_seq_young_combined_df, sep='\t',
            file = 'ATAC_seq_young_unsorted.bed', quote = F, row.names = F, col.names = F)



system2('sort', args =  
        '-k1,1V -k2,2n -k3,3n ATAC_seq_young_unsorted.bed > ATAC_seq_young_sorted.bed')

system2('/data1/00.software/00.common/bedtools/bedtools2/bin/bedtools', 
    args = 'merge -i ATAC_seq_young_sorted.bed > ATAC_seq_young_intersect_final.bed')

system2('/data1/00.software/00.common/bedtools/bedtools2/bin/bedtools', 
    args = 'intersect -a ATAC_seq_young_intersect_final.bed -b K77/processed_You_005.rmbl.bed -wa -wb > young_HSC_combined.bed')

young_HSC_combined <- read.table('young_HSC_combined.bed')


row_1 <- young_HSC_combined[1, ]
chr_col <- grep('chr', row_1)
chr_col1 <- chr_col[1]
chr_start1 <- chr_col[1]+1
chr_end1 <- chr_col[1]+2

chr_col2 <- chr_col[2]
chr_start2 <- chr_col[2]+1
chr_end2 <- chr_col[2]+2

young_HSC_intersect_final <- 
data.frame('seqnames' = young_HSC_combined[, chr_col1],
           'start' = pmin(young_HSC_combined[,chr_start1], young_HSC_combined[,chr_start2]),
           'end' = pmax(young_HSC_combined[,chr_end1], young_HSC_combined[,chr_end2]))

write.table(young_HSC_intersect_final, file = 'young_HSC_intersect_final.bed', sep = '\t',
            quote = F, row.names = F, col.names = F)

young_HSC_intersect_final_anno <- as.data.frame(
    annotatePeak('young_HSC_intersect_final.bed',tssRegion=c(-3000, 3000),
                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                 annoDb="org.Mm.eg.db"))

write.table(young_HSC_intersect_final_anno, file = 'young_HSC_intersect_final_anno.bed')

ATAC_seq_old_list <- get_pairs(ATAC_seq_files[grep('627|628|629|630', ATAC_seq_files)])


ATAC_seq_old_combined <- lapply(ATAC_seq_old_list, function(tmp){
    system2('/data1/00.software/00.common/bedtools/bedtools2/bin/bedtools', 
        args = paste0('intersect -a ', tmp[1], ' -b ', 
                     tmp[2], ' -wa -wb > tmp.bed'))
    df <- read.table('tmp.bed')
    row_1 <- df[1, ]
    chr_col <- grep('chr', row_1)
    chr_col1 <- chr_col[1]
    chr_start1 <- chr_col[1]+1
    chr_end1 <- chr_col[1]+2
    
    chr_col2 <- chr_col[2]
    chr_start2 <- chr_col[2]+1
    chr_end2 <- chr_col[2]+2
    
    new_df <- data.frame('seqnames' = df[, chr_col1],'start' = pmin(df[,chr_start1], df[,chr_start2]),
                         'end' = pmax(df[,chr_end1], df[,chr_end2]))
    new_df
})



ATAC_seq_old_combined_df <- do.call(rbind, ATAC_seq_old_combined)

write.table(ATAC_seq_old_combined_df, sep='\t',
            file = 'ATAC_seq_old_unsorted.bed', quote = F, row.names = F, col.names = F)



system2('sort', args =  
        '-k1,1V -k2,2n -k3,3n ATAC_seq_old_unsorted.bed > ATAC_seq_old_sorted.bed')

system2('/data1/00.software/00.common/bedtools/bedtools2/bin/bedtools', 
    args = 'merge -i ATAC_seq_old_sorted.bed > ATAC_seq_old_intersect_final.bed')

system2('/data1/00.software/00.common/bedtools/bedtools2/bin/bedtools', 
    args = 'intersect -a ATAC_seq_old_intersect_final.bed -b K77/processed_Old_005.rmbl.bed -wa -wb > old_HSC_combined.bed')

old_HSC_combined <- read.table('old_HSC_combined.bed')


row_1 <- old_HSC_combined[1, ]
chr_col <- grep('chr', row_1)
chr_col1 <- chr_col[1]
chr_start1 <- chr_col[1]+1
chr_end1 <- chr_col[1]+2

chr_col2 <- chr_col[2]
chr_start2 <- chr_col[2]+1
chr_end2 <- chr_col[2]+2

old_HSC_intersect_final <- 
data.frame('seqnames' = old_HSC_combined[, chr_col1],
           'start' = pmin(old_HSC_combined[,chr_start1], old_HSC_combined[,chr_start2]),
           'end' = pmax(old_HSC_combined[,chr_end1], old_HSC_combined[,chr_end2]))

write.table(old_HSC_intersect_final, file = 'old_HSC_intersect_final.bed', sep = '\t',
            quote = F, row.names = F, col.names = F)

old_HSC_intersect_final_anno <- as.data.frame(
    annotatePeak('old_HSC_intersect_final.bed',tssRegion=c(-3000, 3000),
                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                 annoDb="org.Mm.eg.db"))

write.table(old_HSC_intersect_final_anno, file = 'old_HSC_intersect_final_anno.bed')

####先找出年轻组最少有两个样品的peaks，再找出年老组最少有两个样品的peaks，在合并H4k77
HSC_intersect_final <- rbind(young_HSC_intersect_final, old_HSC_intersect_final)

write.table(HSC_intersect_final, sep='\t',
            file = 'HSC_intersect_final_unsorted.bed', quote = F, row.names = F, col.names = F)

system2('sort', args =  
        '-k1,1V -k2,2n -k3,3n HSC_intersect_final_unsorted.bed > HSC_intersect_final_sorted.bed')

system2('/data1/00.software/00.common/bedtools/bedtools2/bin/bedtools', 
    args = 'merge -i HSC_intersect_final_sorted.bed > combined_HSC_final.bed')


young_H4K77 <- read.table('K77//processed_You_005.rmbl.bed')

old_H4K77 <- read.table('K77//processed_Old_005.rmbl.bed')

combined_H4K77 <- rbind(young_H4K77, old_H4K77)

write.table(combined_H4K77, sep='\t',
            file = 'combined_H4K77_unsorted.bed', quote = F, row.names = F, col.names = F)

system2('sort', args =  
        '-k1,1V -k2,2n -k3,3n combined_H4K77_unsorted.bed > combined_H4K77_sorted.bed')

system2('/data1/00.software/00.common/bedtools/bedtools2/bin/bedtools', 
    args = 'merge -i combined_H4K77_sorted.bed > combined_H4K77_final.bed')


system2('/data1/00.software/00.common/bedtools/bedtools2/bin/bedtools', 
    args = 'intersect -a young_HSC_intersect_final.bed -b old_HSC_intersect_final.bed -v > only_young_final.bed')

system2('/data1/00.software/00.common/bedtools/bedtools2/bin/bedtools', 
    args = 'intersect -a old_HSC_intersect_final.bed -b young_HSC_intersect_final.bed -v > only_old_final.bed')

library(DiffBind)

library(plyr)

ATAC_seq_bam <- paste0(ATAC_seq_dir, dir(ATAC_seq_dir))
ATAC_seq_bam <- ATAC_seq_bam[grepl('6378.bam$', ATAC_seq_bam)]

data_info <- data.frame(Tissue = c(rep('Y',4), rep('O',4)), 
                        Factor = 'HSC', 
                        Condition = c(rep('Young',4), rep('Old',4)), 
                        Treatment  ='no', 
                        Replicate= c(1:4, 1:4), 
                        bamReads = ATAC_seq_bam, 
                        ControlID = '', 
                        bamControl = '', 
                        Peaks = ATAC_seq_files, 
                        PeakCaller = 'bed'
                        )


GSE162551_dba <- dba(sampleSheet=data_info) 
GSE162551_dba <- dba.count(GSE162551_dba, bUseSummarizeOverlaps=TRUE) 


GSE162551_dba <- dba.contrast(GSE162551_dba, categories = DBA_CONDITION, minMembers = 2)
GSE162551_dba <- dba.analyze(GSE162551_dba, method = DBA_DESEQ2)

GSE162551_dba_readscores <-dba.plotHeatmap(GSE162551_dba, contrast=1,
                                           correlations=FALSE, th=1,
                                           scale="row", colScheme = colorRampPalette(c("#00C6C6","#FF5353"))(n=50))

#normCounts <- dba.peakset(GSE162551_dba, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)

result <- data.frame(dba.report(GSE162551_dba, th = 1))

write.table(result, file = 'dba_result_unsorted.bed', quote=F, row.names = F, sep = '\t')

system2('sort', args =  
        '-k1,1V -k2,2n -k3,3n dba_result_unsorted.bed > dba_result_sorted.bed')

system2('/data1/00.software/00.common/bedtools/bedtools2/bin/bedtools', 
    args = 'intersect -a dba_result_sorted.bed -b combined_HSC_final.bed -wa -wb > dba_result_combined_HSC.bed')


dba_result_combined_HSC <- read.table('dba_result_combined_HSC.bed')
colnames(dba_result_combined_HSC)[1:11] <- colnames(result)
row_1 <- dba_result_combined_HSC[1, ]
chr_col <- grep('chr', row_1)
chr_col1 <- chr_col[1]
chr_start1 <- chr_col[1]+1
chr_end1 <- chr_col[1]+2

chr_col2 <- chr_col[2]
chr_start2 <- chr_col[2]+1
chr_end2 <- chr_col[2]+2

new_dba_result_combined_HSC <- 
data.frame('seqnames' = dba_result_combined_HSC[, chr_col1],
           'start' = pmin(dba_result_combined_HSC[,chr_start1], dba_result_combined_HSC[,chr_start2]),
           'end' = pmax(dba_result_combined_HSC[,chr_end1], dba_result_combined_HSC[,chr_end2]))
# new_dba_result_combined_HSC <- cbind(new_dba_result_combined_HSC, dba_result_combined_HSC[, 4:11])
# colnames(new_dba_result_combined_HSC) <- colnames(result)

new_dba_result_combined_HSC_anno <- as.data.frame(
    annotatePeak(makeGRangesFromDataFrame(new_dba_result_combined_HSC),tssRegion=c(-3000, 3000),
                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                 annoDb="org.Mm.eg.db"))

new_dba_result_combined_HSC_anno <- cbind(new_dba_result_combined_HSC_anno, dba_result_combined_HSC[, 4:11])

write.table(new_dba_result_combined_HSC_anno, file = 'new_dba_result_combined_HSC_anno.bed',
            quote=F, row.names = F, sep = '\t')

tmp_df <- new_dba_result_combined_HSC_anno[, 1:3]
tmp_df <- tmp_df[!duplicated(tmp_df), ]###有一些peaks在合并的过程中重复了
write.table(tmp_df, file = 'new_dba_result_combined_HSC.bed',
            quote=F, row.names = F, sep = '\t')

system2('awk', 
    args = "'($3-$2)>6 {print $1,$2,$3}' OFS='\t' new_dba_result_combined_HSC.bed > new_dba_result_combined_HSC_for_MEME.bed")

system2('/data1/00.software/00.common/bedtools/bedtools2/bin/bedtools', 
    args = "getfasta -fi /data1/02.private/dengyj/analysis/H4K77/support/mm10_rmbl_genome.fa -bed new_dba_result_combined_HSC_for_MEME.bed -fo file_for_MEME.fa")

filtered_new_dba_result_combined_HSC_anno <- 
ddply(new_dba_result_combined_HSC_anno, 'SYMBOL', summarise, Fold = mean(Fold))

filtered_new_dba_result_combined_HSC_anno <- 
filtered_new_dba_result_combined_HSC_anno[!is.na(filtered_new_dba_result_combined_HSC_anno$SYMBOL), ]
filtered_new_dba_result_combined_HSC_anno <- filtered_new_dba_result_combined_HSC_anno[
    order(filtered_new_dba_result_combined_HSC_anno$Fold, decreasing = T),]

write.table(filtered_new_dba_result_combined_HSC_anno[, c("SYMBOL", 'Fold')], 
            file = 'combined_HSC.rnk',
            col.names = F, quote = F, sep = '\t', row.names = F )

only_young_anno <- as.data.frame(
    annotatePeak('only_young_final.bed',tssRegion=c(-3000, 3000),
                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                 annoDb="org.Mm.eg.db"))

only_old_anno <- as.data.frame(
    annotatePeak('only_old_final.bed',tssRegion=c(-3000, 3000),
                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                 annoDb="org.Mm.eg.db"))

go <- read.gmt('/data1/02.private/dengyj/analysis/database/GSEA/c5.all.v7.0.symbols_mouse.gmt')
kegg <- read.gmt('/data1/02.private/dengyj/analysis/database/GSEA/c2.cp.kegg.v7.0.symbols_mouse.gmt')
reactome <- read.gmt('/data1/02.private/dengyj/analysis/database/GSEA/ReactomePathways_mouse.gmt')

database <- rbind(go, kegg, reactome)

only_young_enricher <- enricher(unique(only_young_anno$SYMBOL), TERM2GENE = database)

only_old_enricher <- enricher(unique(only_old_anno$SYMBOL), TERM2GENE = database)


wb <- createWorkbook()


addWorksheet(wb, 'only_young_enricher')
writeData(wb, sheet = 'only_young_enricher', rowNames = T, x = only_young_enricher)

addWorksheet(wb, 'only_old_enricher')
writeData(wb, sheet = 'only_old_enricher', rowNames = T, x = only_old_enricher)
saveWorkbook(wb, file = paste0('only_young_or_old_enricher', ".xlsx"), overwrite = TRUE)



ls('package:clusterProfiler')

K77_young <- readPeakFile("K77//You_005.rmbl.bed")

test <- readPeakFile("ATAC-seq data/GSM4954092_1_Y-HSC_1_peaks.narrowPeak")

for(i in ATAC_seq_files){
    df <- read.table(i)
    print(unique(df[, 1]))
}

class(TxDb.Mmusculus.UCSC.mm10.knownGene)

rownames(df2) <- rownames(df3) <- NULL









#########################################

##########################################
##########################################
##########################################
##########################################
######################找到H4K77ac与ATAC-seq的common peaks之后 基于ATAC进行DiffBind
###########################################
########################################
######################################

##########涉及Fig.4k
####运行位置：/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_seq

library(DiffBind)

library(plyr)

ATAC_seq_bam <- paste0(ATAC_seq_dir, dir(ATAC_seq_dir))
ATAC_seq_bam <- ATAC_seq_bam[grepl('6378.bam$', ATAC_seq_bam)]

data_info <- data.frame(Tissue = c(rep('Y',4), rep('O',4)), 
                        Factor = 'HSC', 
                        Condition = c(rep('Young',4), rep('Old',4)), 
                        Treatment  ='no', 
                        Replicate= c(1:4, 1:4), 
                        bamReads = ATAC_seq_bam, 
                        ControlID = '', 
                        bamControl = '', 
                        Peaks = ATAC_seq_files, 
                        PeakCaller = 'bed'
                        )


GSE162551_dba <- dba(sampleSheet=data_info) 
GSE162551_dba <- dba.count(GSE162551_dba, bUseSummarizeOverlaps=TRUE) 


GSE162551_dba <- dba.contrast(GSE162551_dba, categories = DBA_CONDITION, minMembers = 2)
GSE162551_dba <- dba.analyze(GSE162551_dba, method = DBA_DESEQ2)

GSE162551_dba_readscores <-dba.plotHeatmap(GSE162551_dba, contrast=1,
                                           correlations=FALSE, th=1,
                                           scale="row", colScheme = colorRampPalette(c("#00C6C6","#FF5353"))(n=50))

result <- data.frame(dba.report(GSE162551_dba, th = 1))

write.table(result, file = 'dba_result_unsorted.bed', quote=F, row.names = F, sep = '\t')

system2('sort', args =  
        '-k1,1V -k2,2n -k3,3n dba_result_unsorted.bed > dba_result_sorted.bed')

###########################################

###涉及Fig.4k
library(DiffBind)
a<-readRDS("/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_seq/GSE162551_dba.rds")
hmap <- colorRampPalette(c("#0080C0","#FF5353"))(n=10)
p<-dba.plotHeatmap(a, contrast=1,
                                           correlations=FALSE, th=1,
                                           scale="row", colScheme =hmap)

saveRDS(as.data.frame(p),"/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/H4K77ac_com_ATAC_DiffBind_data.rds")

df<-as.data.frame(dba.report(a,th=1))

saveRDS(df,"/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/DiffBind_df.rds")

#####保存对象
####saveRDS("/data1/02.private/other/GSE162551_dba.rds")

####涉及Fig.4l
#########################计算淋系和髓系发育相关基因的reads 在young和Aged的分布  
###############

############由于是基于ATAC的reads的分布计算的 因此需要生成bw文件 统计reads数

#########基于Linux base 环境运行！！！！！！！！！！！！！！！！！！！！ 

##运行环境：base
###运行位置：/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data

#!/bin/bash


SAMPLE_NAME='SRR13186609 SRR13186610 SRR13186611 SRR13186612 SRR13186627 SRR13186628 SRR13186629 SRR13186630'
for name in $SAMPLE_NAME
do
	bamCoverage -p 5 --normalizeUsing CPM -b /data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/sorted_${name}_uM.bam -o ./bwfile/sorted_${name}_uM.bw --skipNAs --effectiveGenomeSize 2652783500 
done

##########注释ATAC-seq和CUT&tag的共定位peaks的结果  并通过GO数据库 生成计算的淋系和髓系的bed文件
##运行环境 R base

common_peaks<-readRDS("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/DiffBind_df.rds")

common_peaks$start<-common_peaks$start+1

write.table(common_peaks[,c("seqnames","start","end")],"/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/H4K77ac_com_ATAC_DiffBind_rawdata.bed",row.names = F,
            col.names = F,quote = F,sep = "\t")

whole_gene_pos<-readRDS("/data1/02.private/guanpb/data/Trim68/cut_tag/whole_gene_annotation.rds")
h2m_gene<-as.data.frame(fread("/data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip"))
gmt_GO <- strsplit(readLines("/data1/02.private/guanpb/GSEA_GMT/new/c5.go.v7.5.1.symbols.gmt"), split = '\t')
geneset_names <- sapply(1:length(gmt_GO), function(x) {gmt_GO[[x]][1]})
gmt_GO <- sapply(1:length(gmt_GO), function(x) {gmt_GO[[x]][-c(1,2)]})
names(gmt_GO) <- geneset_names
for (i in (1:length(gmt_GO))){
    gmt_GO[[i]]<-h2m_gene$"Probe Set ID"[h2m_gene$"Gene Symbol"%in%gmt_GO[[i]]]
}

B_cell_genes<-gmt_GO$GOBP_REGULATION_OF_B_CELL_DIFFERENTIATION
T_cell_genes<-gmt_GO$GOBP_REGULATION_OF_T_CELL_DIFFERENTIATION
Myeloid_Cell_genes<-gmt_GO$GOBP_POSITIVE_REGULATION_OF_MYELOID_CELL_DIFFERENTIATION

all_gene_B<-whole_gene_pos[whole_gene_pos$SYMBOL%in%B_cell_genes,c("seqnames","start","end","SYMBOL")]

all_gene_T<-whole_gene_pos[whole_gene_pos$SYMBOL%in%T_cell_genes,c("seqnames","start","end","SYMBOL")]

all_gene_My<-whole_gene_pos[whole_gene_pos$SYMBOL%in%Myeloid_Cell_genes,c("seqnames","start","end","SYMBOL")]

all_gene_Ly_df<-rbind(all_gene_B,all_gene_T)

write.table(all_gene_My,"/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/Myeloid_Cell_genes.bed",row.names = F,
            col.names = F,quote = F,sep = "\t")
write.table(all_gene_T,"/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/T_Cell_genes.bed",row.names = F,
            col.names = F,quote = F,sep = "\t")
write.table(all_gene_B,"/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/B_Cell_genes.bed",row.names = F,
            col.names = F,quote = F,sep = "\t")
write.table(unique(all_gene_Ly_df),"/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/Lymphoid_Cell_genes.bed",row.names = F,
            col.names = F,quote = F,sep = "\t")

library(data.table)
myeloid_df<-as.data.frame(fread("/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/Myeloid_Cell_genes.bed"))
Lymphoid_df<-as.data.frame(fread("/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/Lymphoid_Cell_genes.bed"))
T_df<-as.data.frame(fread("/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/T_Cell_genes.bed"))
B_df<-as.data.frame(fread("/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/B_Cell_genes.bed"))

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

whole_gene_pos<-readRDS("/data1/02.private/guanpb/data/Trim68/cut_tag/whole_gene_annotation.rds")

####

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
data<-readPeakFile("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/H4K77ac_com_ATAC_DiffBind_rawdata.bed")
data <- annotatePeak(data, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
anno_data_3Kb<-as.data.frame(data)
anno_data_3Kb$anno<-""
for (i in 1:length(rownames(anno_data_3Kb))){
    a<-anno_data_3Kb$annotation[i]
    b<-unlist(strsplit(a,''))
    n<-length(which(b=="("))
    if (n==0){
        anno_data_3Kb$anno[i]<-a
    }else{
        n1<-which(b=="(")
        c<-b[1:n1-1]
        anno_data_3Kb$anno[i]<-c(paste(c,collapse=''))
    }
}
anno_data_3Kb_allgene<-anno_data_3Kb$SYMBOL
anno_data_3Kb_TSSgene<-anno_data_3Kb$SYMBOL[anno_data_3Kb$anno=="Promoter "]
###
common_peak_3Kb_TSSgene_df_Lym<-whole_gene_pos[whole_gene_pos$SYMBOL%in%intersect(anno_data_3Kb_TSSgene,Lymphoid_df$V4),c("seqnames","start","end","SYMBOL")]

common_peak_3Kb_TSSgene_df_Mye<-whole_gene_pos[whole_gene_pos$SYMBOL%in%intersect(anno_data_3Kb_TSSgene,myeloid_df$V4),c("seqnames","start","end","SYMBOL")]

write.table(common_peak_3Kb_TSSgene_df_Mye,"/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/common_peaks_3Kb_TSS_Myeloid_Cell_genes.bed",row.names = F,
            col.names = F,quote = F,sep = "\t")
write.table(common_peak_3Kb_TSSgene_df_Lym,"/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/common_peaks_3Kb_TSS_Lymphoid_Cell_genes.bed",row.names = F,
            col.names = F,quote = F,sep = "\t")
###

#####
Mye_3Kb_TSSgene_common_peaks<-anno_data_3Kb[anno_data_3Kb$anno=="Promoter "&anno_data_3Kb$SYMBOL%in%myeloid_df$V4 ,
               c("seqnames","start","end","SYMBOL")]
Lym_3Kb_TSSgene_common_peaks<-anno_data_3Kb[anno_data_3Kb$anno=="Promoter "&anno_data_3Kb$SYMBOL%in%Lymphoid_df$V4 ,
               c("seqnames","start","end","SYMBOL")]

#write.table(Mye_3Kb_TSSgene_common_peaks,"/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/Mye_3Kb_TSSgene_common_peaks.bed",row.names = F,
            col.names = F,quote = F,sep = "\t")
#write.table(Lym_3Kb_TSSgene_common_peaks,"/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/Lym_3Kb_TSSgene_common_peaks.bed",row.names = F,
            col.names = F,quote = F,sep = "\t")

##############计算reads分布

#########基于Linux base 环境运行！！！！！！！！！！！！！！！！！！！！ 

##运行环境：base

###运行位置：/data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data

#!/bin/bash

mkdir -p heatmap/mat/
mkdir -p heatmap/plot/

computeMatrix reference-point \
	--referencePoint TSS \
	-b 3000 -a 3000  -bs 100 \
	-R /data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/common_peaks_3Kb_TSS_Myeloid_Cell_genes.bed  \
	-S bwfile/sorted_SRR13186609_uM.bw bwfile/sorted_SRR13186610_uM.bw bwfile/sorted_SRR13186611_uM.bw bwfile/sorted_SRR13186612_uM.bw bwfile/sorted_SRR13186627_uM.bw bwfile/sorted_SRR13186628_uM.bw bwfile/sorted_SRR13186629_uM.bw bwfile/sorted_SRR13186630_uM.bw \
	-o heatmap/mat/Myeloid_cell_3Kb.mat.gz -p 8 --missingDataAsZero \
	--outFileSortedRegions heatmap/mat/regions_Myeloid_cell_3Kb.bed


	plotProfile -m heatmap/mat/Myeloid_cell_3Kb.mat.gz  \
	-out heatmap/plot/Myeloid_cell_3Kb.png \
	--dpi 600 --plotWidth 8 \
	--outFileNameData heatmap/mat/Myeloid_cell_3Kb.mat_x.tab


computeMatrix reference-point \
	--referencePoint TSS \
	-b 3000 -a 3000  -bs 100 \
	-R /data1/02.private/guanpb/GSE162551_Mm_HSC_ATAC_raw_data/common_peaks_3Kb_TSS_Lymphoid_Cell_genes.bed  \
	-S bwfile/sorted_SRR13186609_uM.bw bwfile/sorted_SRR13186610_uM.bw bwfile/sorted_SRR13186611_uM.bw bwfile/sorted_SRR13186612_uM.bw bwfile/sorted_SRR13186627_uM.bw bwfile/sorted_SRR13186628_uM.bw bwfile/sorted_SRR13186629_uM.bw bwfile/sorted_SRR13186630_uM.bw \
	-o heatmap/mat/Lymphoid_cell_3Kb.mat.gz -p 8 --missingDataAsZero \
	--outFileSortedRegions heatmap/mat/regions_Lymphoid_cell_3Kb.bed


	plotProfile -m heatmap/mat/Lymphoid_cell_3Kb.mat.gz  \
	-out heatmap/plot/Lymphoid_cell_3Kb.png \
	--dpi 600 --plotWidth 8 \
	--outFileNameData heatmap/mat/Lymphoid_cell_3Kb.mat_x.tab  


######
