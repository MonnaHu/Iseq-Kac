library(data.table)

df_Y<-fread("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/ct_SEACR005_H4K77ac_only_Y_GSEA_matrix.txt")

df_Y

df_O<-fread("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/ct_SEACR005_H4K77ac_only_O_GSEA_matrix.txt")

df_O

library(Seurat)

library(Seurat)
HSC<-readRDS("/data1/02.private/other/11_5/object2.rds")

Idents(HSC)<-"time"

plan("multisession", workers = 16)
options(future.globals.maxSize= 4*1024^4)
markers_features_Y<-FindMarkers(HSC,`ident.1` = "2M",`ident.2` = "28M",logfc.threshold =log2(1),min.pct = 0,pseudocount.use = 0,features =df_Y$Name )

markers_features_Y

plan("multisession", workers = 16)
options(future.globals.maxSize= 4*1024^4)
markers_features_O<-FindMarkers(HSC,`ident.1` = "2M",`ident.2` = "28M",logfc.threshold =log2(1),min.pct = 0,pseudocount.use = 0,features =df_O$Name )

plan("multisession", workers = 16)
options(future.globals.maxSize= 4*1024^4)
markers_features<-FindAllMarkers(HSC,,logfc.threshold =log2(1),min.pct = 0,pseudocount.use = 0 )

markers_features

plan("sequential")

all_Y_up<-markers_features[markers_features$p_val_adj<0.05 &
                          markers_features$avg_log2FC>log2(1.5)&
                           markers_features$cluster=="2M",]

all_Y_up

all_O_Up<-markers_features[markers_features$p_val_adj<0.05 &
                          markers_features$avg_log2FC>log2(1.5)&
                           markers_features$cluster=="28M",]

all_Y_up

all_O_Up

markers_features

gene_Y<-rownames(markers_features_Y[markers_features_Y$p_val_adj<0.05 &
                          markers_features_Y$avg_log2FC>log2(1.5),])

length(gene_Y)

gene_O<-rownames(markers_features_O[markers_features_O$p_val_adj<0.05 &
                          markers_features_O$avg_log2FC<(-log2(1.5)),])

length(gene_O)

sample<-c("gene_Y","gene_O")

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
library(data.table)

library(openxlsx)
wb <- createWorkbook()
for (i in sample){
    gene_list<-get(i)
    transID = bitr(gene_list,
                    fromType="SYMBOL",
                    toType=c("ENSEMBL", "ENTREZID"),
                    OrgDb="org.Mm.eg.db")
    ego <- enrichGO(transID$ENTREZID, 
                   keyType = "ENTREZID", 
						   OrgDb = org.Mm.eg.db, 
                    	   ont = "BP", 
                    	   pAdjustMethod = "BH", 
                    	   pvalueCutoff = 0.1,
                    	   qvalueCutoff = 0.1,, 
                    	   readable = TRUE)
    name1<-paste0(i,"_GO")
    addWorksheet(wb,name1)
    writeData(wb, sheet =name1, 
				x = as.data.frame(ego))
#dotplot(ego,showCategory=50)
#write.table(as.data.frame(ego),file=paste0("/data4/cut_TAG/HIF2_Gata3/annotext/pathway319/",i,"_go.xlsx"),sep="\t",col.names=T,row.names =F,quote = F)
    kegg <- enrichKEGG(transID$ENTREZID,
								organism ="mmu",
								pvalueCutoff = 0.1,
								qvalueCutoff = 0.1,
								use_internal_data =FALSE)
    KEGG<- setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
    name2<-paste0(i,"_KEGG")
    addWorksheet(wb, name2)
    writeData(wb, sheet = name2, 
				x = as.data.frame(KEGG))
#dotplot(kegg, showCategory = 20)
#write.table(as.data.frame(KEGG),file=paste0("/data4/cut_TAG/HIF2_Gata3/annotext/pathway319/",i,"_kegg.xlsx"),sep="\t",col.names=T,row.names =F,quote = F)
    Reactome<-enrichPathway(transID$ENTREZID,
										readable=T,
										pvalueCutoff = 0.1,
										qvalueCutoff = 0.1,
										maxGSSize = 1000,
                       organism = "mouse")
    name3<-paste0(i,"_Rtm")
    addWorksheet(wb, name3)
    writeData(wb, sheet = name3, 
				x = as.data.frame(Reactome))
#dotplot(Reactome,showCategory=20)
#write.table(as.data.frame(Reactome),file=paste0("/data4/cut_TAG/HIF2_Gata3/annotext/pathway319/",i,"_Rtm.xlsx"),sep="\t",col.names=T,row.names =F,quote = F)

}

saveWorkbook(wb, file = "ORA_0310.xlsx", overwrite = TRUE)







df_pathway

df_pathway$log10_P_value<-""
for( i in (1:nrow(df_pathway))){
    if(df_pathway$class[i]=="Young"){
            df_pathway$log10_P_value<-
    }else{
        df_pathway$log10_P_value<-
    }
}

NAME<-c("lymphocyte differentiation","T cell differentiation","")

Young_name<-c('lymphocyte differentiation','T cell differentiation','regulation of T cell activation','positive regulation of T cell proliferation','B cell differentiation','B cell activation','regulation of Notch signaling pathway','Wnt signaling pathway','TOR signaling',
              'phosphatidylinositol 3-kinase signaling','phosphatidylinositol-mediated signaling','positive regulation of kinase activity','regulation of GTPase activity',
              'protein acetylation','histone modification','regulation of DNA-templated transcription elongation','transcription by RNA polymerase I','double-strand break repair','telomere organization','nucleocytoplasmic transport')

rev(df_pathway$name[df_pathway$class=="Young"])
rev(df_pathway$name[df_pathway$class=="Aged"])

Aged_name<-c('regulation of calcium ion transmembrane transport',
'regulation of metal ion transport','regulation of actin cytoskeleton organization',
'carbohydrate derivative catabolic process',
'membrane lipid metabolic process','glutathione metabolic process',
'histone H4-K20 methylation',
'regulation of protein-containing complex assembly',
'regulation of response to DNA damage stimulus',
'immunoglobulin production involved in immunoglobulin-mediated immune response')

df_pathway

df_pathway[order(df_pathway$),]

##############
library(openxlsx)
library(ggplot2)
df_pathway<-read.xlsx("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/ORA analysis_z_choose.xlsx",sheet="Sheet2")

df_pathway$GR<-""
for( i in (1:nrow(df_pathway))){
    if(df_pathway$class[i]=="Young"){
            df_pathway$GR[i]<-(eval(parse(text = df_pathway$GeneRatio[i])))
    }else{
        df_pathway$GR[i]<-(-eval(parse(text = df_pathway$GeneRatio[i])))
    }
}
df_pathway$GR<-as.numeric(df_pathway$GR)
df_pathway$class<-factor(df_pathway$class,levels=c("Young","Aged"))
options(repr.plot.width=20, repr.plot.height=12)
ggplot(df_pathway,aes(x=GR,y=name))+
    geom_bar(stat='identity',aes(fill=class))+
scale_y_discrete(limits = rev(c(Young_name,
                                Aged_name)), 
                 labels = rev(c(Young_name,
                                Aged_name)))+
#scale_colour_gradientn(colours = c("#D26B44","white","blue"),breaks=c(0,0.01,0.05))+
theme(axis.text.x = element_text(size=20,angle = 60,hjust = 1),
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
    xlab("GeneRatio")+
 scale_fill_manual(values = c("Young"="#4bc1ee","Aged"="#f29438"))+
#scale_size_continuous(range = c(9,5),guide = FALSE,limits=c(0,0.1),breaks=c(0,0.05,0.1))+
    #scale_x_continuous(limits =c(0,0.06) , breaks=c(0,75,150,225))+
NULL
ggsave("/data1/02.private/guanpb/Jupyter_R_4.2.0/old_note/old_analysis/cut_tag/H4K77_1st/new analysis/NChmB_rv/plot/DEGs_H4K77ac_ORA_barplot.svg",height =12 ,width = 22,dpi=600)


