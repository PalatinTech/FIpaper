# This R script is for processing snRNA-seq data for Pl8177 50ug, 100ug, placebo, mesalazine and sham
# experiment one

# Pl8177 50ug vs Placebo (new +old+repeat)

# Pl8177 50ug vs Placebo (new+old+repeat) and Sham  (round1)

# Mesalazine vs Placebo and Sham (round 1+round2)

library(snRNAanalysis)
library(clusterProfiler)

#read new and old objects 
setwd("/mnt2/rat_colon_scdata_UC/codes/FI_paper/")

readobject<-function(input,inptype,dir,subdir)
{
  DR=NULL
  input.data=NULL
  data=NULL
  
  # Load the dataset
  path=paste(dir,subdir,"/filtered_feature_bc_matrix/",sep="")
  input.data<-Read10X(data.dir=path)
  index<-grep(rownames(input.data),pattern="^mr")
  tmp<-data.frame(name=rownames(input.data)[index])
  mtgenes<-read.table("/mnt2/rat_colon_scdata_UC/codes/Rat_Mtgenes_RGD.txt")
  colnames(mtgenes)<-c("name")
  tmp2<-data.frame(name=toupper(mtgenes$name))
  #mtgenes$name<-toupper(mtgenes$name)
  mtgenes<-do.call("rbind",list(mtgenes,tmp,tmp2))
  #mtgenes$name<-toupper(mtgenes$name)
  #mtgenes<-rbind(mtgenes,tmp)
  index<-which(rownames(input.data)%in%mtgenes$name)
  data<-input.data[-index,]
  
  DR <- CreateSeuratObject(counts = data, project = as.character(input), min.cells = 3)
  cat("NUmber of features(genes) and UMI (transcripts) per nuclei ",as.character(input),"\t: ",dim(DR)[1],"\t",dim(DR)[2],"\n:")
  DR<-AddMetaData(DR,metadata=as.character(input),col.name='sample' )
  DR<-AddMetaData(DR,metadata=as.character(inptype),col.name='type' )
  
  #plots
  pdf(paste("plots/",as.character(input),"_VlnPlot.pdf",sep=""),width = 12,height=8)
  plot1<-VlnPlot(DR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(plot1)
    dev.off()
  
  pdf(paste("plots/",as.character(input),"_Scatter.pdf",sep=""),width=12,height=8)
  plot1 <- FeatureScatter(DR, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1)
  dev.off()
  
  DR1 <- subset(DR, subset = nFeature_RNA > 200 )
  
  cat("Number of gene and cells before filtering\t",as.character(input),"\t: ",dim(DR)[1],"\t",dim(DR)[2],
      "\tNumber of genes and cells with >200 genes\t",dim(DR1)[1],"\t",dim(DR1)[2],"\n")
  return(DR1)
}

#sink("analysis_log.txt")
sampledesc<-read.table("/mnt2/rat_colon_scdata_UC/round2/codes/oldnew_sample_decsription.txt",header=T,sep="\t")
listall=sampledesc$Sample.ID[1]
assign(sampledesc$Sample[1],readobject(sampledesc$ID[1],sampledesc$Treatment[1],sampledesc$dir[1],sampledesc$subdir[1]))
for(i in 2:nrow(sampledesc))
{
  assign(sampledesc$Sample[i],readobject(sampledesc$ID[i],sampledesc$Treatment[i],sampledesc$dir[i],sampledesc$subdir[i]))
  listall=append(listall,sampledesc$Sample[i])
}
#sink()
#closeAllConnections()


# Placebo
placebo<-merge(GR2_2,y=c(GR2_3,GR2_4),add.cell.ids=c("GR22","GR23","GR24"))
saveRDS(placebo,"Rdata/placebo.RData")

# PL8177 100ug
PL8177_100<-merge(GR5_1,y=c(GR5_3),add.cell.ids=c("GR51","GR53"))

# PL8177 50ug
PL8177_50<-merge(GR4_1,y=c(GR4_3,GR4_4),add.cell.ids=c("GR41","GR43","GR44"))
saveRDS(PL8177_50,"Rdata/PL8177_50.RData")

# Mesalazine
mesal<-merge(GR7_7,y=c(GRF7_5,GRF7_6),add.cell.ids=c("GR73","GR75","GR76"))

#sham
sham<-merge(GR1_2,y=c(GR1_3,GR1_5),add.cell.ids=c("GR12","GR13","GR15"))
saveRDS(sham,"Rdata/sham.RData")

#normalization/integration and visualization
#placebo_norm<-normobj(placebo)
placebo_int<-intanalysis(placebo,"placebo_int")
DimPlot(placebo_int,split.by = "sample",reduction = "umap",label=TRUE)


#PL8177_100_norm<-normobj(PL8177_100)
PL8177_100_int<-intanalysis(PL8177_100,"PL8177_100_int")
DimPlot(PL8177_100_int,split.by = "sample",reduction = "umap",label=TRUE)

#PL8177_50_norm<-normobj(PL8177_50)
PL8177_50_int<-intanalysis(PL8177_50,"PL8177_50_int")
DimPlot(PL8177_50_int,split.by = "sample",reduction = "umap",label=TRUE)

#mesal_norm<-normobj(mesal)
mesal_int<-intanalysis(mesal,"mesal_int")
DimPlot(mesal_int,split.by = "sample",reduction = "umap",label=TRUE)

#sham
sham_norm<-normobj(sham)
sham_int<-intanalysis(sham,"sham_int")
DimPlot(sham_int,split.by = "sample",reduction = "umap",label=TRUE)

rm(list=ls(pattern="^GR"))


#merge treatment groups
#PL8177 50ug, placebo and sham

shamplacebo8177_50<-merge(sham_int, y=c(placebo_int,PL8177_50_int),add.cell.ids=c("Sham","Placebo","PL8177_50"))
shamplacebo8177_50<-normobj(shamplacebo8177_50)
DimPlot(shamplacebo8177_50,group.by = "type",label = TRUE)
DimPlot(shamplacebo8177_50,split.by = "sample",label=TRUE)

DefaultAssay(shamplacebo8177_50)<-"RNA"
# find markers in the seurat object
#markers for enterocyte (https://www.sciencedirect.com/science/article/pii/S2352345X21000266)
VlnPlot(shamplacebo8177_50,c("Epcam","Fabp1","Cldn3","Krt8","Ifi27","Phgr1","Cldn7","Krt19"))

#markers for Tuft cells
VlnPlot(shamplacebo8177_50,c("Lrmp","Dclk1","Cd24a","Trpm5"))

#markers for Goblet: Elf3, Cldn4, Itln1, Spink1, Tff3, Manf, Ccl9, Clca3, Fcgbp
VlnPlot(shamplacebo8177_50,c("Cldn4","Fcgbp","Slc12a8","Spdef"))

#markers for enterocyte: Mep1a, Fgf15, Cbr1, Ephx2, Clec2h

#markers for fibroblast:
VlnPlot(shamplacebo8177_50,c("Col1a1","Col1a2","Col6a1","Col6a2"))

VlnPlot(shamplacebo8177_50,c("Ccl5","Gzmb","Klrd1","Gzma","Cd8b","Cd8a"))

#markers for Cd8+  T cells:
VlnPlot(shamplacebo8177_50,c("Nkg7","Ccl5","Skap1"))

#markers for naive B
VlnPlot(shamplacebo8177_50,c("Cd79a","Cd19","Cd38","Cxcr4"))

#markers for Plasma
VlnPlot(shamplacebo8177_50,c("Cd79a","Cd27","Spdc1","Igj"))

#markers for mast cells
VlnPlot(shamplacebo8177_50,c("Cd69","Cd44"))

#markers for monocytes
VlnPlot(shamplacebo8177_50,c("Lyz","Cd14","Hla-dra","Tpsab1"))

#markers for CD4+ memory cells
VlnPlot(shamplacebo8177_50,c("Il7r"))
#markers for Glial cells
VlnPlot(shamplacebo8177_50,c("S100b"))

#Cenpf: t cell
#Ptprc: innate lymphoid cells
#Arf4, Tpd52: Plasma
#Ebf1 - b cells
#Lrmp- Germinal center B cells
#Cadps- mast cells not sure
#Syt1- follicular

markerfunc<-function(input,type)
{
  #find marker genes and annotate clusters
  DefaultAssay(input) <- "RNA"
  input.markers <- FindAllMarkers(input,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  return(input.markers)
}
shamplacebo8177_50.marker<-markerfunc(shamplacebo8177_50,"shamplacebo8177_50")

#cluster 13: Myh11, Actg2 Smooth muscle cells
#cluster 18: Igfbp6, Ccdc80: fibriblast
#cluster17:Chga Enteroendocrine
#cluster 7: Samd9, Herc6, Stat1 : used protein Atlas for definign cell types as T cells
#cluster12: Prox1: Enteroendocrine progenitors
#cluster 4: Krt20: enterocytes

#immune markers
immune.markers<-read.delim("/mnt2/rat_colon_scdata_UC/round2/codes/marker/immune_marker_PMID31348891.txt",sep="\t")
shamplacebo8177_50.marker$genename<-toupper(shamplacebo8177_50.marker$gene)
shamplacebo8177_50.immune.markers<-merge(shamplacebo8177_50.marker,immune.markers,by.x=c("genename"),by.y=c("gene"))



#define cell types
shamplacebo8177_50.label<-RenameIdents(shamplacebo8177_50,
                                       "0"="Unknown",
                                       "1"="Enterocyte",
                                       "2"="Enterocyte",
                                       "3"="Enterocyte",
                                       "4"="Goblet",
                                       "5"="Enterocyte",
                                       "6"="Unknown",
                                       "7"="Enterocyte",
                                       "8"="T cells",
                                       "9"="Goblet",
                                       "10"="T cells",
                                       "11"="EEP",
                                       "12"="SM",
                                       "13"="Fibroblast",
                                       "14"="Tuft",
                                       "15"="B cells",
                                       "16"="Unknown",
                                       "17"="EE",
                                       "18"="Fibroblast")
shamplacebo8177_50.label$newcluster="NA"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="0",]$newcluster="Unknown"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="1",]$newcluster="Enterocyte"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="2",]$newcluster="Enterocyte"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="3",]$newcluster="Enterocyte"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="4",]$newcluster="Goblet"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="5",]$newcluster="Enterocyte"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="6",]$newcluster="Unknown"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="7",]$newcluster="Enterocyte"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="8",]$newcluster="T cells"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="9",]$newcluster="Goblet"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="10",]$newcluster="T cells"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="11",]$newcluster="EEP"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="12",]$newcluster="SM"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="13",]$newcluster="Fibroblast"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="14",]$newcluster="Tuft"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="15",]$newcluster="B cells"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="16",]$newcluster="Unknown"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="17",]$newcluster="EE"
shamplacebo8177_50.label@meta.data[shamplacebo8177_50.label$seurat_clusters=="18",]$newcluster="Fibroblast"





##################    FIGURE 1
pdf("plots/Fig1_Umap.pdf")
DimPlot(shamplacebo8177_50.label,
        group.by = "type",
        reduction = "umap",
        cols = c("#58409A","red","#339966"))+
  guides(color = guide_legend(override.aes = list(size=8), ncol=1,label.theme = element_text(size=16)))
dev.off()

pdf("plots/Fig1_UMAP_labelled.pdf")
DimPlot(shamplacebo8177_50.label,label.size = 8,pt.size = 1,label=TRUE,shuffle = TRUE)+NoLegend()
dev.off()

################   FIGURE 2
genes<-c("Epcam","Krt8","Cldn7","Aqp8",
         "Slc12a8","Spdef",
         "Gbp2",
         "Skap1","Ptprc","Cd247",
         "Prox1",
         "Acta2","Myh11",
         "Col6a2",
         "Lrmp","Dclk1","Trpm5",
         "Cd19","Cd74","Ms4a1","Cd22",
         "Chga","Chgb")
pdf("plots/Fig2_DotPlot.pdf")
DotPlot(shamplacebo8177_50.label,
        features = genes,
        scale.by = "size",
        col.min = 0.1,
        dot.min = 0.1)+
  theme(axis.text.x = element_text(angle = 90,size=16),
        axis.text.y = element_text(size=16))+
  xlab("")+ylab("")+scale_size(range = c(1, 8))   
dev.off()

################  FIGURE 3
toplot<-prop.table(table(Idents(shamplacebo8177_50.label),shamplacebo8177_50.label$type),2)*100
tmp<-reshape2::melt(toplot)
pdf("plots/Fig3_barplot.pdf")
ggplot(data=tmp,aes(x=Var1,y=value,fill=Var2))+
  geom_bar(stat = "identity",position = position_dodge())+scale_fill_manual(values=c("#58409A","red","#339966"))+
  theme(axis.text=element_text(size=14,color = "black",angle=90),
        panel.background = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        axis.line=element_line(color = "black"))+
  xlab("")+
  ylim(c(0,100))+
  ylab("Relative percentage of cells")
dev.off()

################ FIGURE 4




saveRDS(shamplacebo8177_50.label,"Rdata/shamplacebo8177_50.label.RData")


#find cluster markers between Sham and Placebo
#HIgh or low in sham as compared to placebo
for( num in seq(1,10) )
{
  cat(num)
  c<-data.frame(cluster=unique(Idents(shamplacebo8177_50.label)))
  t=as.character(c$cluster[num])
  name=unlist(strsplit(as.character(t)," "))[[1]]
  DefaultAssay(shamplacebo8177_50.label)<-"RNA"
  assign(paste("clusterSham",name,sep=""),FindMarkers(shamplacebo8177_50.label,
                                                      ident.1 = 'Sham',
                                                      ident.2='Placebo',
                                                      min.pct = 0.2,
                                                      logfc.threshold = 0,
                                                      #min.pct=0.25,
                                                      group.by = "type",
                                                      subset.ident = t,
                                                      test.use = "MAST"))
  fvar<-get(paste("clusterSham",name,sep=""))
  write.table(fvar,file=paste("DEG/clusterSham",name,".txt",sep=""),
              quote=F, sep="\t")
  rm(fvar)
  
}

for( num in seq(1,10) )
{
  cat(num)
  c<-data.frame(cluster=unique(Idents(shamplacebo8177_50.label)))
  t=as.character(c$cluster[num])
  name=unlist(strsplit(as.character(t)," "))[[1]]
  DefaultAssay(shamplacebo8177_50.label)<-"RNA"
  assign(paste("clusterPL8177",name,sep=""),FindMarkers(shamplacebo8177_50.label,
                                                      ident.1 = 'PL8177-50ug',
                                                      ident.2='Placebo',
                                                      min.pct = 0.2,
                                                      logfc.threshold = 0,
                                                      #min.pct=0.25,
                                                      group.by = "type",
                                                      subset.ident = t,
                                                      test.use = "MAST"))
  fvar<-get(paste("clusterPL8177",name,sep=""))
  write.table(fvar,file=paste("DEG/clusterPL8177",name,".txt",sep=""),
              quote=F, sep="\t")
  rm(fvar)
  
}
#FIND COMMON DEG IN PL8177 AND SHAM VS PLACEBO t CELLS
clusterShamT$gene<-rownames(clusterShamT)
clusterPL8177T$gene<-rownames(clusterPL8177T)

nrow(clusterShamT[abs(clusterShamT$avg_log2FC)> 1.5 & clusterShamT$p_val_adj<0.05,])
nrow(clusterPL8177T[abs(clusterPL8177T$avg_log2FC)> 1.5 & clusterPL8177T$p_val_adj<0.05,])
a<-merge(clustershamT,clusterPL8177T,by=c("gene"))%>%filter(sign(avg_log2FC.x)==sign(avg_log2FC.y))
nrow(a[a$avg_log2FC.x< (-1.5) & a$avg_log2FC.y< (-1.5) & a$p_val_adj.x<0.05 & a$p_val_adj.y<0.05,])
COMMON<-a[a$avg_log2FC.x< (-1.5) & a$avg_log2FC.y< (-1.5) & a$p_val_adj.x<0.05 & a$p_val_adj.y<0.05,]

#FIND COMMON DEG IN PL8177 AND SHAM VS PLACEBO t CELLS
clusterShamEnterocyte$gene<-rownames(clusterShamEnterocyte)
clusterPL8177Enterocyte$gene<-rownames(clusterPL8177Enterocyte)

nrow(clusterShamEnterocyte[abs(clusterShamEnterocyte$avg_log2FC)> 0.5 & clusterShamEnterocyte$p_val_adj<0.05,])
nrow(clusterPL8177Enterocyte[abs(clusterPL8177Enterocyte$avg_log2FC)> 0.5 & clusterPL8177Enterocyte$p_val_adj<0.05,])
a<-merge(clusterShamEnterocyte,clusterPL8177Enterocyte,by=c("gene"))%>%filter(sign(avg_log2FC.x)==sign(avg_log2FC.y))
nrow(a[a$avg_log2FC.x>0.5 & a$avg_log2FC.y>0.5 & a$p_val_adj.x<0.05 & a$p_val_adj.y<0.05,])
COMMON<-a[a$avg_log2FC.x>0.5 & a$avg_log2FC.y>0.5 & a$p_val_adj.x<0.05 & a$p_val_adj.y<0.05,]


plotvolcano<-function(input,title)
{
  p1<-EnhancedVolcano(input,x="avg_log2FC",y="p_val",
                      lab=rownames(input),
                      FCcutoff = 0.5,
                      pCutoff = 10e-3,  
                      xlim = c(-4,4),
                      title=title,
                      legendPosition = "",
                      subtitle = "",
                      cutoffLineType = "blank",
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,labSize = 5,
                      drawConnectors = TRUE,max.overlaps = 10,
                      col=c('black', 'black', 'black', 'red3') )
  return(p1)
}
#remove row with gene names starting with AABR*
#Volcano plot for T cells
tmp<-grep(rownames(clusterShamT),pattern="AABR")
clusterShamT<-clusterShamT[-tmp,]

tmp<-NULL
tmp<-grep(rownames(clusterPL8177T),pattern="AABR")
clusterPL8177T<-clusterPL8177T[-tmp,]

pdf("plots/DE_clusterShamTcell.pdf")
plotvolcano(clusterShamT,"T cells: Sham vs placebo")
dev.off()

pdf("plots/DE_clusterPL8177Tcell.pdf")
plotvolcano(clusterPL8177T,"T cells: PL8177-50ug vs placebo")
dev.off()

#FCutoff 1 and pvalue 10-e3
pdf("plots/DE_clusterShamEnterocyte.pdf")
plotvolcano(clusterShamEnterocyte,"Enterocyte cells: Sham vs placebo")
dev.off()

pdf("plots/DE_clusterPL8177Enterocyte.pdf")
plotvolcano(clusterPL8177Enterocyte,"Enterocyte cells: PL8177-50ug vs placebo")
dev.off()


######################################################
#######  FIGURE 5
geneboxplot<-function(object,gene,cluster,title,cols,type=NULL,dir)
{
  DefaultAssay(object)<-"RNA"
  a<-FetchData(object,vars=c(gene,cols))
  b<-a[a$newcluster==as.character(cluster),]
  tmp<-b[b$Samd9>0,]
  p<-ggplot(data=tmp,aes(x=type,y=tmp[[gene]],fill=type))+
    geom_boxplot()+
   # geom_point(aes(fill = type), size = 3, shape = 21 )+  #position = position_jitterdodge()
    scale_fill_manual(values = c("#58409A","red","#339966"))+
    theme(panel.background = element_blank(),
          axis.title = element_text(size=14),
          axis.line = element_line(color="black"),
          axis.text = element_text(color="black",size=14),
          plot.title = element_text(hjust = 0.5))+
    ylab("Expression level")+
    ggtitle(title)
  ggsave(filename=paste(as.character(dir),"/",gene,"_",cluster,"boxplot.pdf",sep=""),plot=p)
  
}

DefaultAssay(shamplacebo8177_50.label)<-"RNA"
a<-FetchData(shamplacebo8177_50.label,vars=c("Samd9","sample","type","newcluster"))
b<-FetchData(shamplacebo8177_50.label,vars=c("Herc6","sample","type","newcluster"))
c<-FetchData(shamplacebo8177_50.label,vars=c("Muc13","sample","type","newcluster"))             
d<-FetchData(shamplacebo8177_50.label,vars=c("Ndrg1","sample","type","newcluster")) 
a$gene<-"Samd9 (T cells)"
a<-a[a$newcluster=="T cells",]
b$gene<-"Herc6 (T cells)"
b<-b[b$newcluster=="T cells",]
c$gene<-"Muc13 (Enterocytes)"
c<-c[c$newcluster=="Enterocyte",]
d$gene<-"Ndrg1 (Enterocyte)"
d<-d[d$newcluster=="Enterocyte",]
colnames(a)<-c("value","sample","type","newcluster","gene")
colnames(b)<-c("value","sample","type","newcluster","gene")
colnames(c)<-c("value","sample","type","newcluster","gene")
colnames(d)<-c("value","sample","type","newcluster","gene")

final<-rbind(a,b,c,d)
final<-rbind(a,b,c)
tmp<-final[final$value>0,]
p<-ggplot(data=tmp,aes(x=type,y=value,fill=type))+
  geom_boxplot()+facet_wrap(~gene,ncol=1)+
  scale_fill_manual(values = c("#58409A","red","#339966"))+
  theme(panel.background = element_blank(),
        axis.title = element_text(size=14),
        axis.line = element_line(color="black"),
        axis.text = element_text(color="black",size=12),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 14),legend.position = "none")+
  ylab("Expression level")+xlab("")

############# FIGURE 6 GSEA
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)


###################################3
#pathway analysis using fgsea
#
##################################
library(fgsea)
library(tibble)

pathanalysis<-function(cluster,name)
{
  fdata<-data.frame(name=toupper(rownames(cluster)),FC=cluster$avg_log2FC)
  ranks<-deframe(fdata)
  write.table(ranks,file=paste("pathway_analysis/data/",name,"_rank.rnk",sep=""),
              quote = F,sep="\t",
              col.names = F)
  
  
  # Load the pathways into a named list
  pathways.hallmark <- gmtPathways("/mnt2/rat_retina_scdata_DR/codes/PL9654/PL9654_w_mito/h.all.v7.4.symbols.gmt")
  #pathways.hallmark <- gmtPathways("c2.all.v7.4.symbols.gmt")
  #pathways.hallmark<-gmtPathways("c6.all.v7.4.symbols.gmt")
  #pathways.hallmark
  
  #run fgsea
  fgseadata<-fgsea(pathways=pathways.hallmark,stats=ranks,nperm=1000)
  
  fgseadata <- fgseadata %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  # Show in a nice table:
  a<-fgseadata %>%
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
    arrange(padj) %>%
    DT::datatable()
  n=paste("pathway_analysis/",as.character(name),"_hallmark.txt",sep="")
  fwrite(fgseadata,file=n)
  
  
}

pathanalysis(clusterShamT,"clusterShamT")

#glist<-clusterShamT[clusterShamT$p_val_adj<0.05,]$avg_log2FC
glist<-clusterShamT$avg_log2FC

#names(glist)<-toupper(rownames(clusterShamT[clusterShamT$p_val_adj<0.05,]))
names(glist)<-toupper(rownames(clusterShamT))
glist<-na.omit(glist)
glist = sort(glist, decreasing = TRUE)
gse<-gseGO(geneList = glist,
           ont="BP",
           keyType = "SYMBOL",
           verbose = TRUE,
           OrgDb = org.Hs.eg.db,
           by="fgsea",
           pAdjustMethod = "BH")
clusterProfiler::dotplot(gse)
cnetplot(gse,foldChange = glist,colorEdge = T,cex_gene=0.5,cex_label_gene=0.5)

#running MsigDB
#The msigdbr R package provides Molecular Signatures Database (MSigDB) gene sets typically used with the Gene Set Enrichment Analysis (GSEA) software:
install.packages("msigdbr")
library(msigdbr)
#all_gene_sets = msigdbr(species = "Rattus norvegicus")
#getting all hallmark gene sets
msigdbr_t2g = msigdbr(species = "Homo sapiens",category = "H")%>%
  dplyr::select(gs_name,gene_symbol)
head(msigdbr_t2g)

#getting all hallmark gene sets 
#gsets<-all_gene_sets %>%
#  dplyr::filter(gs_cat == "H")
#msigdbr_t2g<-gsets%>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

#Using Enrichr for gene set enrichment analysis
em<-enricher(names(glist),TERM2GENE = msigdbr_t2g)
clusterProfiler::dotplot(em)
cnetplot(em,foldChange = glist,colorEdge = T,cex_gene=0.5,cex_label_gene=0.5)

#Using GSEA for gene set analysis
#Sham vs Placebo cluster T cells
em2<-GSEA(glist,TERM2GENE = msigdbr_t2g,pAdjustMethod = "none")
clusterProfiler::dotplot(em2)
clusterProfiler::dotplot(em2,font.size=11,title="Enriched Pathways in Tcell (Sham vs Placebo)")
a<-cnetplot(em2,foldChange = glist,colorEdge = T,cex_gene=1,cex_label_gene=0.6,cex_label_category=0.8)
em<-em2[em2$Description=="HALLMARK_INTERFERON_ALPHA_RESPONSE" | 
          em2$Description=="HALLMARK_TNFA_SIGNALING_VIA_NFKB" |
          em2$Description=="HALLMARK_INTERFERON_GAMMA_RESPONSE",asis=TRUE]
cnetplot(em,foldChange = glist,colorEdge = T,cex_gene=1,cex_label_gene=0.6,cex_label_category=0.8)
png("clusterShamT_cnet.png")
print(a)
dev.off()

#Pl8177 50ug vs Placebo T cells
glist<-clusterPL8177T$avg_log2FC
names(glist)<-toupper(rownames(clusterPL8177T))
glist<-na.omit(glist)
glist = sort(glist, decreasing = TRUE)
em2<-GSEA(glist,TERM2GENE = msigdbr_t2g,pAdjustMethod = "none")

clusterProfiler::dotplot(em2,font.size=11,title="Enriched Pathways in Tcell (PL8177-50ug vs Placebo)")
b<-cnetplot(em2,foldChange = glist,colorEdge = T,cex_gene=1,cex_label_gene=0.6,cex_label_category=0.8,shadowtext="gene")
em<-em2[em2$Description=="HALLMARK_INTERFERON_ALPHA_RESPONSE" | 
          em2$Description=="HALLMARK_TNFA_SIGNALING_VIA_NFKB" |
          em2$Description=="HALLMARK_INTERFERON_GAMMA_RESPONSE",asis=TRUE]
cnetplot(em,foldChange = glist,colorEdge = T,cex_gene=1,cex_label_gene=0.6,cex_label_category=0.8)+guides(edge_color = "none")


#plot disease vs healthy and 8177 vs disease
clusterPlaceboT<-FindMarkers(shamplacebo8177_50.label,
            ident.1 = 'Placebo',
            ident.2='Sham',
            min.pct = 0.2,
            logfc.threshold = 0,
            #min.pct=0.25,
            group.by = "type",
            subset.ident = "T cells",
            test.use = "MAST")

glist2<-clusterPlaceboT$avg_log2FC
names(glist2)<-toupper(rownames(clusterPlaceboT))
glist2<-na.omit(glist2)
glist2 = sort(glist2, decreasing = TRUE)
em3<-GSEA(glist2,TERM2GENE = msigdbr_t2g,pAdjustMethod = "none")

clusterProfiler::dotplot(em3,font.size=11,title="Enriched Pathways in Tcell (PL8177-50ug vs Placebo)")
c<-cnetplot(em3,foldChange = glist2,colorEdge = T,cex_gene=1,cex_label_gene=0.6,cex_label_category=0.8,shadowtext="gene")
em<-em3[em3$Description=="HALLMARK_INTERFERON_ALPHA_RESPONSE" | 
          em3$Description=="HALLMARK_TNFA_SIGNALING_VIA_NFKB" |
          em3$Description=="HALLMARK_INTERFERON_GAMMA_RESPONSE",asis=TRUE]
cnetplot(em,foldChange = glist2,colorEdge = T,cex_gene=1,cex_label_gene=0.6,cex_label_category=0.8)+guides(edge_color = "none")

pathanalysis(clusterLeukocytes,"clusterLeukocytes")
pathanalysis(clusterPT,"clusterPT")
pathanalysis(clusterMesangial,"clusterMesangial")
pathanalysis(clusterEndo,"clusterEndo")

#DOTPLOT FOR ENTEROCYTES SHAM vs PLACEBO and PL8177 vs PLACEBO
glist2<-clusterShamEnterocyte$avg_log2FC
names(glist2)<-toupper(rownames(clusterShamEnterocyte))
glist2<-na.omit(glist2)
glist2 = sort(glist2, decreasing = TRUE)
em3<-GSEA(glist2,TERM2GENE = msigdbr_t2g,pAdjustMethod = "none")

clusterProfiler::dotplot(em3,font.size=11,title="Enriched Pathways in Enterocytes (PL8177-50ug vs Placebo)")

#STZ vs control
pathanalysis(clusterSTZctrlPodocytes,"clusterSTZcntPodocytes")
pathanalysis(clusterSTZctrlLeukocytes,"clusterSTZcntLeukocytes")
pathanalysis(clusterSTZctrlPT,"clusterSTZcntPT")
pathanalysis(clusterSTZctrlMesangial,"clusterSTZcntMesangial")
pathanalysis(clusterSTZctrlEndo,"clusterSTZcntEndo")

#pathway plots
pathplot<-function(inpfile,outfile,cell)
{
  plotclu<-read.delim(inpfile,sep=",",header=T,check.names = F,row.names = NULL)
  pdf(paste("pathway_analysis/",outfile,".pdf",sep=""))
  p<-ggplot(plotclu[plotclu$pval<0.05,], aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste("Pathways enriched in ",cell,sep="")) +
    theme_minimal()
  print(p)
  dev.off()
}

         









    
VlnPlot(shamplacebo8177_50.label,c("Slfn4"),split.by = "type",cols=c("#7030A0","red","blue"))
geneboxplot(shamplacebo8177_50.label,"Slfn4","Tcells","Slfn4 expression in T cells",cols = c("sample","type","newcluster"),type="type",dir="plots")
geneboxplot(shamplacebo8177_50.label,"Samd9","T cells","Samd9 expression in T cells",cols = c("sample","type","newcluster"),type="type",dir="plots")
geneboxplot(shamplacebo8177_50.label,"Stat1","T cells","Stat1 expression in T cells",cols = c("sample","type","newcluster"),type="type",dir="plots")
geneboxplot(shamplacebo8177_50.label,"Herc6","Tcells","Herc6 expression in T cells",cols = c("sample","type","newcluster"),type="type",dir="plots")
geneboxplot(shamplacebo8177_50.label,"Muc13","Enterocyte","Muc13 expression in Enterocyte cells",cols = c("sample","type","newcluster"),type="type",dir="plots")
geneboxplot(shamplacebo8177_50.label,"Cd74","B cells","Cd74 expression in B cells",cols = c("sample","type","newcluster"),type="type",dir="plots")

#check overall CD24 expression
DefaultAssay(shamplacebo8177_50.label)<-"RNA"
a<-FetchData(shamplacebo8177_50.label,vars=c("Stat1","type","newcluster"))
p<-ggplot(data=a,aes(x=type,y=a$Cd24,fill=type))+
  geom_boxplot(alpha=0.08)+
  geom_point(aes(fill = type), size = 3, shape = 21, position = position_jitterdodge())+
  scale_fill_manual(values = c("#423074", "#1DA5CB","red"))+
  theme(panel.background = element_blank(),
        axis.title = element_text(size=14),
        axis.line = element_line(color="black"),
        axis.text = element_text(color="black",size=14),
        plot.title = element_text(hjust = 0.5))+
  ylab("Normalized expression ")+
  ggtitle("Overall Cd24 expression")

library(dittoSeq)
dittoSeq::dittoDotPlot(shamplacebo8177_50.label,vars=c("Stat1","Samd9","Slfn4","Herc6","Gbp3","Hk2","Slc9a3"),split.by = "type",group.by = "sample")
dittoSeq::dittoDotPlot(shamplacebo8177_50.label,vars=c("Stat1","Samd9","Slfn4","Herc6","Gbp3","Hk2","Slc9a3"),split.by = "newcluster",group.by = "type",min.color = "grey")

dittoSeq::dittoDotPlot(shamplacebo8177_50.label,vars=c("Gpx2","Muc13","Krt20","Gbp3","Stat1","Samd9","Slfn4","Herc6","Hk2","Slc9a3"),
                       split.by = "newcluster",group.by = "type",min.color = "grey")

#Pathway analysis Tcells
snRNAanalysis::pathanalysis(clusterPL8177Tcells,"Tcells8177")
snRNAanalysis::pathwayplot("pathway_analysis/Tcells8177_hallmark.txt","Tcells_8177","Tcells")

snRNAanalysis::pathanalysis(clusterShamTcells,"Tcellssham")
snRNAanalysis::pathwayplot("pathway_analysis/Tcellssham_hallmark.txt","Tcells_Sham","Tcells")


#Pathway analysis Enetrocyte
snRNAanalysis::pathanalysis(clusterPL8177Enterocyte,"Enterocyte8177")
snRNAanalysis::pathwayplot("pathway_analysis/Enterocyte8177_hallmark.txt","Enterocyte_8177","Enterocyte")

snRNAanalysis::pathanalysis(clusterShamEnterocyte,"Enterocytesham")
snRNAanalysis::pathwayplot("pathway_analysis/Enterocytesham_hallmark.txt","Enterocyte_Sham","Enterocyte")


iris %>% group_by(Species) %>% 
  summarise(length=max(Sepal.Length),
            Sepal.Width=Sepal.Width[which.max(Sepal.Length)])

