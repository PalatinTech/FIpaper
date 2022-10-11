#Psuedo bulk RNA -seq
library(SingleCellExperiment)
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(snRNAanalysis)

#normalization function
normobj<-function(input)
{
  input <- SCTransform(input,  verbose = FALSE)
  # These are now standard steps in the Seurat workflow for visualization and clustering
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(input), 20)

  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(input)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2

  input<- RunPCA(input, verbose = FALSE)
  DimHeatmap(input, dims = 1:30, cells = 500, balanced = TRUE)

  ElbowPlot(input,ndims = 40)
  # NOTE: This process can take a long time for big datasets, comment out for expediency. More


  input <- RunUMAP(input, dims = 1:20, verbose = FALSE) #used 20 dim for PL9654
  input<-RunTSNE(input,dims=1:20,verbose=FALSE)

  input<- FindNeighbors(input, dims = 1:20, verbose = FALSE)
  input <- FindClusters(input, verbose = FALSE,resolution=0.6)

  return(input)

}


Sham<-readRDS("/mnt2/rat_colon_scdata_UC/codes/FI_paper/Rdata/sham.RData")
Placebo<-readRDS("/mnt2/rat_colon_scdata_UC/codes/FI_paper/Rdata/placebo.RData")

Pl8177_50ug<-readRDS("/mnt2/rat_colon_scdata_UC/codes/FI_paper/Rdata/PL8177_50.RData")
#PL8177_50<-readRDS("/mnt2/rat_colon_scdata_UC/round2/codes/Rdata/PL")

#combine all the three treatments
PL8177_sham_placebo<-merge(Pl8177_50ug,y=c(Sham,Placebo),
                       add.cell.ids=c("PL8177 50ug","Sham","Placebo"))

#Extract raw counts and meta data
counts<-PL8177_sham_placebo@assays$RNA@counts
metadata<-PL8177_sham_placebo@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(PL8177_sham_placebo@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "sample")]
groups <- colData(sce)[, c( "sample")]
groups

# Explore the raw counts for the dataset

## Check the assays present
assays(sce)

## Explore the raw counts for the dataset
dim(counts(sce))

counts(sce)[1:6, 1:6]

## Explore the cellular metadata for the dataset
dim(colData(sce))
head(colData(sce))

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei<-as.data.frame(PL8177_sham_placebo@meta.data %>%
  group_by(PL8177_sham_placebo@meta.data$type,PL8177_sham_placebo@meta.data$sample)%>%
  summarise(n()))
colnames(ei)<-c("group_id","sample_id","n_Cells")

ei%>%group_by(group_id)%>%summarise(sum(n_Cells))
dim(sce)
#remove genes which are expressed in less than 3 cells 
sce<-sce[ncol(counts(sce))-rowCounts(counts(sce),value=0) >=300,] 

dim(sce)

groups <- colData(sce)[, c( "sample","type")]
# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)),
                       groupings = groups, fun = "sum")
#now we have genes which are expressed in more than 3 cells and cells with minimium 200 genes (Seurat limit)
dim(pb)

pb[1:6, 1:6]

#Use DESEq2 to normalize the count data to account for difference in library size and RNA composition between samples
#Use normalize counts ot make some QC plots at the gene and sample level
#find DE expressed genes

#create count matrix for DE object
count_matrix<-data.frame(t(pb))


#create metadata for the count matrix
meta<-unique(data.frame(metadata$sample,metadata$type))
colnames(meta)<-c("sample","type")
rownames(meta)<-paste(meta$sample,"_",meta$type,sep="")
rownames(meta)<-str_replace(rownames(meta),"-",".")
rownames(meta)<-str_replace(rownames(meta)," ",".")
rownames(meta)<-str_replace(rownames(meta),",",".")
meta<-meta[order(meta$type,decreasing = TRUE),]

count_matrix<-count_matrix[,match(rownames(meta),colnames(count_matrix))]

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(meta) == colnames(count_matrix))

#create DESEq2 object
dds <- DESeqDataSetFromMatrix(count_matrix,
                              colData = meta,
                              design = ~ type)
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA

p<-DESeq2::plotPCA(rld, intgroup = "type",ntop=2000)


ggplot(data=p$data,aes(x=PC1,y=PC2,col=type,label=name))+
  geom_point(size=5)+geom_text_repel(size=6)+
  scale_color_manual(values=c("#58409A","red","#339966"))+
  theme(panel.background = element_blank(),
  axis.title = element_text(size=16),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(size=16,colour = "black"),
  legend.text = element_text(size=14),
  legend.position = "bottom")+
  xlab("PC1: 41% variance")+ylab("PC2: 19% variance")
