library(Seurat)
library(harmony)
library(Rcpp)
library(clustree)
library(ggraph)
library(RColorBrewer)
library(patchwork)
library(R.methodsS3)
library(R.utils)
library(tidyverse) 


setwd("/Users/wangjingyu/desktop/aging/scRNA")
list.files('./')
getwd()

rm(list = ls());gc()

set.seed(123) 

sce.all<-readRDS('scRNA.rds')

#
sce.all=PercentageFeatureSet(sce.all,'^mt-',col.name = "percent.mt")
sce.all<-subset(sce.all,subset=nFeature_RNA>200 & percent.mt<5)
#
test.seu<-sce.all
test.seu<-test.seu%>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData()
test.seu <- RunPCA(test.seu, npcs = 50, verbose = FALSE)

test.seu=test.seu %>% RunHarmony("orig.ident", plot_convergence = TRUE)
test.seu<-RunUMAP(test.seu, dim=1:30,
                  reduction ="harmony")
DimPlot(test.seu,reduction = "umap")

sce_cluster=test.seu

sce_cluster<-FindNeighbors(sce_cluster,reduction="harmony",
                           dims=1:30)
scRNAseq<-sce_cluster
scRNAseq <- FindClusters(scRNAseq, resolution = 0.5)
DimPlot(scRNAseq, reduction = "umap")

markers_genes<-FindAllMarkers(scRNAseq,
                              logfc.threshold = 0.25,
                              #test.use="wilcox",
                              min.pct = 0.25,
                              #min.diff.pct=0.2,
)
write.csv(markers_genes,"markers_genes.csv")

bfreaname.scRNAseq <- scRNAseq
new.cluster.ids <- c("Astrocyte","Microglia","Oligodendrocyte","Oligodendrocyte","Fibroblast","Fibroblast","Mixed cells","Fibroblast", "EC","OPC","SMC",'Mixed cells',"Ependymal cells","Neuron")
names(new.cluster.ids) <- levels(scRNAseq)
scRNAseq <- RenameIdents(scRNAseq, new.cluster.ids)
DimPlot(scRNAseq,reduction = "umap", pt.size = 0.01) 

DimPlot(scRNAseq, split.by = "orig.ident", pt.size = 0.01)

top5<-markers_genes%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
DotPlot(scRNAseq,features=unique(top5$gene),assay="RNA")+RotatedAxis()+ggplot2:::coord_flip()

VlnPlot(scRNAseq, features = c("FTH1"),pt.size=0,flip = T,add.noise = T,split.by = 'orig.ident',split.plot = T)+
  
  theme(axis.text.y = element_blank(),
        
        axis.ticks.y = element_blank(),
        
        axis.title = element_blank(),
        
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        
        legend.position = 'none')

scRNAseq$celltype <- Idents(scRNAseq)
scRNAseq$celltype.orig.ident <- paste(scRNAseq$celltype, scRNAseq$orig.ident, sep = "_")
scRNAseq$celltype <- Idents(scRNAseq)
Idents(scRNAseq) <- "celltype.orig.ident"
I1VSS1 <- FindMarkers(scRNAseq, ident.1 = 'Astrocyte_OF1', ident.2 = 'Astrocyte_YF1', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
write.csv(I1VSS1,"Astrocytedeg.csv")

I1VSS1 <- FindMarkers(scRNAseq, ident.1 = 'Microglia_OF1', ident.2 = 'Microglia_YF1', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
write.csv(I1VSS1,"Microgliadeg.csv")

I1VSS1 <- FindMarkers(scRNAseq, ident.1 = 'Oligodendrocyte_OF1', ident.2 = 'Oligodendrocyte_YF1', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
write.csv(I1VSS1,"Oligodendrocytedeg.csv")

I1VSS1 <- FindMarkers(scRNAseq, ident.1 = 'Fibroblast_OF1', ident.2 = 'Fibroblast_YF1', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
write.csv(I1VSS1,"Fibroblastdeg.csv")

I1VSS1 <- FindMarkers(scRNAseq, ident.1 = 'EC_OF1', ident.2 = 'EC_YF1', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
write.csv(I1VSS1,"ECdeg.csv")

I1VSS1 <- FindMarkers(scRNAseq, ident.1 = 'OPC_OF1', ident.2 = 'OPC_YF1', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
write.csv(I1VSS1,"OPCdeg.csv")

I1VSS1 <- FindMarkers(scRNAseq, ident.1 = 'SMC_OF1', ident.2 = 'SMC_YF1', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
write.csv(I1VSS1,"SMCdeg.csv")

I1VSS1 <- FindMarkers(scRNAseq, ident.1 = 'Ependymal cells_OF1', ident.2 = 'Ependymal cells_YF1', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
write.csv(I1VSS1,"Ependymal cellsdeg.csv")

I1VSS1 <- FindMarkers(scRNAseq, ident.1 = 'Neuron_OF1', ident.2 = 'Neuron_YF1', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
write.csv(I1VSS1,"Neurondeg.csv")

table(scRNAseq$celltype.orig.ident)











