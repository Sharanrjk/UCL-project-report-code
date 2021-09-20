setwd("C:/Users/Sharan/Desktop/root/UCL files/proj Melissa 25 nov/")
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)


mdata = readRDS("singlets_24Nov20.rds")

mdata = readRDS("mdata_final2.rds")

mdata = readRDS(file="CiteSeq_24Nov20.rds") # This includes doublets and negatives as well

# Analyzing metadata
mdata[["percent.mt"]] = PercentageFeatureSet(mdata, pattern = "^Mt-")

VlnPlot(mdata, features = c("nFeature_RNA", "nCount_RNA", "nCount_HTO","percent.mt"), ncol = 4)
hist(mdata@meta.data$HTO_classification)
plot(mdata@meta.data$HTO_classification )
which(mdata$HTO_classification == "HTO1") %>% length()
which(mdata$HTO_classification.global == "Singlet") %>% length()


levels(mdata$HTO_classification)
mdata$HTO_classification = factor(mdata$HTO_classification)
levels(mdata$HTO_classification) = c("TICH-1","TICH-2","LPS","Serum")

FeatureScatter(mdata, feature1 = "nCount_RNA" , feature2 ="nFeature_RNA")
plot(mdata$HTO_classification, mdata$nCount_RNA, main = "nCount distribution / HTO")
plot(mdata$HTO_classification, mdata$nFeature_RNA, main = "nFeature distribution / HTO")
plot(mdata$HTO_classification, mdata$percent.mt, main = "mito% distribution / HTO")
plot(mdata$HTO_classification, mdata$nCount_HTO, main = "nCount_HTO distribution / HTO")

VlnPlot(mdata, features = c("Itgam") ,split.by = "HTO_classification", 
                pt.size = 0.1, cols = c("blue","red")) +aes(color=mdata$HTO_classification)


##################
mdata = subset(mdata, subset = HTO_classification != c("HTO3"))   # Run this if you want to subset metadata and then analyze it
##################

mdata <- NormalizeData(mdata, normalization.method = "LogNormalize", scale.factor = 10000)
mdata <- FindVariableFeatures(mdata, nfeatures = 4000, )
top10 <- head(VariableFeatures(mdata), 10)

plot1 <- VariableFeaturePlot(mdata, cols = c('black','black'))
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
plot1 

features_table = mdata@assays[["RNA"]]@meta.features # The variable features are stored here  as vst. values
detected_genes = features_table[features_table$vst.mean > 0,] %>% row.names()

# Scaling data
mdata <- ScaleData(mdata)

# sctransform used sometimes? check it

# Linear dimensional reduction
mdata <- RunPCA(mdata, npcs = 25, ndims.print = 1:5, nfeatures.print = 5)
ElbowPlot(mdata, ndims = 40)

# Examine and visualize PCA results a few different ways
print(mdata[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mdata, dims = 1:2, reduction = "pca")
DimPlot(mdata, reduction = "pca", dims = c(3,5))
DimHeatmap(mdata, dims = c(1:10), cells = 500, balanced = TRUE)

FeaturePlot(mdata, reduction = "tsne", features = "nFeature_RNA", pt.size = 1.5)


#Analyzing cell cycle (genes correspond to human genes?)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(stringr)
s.genes = tolower(s.genes) %>% str_to_title()
g2m.genes = tolower(g2m.genes) %>% str_to_title()

mdata = CellCycleScoring(mdata, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(mdata, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)

DimPlot(mdata,
        reduction = "pca", pt.size = 1.5,
        split.by= "HTO_classification",
        dims = c(1,5))

mdata$Phase = as.factor(mdata$Phase)

plot(mdata$HTO_classification, mdata$Phase)
levels(mdata$Phase)
FeaturePlot(mdata, features = "Apoe", reduction = "pca")

# Add module score for certain cell types:
microglia_gene.list = c("Tmem119","Sall1","Slc2a5","Itgam", "Aif1",'P2ry12','Cx3cr1')
macrophage_gene.list = c("Ms4a3","Ccr2","Lyz2","Aif1","Itgal","Cd14")
DAM_stage1.list = c("Ctsb","Ctsd","Apoe","B2m","Fth1","Lyz2","Tyrobp","Timp2")
DAM_stage2.list = c( "Axl", "Cst7","Ctsl", "Lpl", "Cd9", "Clec7a" , "Csf1" , "Ccl6" , "Itgax" , "Spp1",'Cd63','Cadm1','Cd68','Ctsz','Ctsa','Cd52','Ccl6','Ank','Serpine2','Gusb','Hif1a')
ARM_gene.list = c("Cst7",'Clec7a','Itgax', 'Cd74','H2-ab1','H2-aa','Ctsb','Ctsd','Spp1','Gpnmb','Dkk2')
senescence_gene.list = c("Cdkn1a ","Cdkn2a", "B4galt1","Cdkn2b")
Homeostatic_gene.list = c('C1qa','C1qb','C1qc','Ctss','Hexb','Olfml3','Csf1r','Cst3')
cellcycle_gene.list = c('Top2a','Cdk1','Mcm6','Mki67')
Interferon_response_gene.list = c('Oas1a','Oasl2','Irf7','Ifit3','Mx1')
Bohlen_inflammation_gene.list = c("Tnf",'Ccl2','Il1b','Cxcl10','Cd14','Il6')
write(cellcycle_gene.list, 'cellcycle.txt')

gene.list = Homeostatic_gene.list
name = "Homeostatic_score"

mdata = AddModuleScore(
  object = mdata,
  features = list(gene.list),
  name = name
)

FeaturePlot(mdata, features = ARM_gene.list, pt.size = 1 ) 

VlnPlot(mdata, features = c('C1qa','C1qb','C1qc','Apoe','Srebf1','Scd2'), pt.size = 0 )

Idents(object = mdata) <- 'HTO_classification'

head(x = mdata[])
# subsetting cells for just G1, will have to run scaling etc again. RUN JUST FOR ANALYSIS- not needed anymore
##  mdata = subset(mdata, subset = Phase == "G1")


# Clustering
mdata <- FindNeighbors(mdata, reduction = "pca", dims = 1:20, nn.eps = 0.5)
mdata <- FindClusters(mdata, resolution = c(0.4, 0.6, 0.8, 1.0, 1.2), n.start = 10)

Idents(object = mdata) <- "monocle_cluster_highres"

#Cluster tree and clustree
mdata = BuildClusterTree(mdata, dims = 1:20)
tree = mdata@tools$BuildClusterTree
tree$tip.label <- paste0("Cluster ", tree$tip.label)

library(ggplot2)
p <- ggtree::ggtree(tree, aes(x, y)) +
  scale_y_reverse() +
  ggtree::geom_tree() +
  ggtree::theme_tree() +
  ggtree::geom_tiplab(offset = 1) +
  ggtree::geom_tippoint(shape = 16, size = 5) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,2.5,0,0), 'cm'))

p


library(clustree)
clustree(mdata, prefix = "RNA_snn_res.")

# for visualization
mdata <- RunTSNE(mdata, dims = 1:15, max_iter = 2000)
mdata <- RunUMAP(mdata, dims = 1:15, min.dist = 0.75)

p1 <- DimPlot(mdata, reduction = "tsne", 
              pt.size = 2, label = TRUE, 
              label.color = "red") + ggtitle(label = "tSNE ")

p2 <- DimPlot(mdata, reduction = "umap", 
              pt.size = 1, label = FALSE, label.color = "red",
              group.by = "monocle_cluster_highres") + ggtitle(label = "UMAP")

p1 <- AugmentPlot(plot = p1)
p2 <- AugmentPlot(plot = p2)
p1+p2  
p1
p2

p3 <- FeaturePlot(mdata, features = c("Cd14","Gch1","Mmp7","Fcer1g"), reduction = "umap",
                  pt.size = 1.5 ,combine = FALSE)
p3
wrap_plots(p3)



#Cluster analysis

cluster1.markers <- FindMarkers(mdata, ident.1 = 1, min.pct = 0.5 )
head(cluster1.markers, n = 10)

cluster2.markers <- FindMarkers(mdata, ident.1 = 2, min.pct = 0.5)
head(cluster2.markers, n = 10)

cluster3.markers <- FindMarkers(mdata, ident.1 = 3, min.pct = 0.5)
head(cluster3.markers, n = 10) 

cluster5.markers <- FindMarkers(mdata, ident.1 = 5, min.pct = 0.5)
head(cluster5.markers, n = 10)

cluster4.markers <- FindMarkers(mdata, ident.1 = 4, min.pct = 0.5)
head(cluster4.markers, n = 10)

cluster6.markers <- FindMarkers(mdata, ident.1 = 6, min.pct = 0.5)
head(cluster6.markers, n = 10)

cluster7.markers <- FindMarkers(mdata, ident.1 = 7, min.pct = 0.5)
head(cluster7.markers, n = 10)

cluster8.markers <- FindMarkers(mdata, ident.1 = 8, min.pct = 0.5)
head(cluster7.markers, n = 10)

cluster9.markers <- FindMarkers(mdata, ident.1 = 9, min.pct = 0.5)
head(cluster9.markers, n = 10)

cluster2to3.markers <- FindMarkers(mdata, ident.1 = 2, ident.2 = 3, min.pct = 0.5)
head(cluster2to3.markers, n = 10)

cluster3_5.markers <- FindMarkers(mdata, ident.1 = c(3,5), min.pct = 0.5)
head(cluster3_5.markers, n = 10)

probe = cluster3.markers
RidgePlot(mdata, features = row.names(head(probe, n = 3)), ncol = 3) 
VlnPlot(mdata, features = row.names(head(probe, n = 5)))

markers = c("Trem2","Apoe","Axl","Lpl","Clec7a","Spp1")
FeaturePlot(mdata, features = markers, reduction = "tsne")

# Volcano plot to represent DEG
res = cluster3.markers
## most DE 
res <- res[order(res$p_val_adj), ]
mostDE <- head(rownames(res),10)
res$mostDE <- rownames(res) %in% mostDE

## volcano plot
ggplot(res,
       aes(x = avg_logFC, y = -log10(p_val_adj), color = mostDE)) +
  geom_point() +
  ggtitle("Volcano plot") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value")



#Inter-HTO analysis (similar format to comparere cells based on metadata)
HTO3to4.markers <- FindMarkers(mdata, ident.1 = "HTO3", ident.2 = "HTO4", group.by = "HTO_classification" ,min.pct = 0.5)
head(HTO3to4.markers, n = 10) 

HTO_Serum_to_TICH.markers <- FindMarkers(mdata, ident.1 = "Serum", ident.2 = c("TICH-1","TICH-2"), 
                               group.by = "HTO_classification" ,min.pct = 0.1, test.use = "wilcox")
head(HTO3to1.markers, n = 20)

HTO4to1.markers <- FindMarkers(mdata, ident.1 = "HTO4", ident.2 = c("HTO1","HTO2"), 
                               group.by = "HTO_classification" ,min.pct = 0.5, test.use = "wilcox")
head(HTO4to1.markers, n = 20)

TICH2.markers <- FindMarkers(mdata, ident.1 = "TICH-2", 
                               group.by = "HTO_classification" ,min.pct = 0.5, test.use = "wilcox")

HTO3to4_conserved.markers <- FindConservedMarkers(mdata, ident.1 = "HTO3", ident.2 = 8, grouping.var = "HTO_classification" ,min.pct = 0.5)
# Not working, find another way

Idents(object = mdata) <- "HTO_classification"
diffexpgenes = head(HTO2to1.markers, n = 12) %>% row.names()

VlnPlot(mdata, features = c("Tnf","Ccl2","Il1b","Cxcl10","Cd14","Il6") , idents = c("TICH-1","TICH-2"), 
        pt.size = 0.1, cols = c("orange","green","blue","red"), ncol = 3)    

# find markers for every cluster compared to all remaining cells, report only the positive ones - then heatmap
mdata.markers <- FindAllMarkers(mdata, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)
mdata.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top10 <- mdata.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top5 <- mdata.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(mdata, features = top5$gene)


# inter metadata proprtion analysis. in this case plotting proportion fo HTOs in each cluster
dat <- data.frame(table(mdata$monocle_cluster_highres,mdata$HTO_classification))
names(dat) <- c("Clusters","HTO","Count")
ggplot(data=dat, aes(x=Clusters, y=Count, fill=HTO)) + geom_bar(stat="identity", position = "stack")  #position takes dodge or stack

# Gene ontology % enrichment analysis
library(gprofiler2)

probes = cluster4.markers[cluster4.markers$avg_logFC > 0,] %>% row.names()

gostres <- gost(query = probes, 
                organism = "rnorvegicus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = detected_genes, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
head(gostres$result, 10)

gostplot(gostres, capped = TRUE, interactive = TRUE)

# Fischers table

features_table = mdata@assays[["RNA"]]@meta.features # The variable features are stored here  as vst. values
detected_genes = features_table[features_table$vst.mean > 0,] %>% row.names()
C_detected_genes = readRDS("Carlo_detected_genes.rds")

ARM.markers = readRDS("ARM_markers.rds")
TRM.markers = readRDS("TRM_markers.rds")
IRM.markers = readRDS("IRM_markers.rds")
H1M.markers = readRDS("H1M_markers.rds")
H2M.markers = readRDS("H2M_markers.rds")
CC.markers = readRDS("proliferation_markers.rds")

Our_cluster = mdatapseudobulk_DESEQtable
Gene_markers = bohlen_DESEQtable


cluster_markers = row.names(Our_cluster[Our_cluster$avg_log_FC < 0.25 ,]) 
gene_list = row.names(Gene_markers[Gene_markers$avg_log_FC > 0.25 ,]) 

non_cluster_markers = detected_genes[ !detected_genes %in% cluster_markers ]
non_list_markers = C_detected_genes[ !C_detected_genes %in% gene_list]
  
  
A = (length(intersect(cluster_markers, gene_list)))
B = (length(intersect(cluster_markers, non_list_markers)))
C = (length(intersect(non_cluster_markers, gene_list)))
D = (length(intersect(non_cluster_markers, non_list_markers)))


fishers.table <- matrix(c(A,B,C,D),ncol=2,byrow=TRUE)
fishers.table <- as.table(fishers.table)
fishers.table
fisher.test(fishers.table)



chisq.test(fishers.table)


###### Monocle pipeline & Trajectory analysis ######################################################

library(SeuratWrappers)
library(monocle3)

mdata.cds = as.cell_data_set(mdata)  # Try this again without wrapper, extract expression data etc from scratch?
mdata.cds = preprocess_cds(mdata.cds)
mdata.cds = reduce_dimension(mdata.cds)
mdata.cds = reduce_dimension(mdata.cds, reduction_method = "tSNE")

mdata.cds = cluster_cells(mdata.cds, reduction_method = "UMAP")
mdata.cds = learn_graph(mdata.cds)

rowData(mdata.cds)$gene_short_name = mdata.cds@rowRanges@partitioning@NAMES
rowData(mdata.cds)$gene_short_name = rowData(mdata.cds)$gene_short_name  %>% as.factor()

plot_cells(mdata.cds, color_cells_by = "cluster",
           cell_size = 1, label_cell_groups = FALSE, reduction_method = "UMAP")



DimPlot(mdata, reduction = "umap", 
        pt.size = 2, label = TRUE, label.color = "red" )

plot_cells(mdata.cds, genes = "Tmem176a")

plot_cells(mdata.cds, color_cells_by = "HTO_classification", label_cell_groups = FALSE)





