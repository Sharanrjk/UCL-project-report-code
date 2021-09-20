setwd("C:/Users/Sharan/Desktop/root/UCL files/proj Andreone/")
library(dplyr)
library(DESeq2)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(limma)

andreone_data = read.csv("GSE144120_dst72_counts_edited.csv")
andreone_data = andreone_data[,-c(1,2, 4:12)]
row.names(andreone_data) = andreone_data[,1]
andreone_data = andreone_data[,-c(1,4,5)]




# Filtering data
any(is.na( andreone_data ))

rm = rowMeans(andreone_data, na.rm = TRUE)
length(rm)
length(rm[rm>0.01])
expressed_genes = which(rm>0.01)
andreone_data = andreone_data[expressed_genes,]

# selecting genes with CV > 5%
rm = rowMeans(andreone_data, na.rm = TRUE)
stdev = sapply(andreone_data, sd, na.rm = TRUE)
cv = stdev/rm
head(cv)
head(cv[cv>0.05])
length(cv[cv > 0.05])  # Very less genes with low co-variance hence not filtering them out

condition = c(rep("IPSC",2), "IPSC_Zymosan","IPSC_Zymosan", "Trem2_KO","Trem2_KO","Trem2_KO_Zymosan","Trem2_KO_Zymosan",
               "Plcg2_KO","Plcg2_KO","Plcg2_KO_Zymosan","Plcg2_KO_Zymosan") %>% as.factor()

condition_plcg = c(rep("non-plcg-KO",8), rep("plcg",4)) %>% as.factor()


batch = c(1,1,2,2,rep(1,10)) %>% as.factor()

columns_metadata = data.frame(condition, condition_plcg, row.names =  colnames(andreone_data))


#Combat for batch effects
library(sva)
adjusted_counts_combat <- ComBat_seq(andreone_data, batch=batch, group=NULL, full_mod=FALSE)

#Limma for batch effects
adjusted_counts_limma <- limma::removeBatchEffect(andreone_data, columns_metadata$batch)
limma::plotMDS(adjusted_counts_limma, labels = columns_metadata$batch, gene.selection = "common")


par(mfrow=c(1,1))
boxplot(as.data.frame(andreone_data),main="Original", ylim = c(0,200))
boxplot(as.data.frame(adjusted_counts_combat),main="Batch corrected combat", ylim=c(0,200))
boxplot(as.data.frame(adjusted_counts_limma),main="Batch corrected limma", ylim=c(0,200))



# creating DeSeq object
dds <- DESeqDataSetFromMatrix(countData =  andreone_data, #Make sure to specify correct input data & colData
                              colData = columns_metadata, 
                              design = ~  condition )


# Transform counts for data visualization & outlier detection (PCA & Heatmap)
rld <- rlog(dds, blind=TRUE)   #normalize and rlog transforms raw counts


mat <- assay(rld)
mat <- limma::removeBatchEffect(mat, rld$batch)
mat_combat <- ComBat_seq(mat, rld$batch)
assay(rld) <- mat

# Plot PCA

DESeq2::plotPCA(rld, intgroup = "condition")

pcaData <- plotPCA(rld, intgroup=c("condition", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap::pheatmap(rld_cor, annotation = columns_metadata[, c("condition"), drop=F])


#Another way to do it:

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition)
colnames(sampleDistMatrix) <- paste(rld$condition)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)




# Run DESeq2 differential expression analysis
dds <- DESeq(dds)
# Plot dispersion estimates
plotDispEsts(dds)

# Contrast
contrast_WT_plcgKO = c("condition","IPSC","Plcg2_KO")
contrast_WTzymosan_PlcgKOzymosan = c("condition","IPSC_Zymosan","Plcg2_KO_Zymosan")
contrast_WT_zymosan = c("condition", "IPSC","IPSC_Zymosan")
contrast_Plcg2_zymosan= c("condition","Plcg2_KO","Plcg2_KO_Zymosan")

res <- results(dds, 
               contrast = contrast_WT_plcgKO,
               alpha = 0.05)


res <- lfcShrink(dds, 
                 contrast =  contrast_WT_plcgKO,
                 res=res,
                 type = "normal")


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output
res_tbl

# Table of results for significant genes
# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res


#Scatterplot of normalized expression of top 20 most significant genes

normalized_counts <- counts(dds, 
                            normalized = TRUE)

## Order results by padj values
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)


top20_sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes)


gathered_top20_sig <- top20_sig_norm %>%
  gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")



ggplot(gathered_top20_sig) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = samplename, 
                 size = 0.5), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes for selected contrast group") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))

# Heatmap of all genes

# Extract normalized counts for only the significant genes
sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res$gene)

# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")
# Run pheatmap using the metadata data frame for the annotation
a <- pheatmap::pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
                   color = heat_colors, # Find a cool color palette
                   cluster_rows = T, 
                   show_rownames = F,
                   annotation = columns_metadata[, c("condition", "condition_plcg")], 
                   border_color = NA, 
                   fontsize = 10, 
                   scale = "row", 
                   fontsize_row = 10, 
                   height = 20,
                   main = "Heatmap of DEGs between selected contrast groups")



res_table_thres <- res_tbl %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.50)


## creating genelabels column for labelling
res_table_thres <- res_table_thres %>% arrange(padj) %>% mutate(genelabels = "")

# For specific genes, get position of them first then insert. 
interest_genes = (c("TREM2","Itgam", "CD14","Itgb1","Itgb2","Itgb3","Itgb5","Fn1","Rac2","Vav","Syk", "Btk","Blnk","Tyrobp") %>% toupper())
Andreone_genes = c("Il1b","Cxcl1","Cxcl8","Ccl2","Il15","Il10","Il6","Ccl20","Il23a","Tnfsf15","Il1rn","Cxcl5","Cxcl3","Cxcl10","Cxcl14")


interest_pos = ( match(interest_genes,res_table_thres$gene) %>% na.omit() )



res_table_thres$genelabels[c(1:20,27)] <- res_table_thres$gene[c(1:20, 27)]

## Volcano plot
ggplot(res_table_thres, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("Volcano plot of DEGs between Plcg2_KO and WT") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +       #(this line makes plot look cleaner but labels dont work if used)
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) + theme(panel.border = element_blank(),
                                                             panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank()) + coord_cartesian(ylim = c(0,50))



# MA plot; requires 3 columns: mean, log fold change, logical table like threshold or significance

ma = res_table_thres[, c("baseMean", "log2FoldChange", "threshold")]
plotMA(ma,ylim=c(-5,5))

# GO
library(gprofiler2)
Gprofiler <- gost(sig_res$gene[which(sig_res$log2FoldChange > 0)], organism = "mmusculus", ordered_query = FALSE,
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                  measure_underrepresentation = FALSE, evcodes = FALSE,
                  user_threshold = 0.05, correction_method = c("bonferroni"),#"g_SCS", "bonferroni","fdr", "false_discovery_rate", "gSCS", "analytical",
                  custom_bg = row.names(andreone_data),
                  numeric_ns = "", sources = NULL) #Gene set enrichment analysis


gostplot(Gprofiler, capped = FALSE, interactive = TRUE)
p <- gostplot(Gprofiler, capped = FALSE, interactive = FALSE) #Plot results

publish_gosttable(Gprofiler, highlight_terms = Gprofiler$result[c(1:2,10,100:102,120,124,125),],
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size"),
                  filename = NULL)

pp <- publish_gostplot(p, highlight_terms = c("GO:0048856","GO:0009653","GO:0006928","GO:0008092","GO:0040011"), 
                       width = NA, height = NA, filename = NULL )

gc <- gconvert(query = c("GO:0009653"), organism = "hsapiens", 
         target="ENSG", mthreshold = Inf, filter_na = TRUE)


GO_genes = intersect(sig_res$gene,gc$name)

GO_genes

