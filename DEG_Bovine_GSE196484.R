###########
### DEG ###
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(RUVSeq)
library(dplyr)

# Define a function to read and process data
read_and_process <- function(file, pattern = "_[^_]*_[^_]*$") {
  ct <- length(read.table(file, sep = "\t", nrows = 1))
  data <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE,
                     colClasses = c("character", rep("numeric", ct - 1)))
  colnames(data) <- gsub(pattern, "", colnames(data))
  return(data)
}

# File names
files <- c("Morula_merged_gene_name_expression.txt",
           "Blastocyst_merged_gene_name_expression.txt",
           "MII_merged_gene_name_expression.txt",
           "8Cells_merged_gene_name_expression.txt",
           "2Cells_merged_gene_name_expression.txt")

# Read and process data
cts_list <- lapply(files, read_and_process)

# Remove 'length' columns if present
cts_list <- lapply(cts_list, function(df) df[, !grepl("^length$", colnames(df))])

# Merge data
merged_data <- Reduce(function(x, y) merge(x, y, by = "name", all = TRUE), cts_list)

# Organize the final data
geneName <- merged_data[, 1]
merged_data <- merged_data[, -1]
rownames(merged_data) <- geneName

# Check sample names and remove duplicates
colnames(merged_data) <- make.names(colnames(merged_data), unique = TRUE)

# Create coldata dataframe
sample_groups <- c("Morula", "Blastocyst", "MII", "8Cells", "2Cells")
sample_counts <- c(length(cts_list[[1]]) - 1, length(cts_list[[2]]) - 1, length(cts_list[[3]]) - 1, length(cts_list[[4]]) - 1, length(cts_list[[5]]) - 1)
coldata <- data.frame(sample = colnames(merged_data),
                      group = rep(sample_groups, sample_counts))

# Ensure row names of coldata are sample names
rownames(coldata) <- coldata$sample

# Ensure coldata and merged_data are in the same order
coldata <- coldata[match(colnames(merged_data), rownames(coldata)), ]

# Verify matching order
stopifnot(all(rownames(coldata) == colnames(merged_data)))


# DESeq2 Analysis

x <- factor(coldata$group)
design <- model.matrix(~0 + x)
colnames(design) <- levels(x)
dds <- DESeqDataSetFromMatrix(countData = merged_data, colData = coldata, design = design)

# Filter lowly expressed genes
cpm <- apply(merged_data, 2, function(x) x / sum(x) * 10^6)
filter <- cpm > 1
keep <- rowSums(filter) >= (ncol(cpm) / 2)
dds <- dds[keep,]
cpm <- cpm[keep,]

# Generate RPKM file
#rpkm <- apply(merged_data, 2, function(x) x / sum(x)) * 10^9 / geneLength

# PCA and heatmap
rld <- vst(dds, blind = FALSE)
pheatmap(cor(assay(rld)), annotation_col = coldata, show_colnames = FALSE, main = "Correlation within samples")

pcaData <- plotPCA(rld, intgroup = c("group"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 1) +
  ggtitle("PCA using DESeq2") +
  geom_text_repel(aes(label = name, fontface = 2), size = 4, box.padding = 0.5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(plot.title = element_text(size = 14, family = "Helvetica", face = "bold", hjust = 0.5),
        text = element_text(size = 12, family = "Helvetica"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))


# Batch correction with RUVSeq

set <- newSeqExpressionSet(as.matrix(merged_data[keep, ]), phenoData = data.frame(x, row.names = colnames(merged_data)))
y <- DGEList(counts = counts(set), group = x)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type = "deviance")
seqUQ <- betweenLaneNormalization(set, which = "upper")
ruv <- RUVr(seqUQ, rownames(set), k = 1, res)
ddsruv <- DESeqDataSetFromMatrix(countData = normCounts(ruv), colData = coldata, design = design)
ddsruv <- DESeq(ddsruv)
cpm2 <- apply(counts(ddsruv), 2, function(x) x / sum(x) * 10^6)

# Heatmap post batch correction
rldruv <- vst(ddsruv, blind = FALSE)
pheatmap(cor(assay(rldruv)), annotation_col = coldata, show_colnames = FALSE, main = "Correlation within samples K=1")

# PCA post batch correction
pcaData <- plotPCA(rldruv, intgroup = c("group"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  ggtitle("PCA using DESeq2 K=1") +
  geom_text_repel(aes(label = name, fontface = 2), size = 4, box.padding = 0.5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(plot.title = element_text(size = 14, family = "Helvetica", face = "bold", hjust = 0.5),
        text = element_text(size = 12, family = "Helvetica"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))


ddsruv <- DESeq(ddsruv)  # Since correction K=1
#plotDispEsts(ddsruv)
resultsNames(ddsruv)
res1 <- results(ddsruv, contrast=list("X2Cells","X8Cells")) #2Cell vs 8CEll 
res2 <- results(ddsruv, contrast=list("X2Cells","MII")) #2Cell vs MII
res3 <- results(ddsruv, contrast=list("X2Cells","Blastocyst")) #2Cell vs Blastocyst
res4 <- results(ddsruv, contrast=list("X2Cells","Morula")) #2Cell vs Morula
res5 <- results(ddsruv, contrast=list("X8Cells","MII")) #8Cell vs MII
res6 <- results(ddsruv, contrast=list("X8Cells","Blastocyst")) #8Cell vs Blastocyst
res7 <- results(ddsruv, contrast=list("X8Cells","Morula")) #8Cell vs Morula
res8 <- results(ddsruv, contrast=list("Blastocyst","Morula")) #Blastocyst vs Morula
res9 <- results(ddsruv, contrast=list("Blastocyst","MII")) #Blastocyst vs MII
res10 <- results(ddsruv, contrast=list("Morula","MII")) #Morula vs MII


res1$padj <- ifelse(is.na(res1$padj), 1, res1$padj)
res2$padj <- ifelse(is.na(res2$padj), 1, res2$padj)
res3$padj <- ifelse(is.na(res3$padj), 1, res3$padj)
res4$padj <- ifelse(is.na(res4$padj), 1, res4$padj)
res5$padj <- ifelse(is.na(res5$padj), 1, res5$padj)
res6$padj <- ifelse(is.na(res6$padj), 1, res6$padj)
res7$padj <- ifelse(is.na(res7$padj), 1, res7$padj)
res8$padj <- ifelse(is.na(res8$padj), 1, res8$padj)
res9$padj <- ifelse(is.na(res9$padj), 1, res9$padj)
res10$padj <- ifelse(is.na(res10$padj), 1, res10$padj)


res1 <- cbind(as.data.frame(res1),cpm2)
res1f <- res1[res1$padj<0.05 & abs(res1$log2FoldChange)>1,]#5135
dim(res2f)
write.table(res1f, "DEG_2cell_8cell.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

res2 <- cbind(as.data.frame(res2),cpm2)
res2f <- res2[res2$padj<0.05 & abs(res2$log2FoldChange)>1,]#1680
write.table(res2f, "DEG_2cell_MII.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

res3 <- cbind(as.data.frame(res3),cpm2)
res3f <- res3[res3$padj<0.05 & abs(res3$log2FoldChange)>1,]#7870
write.table(res3f, "DEG_2cell_blastocyst.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

res4 <- cbind(as.data.frame(res4),cpm2)
res4f <- res4[res4$padj<0.05 & abs(res4$log2FoldChange)>1,]#7588
write.table(res4f, "DEG_2cell_morula.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

res5 <- cbind(as.data.frame(res5),cpm2)
res5f <- res5[res5$padj<0.05 & abs(res5$log2FoldChange)>1,]#5094
write.table(res5f, "DEG_8cell_MII.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

res6 <- cbind(as.data.frame(res6),cpm2)
res6f <- res6[res6$padj<0.05 & abs(res6$log2FoldChange)>1,]#6343
write.table(res6f, "DEG_8cell_blastocyst.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

res7 <- cbind(as.data.frame(res7),cpm2)
res7f <- res7[res7$padj<0.05 & abs(res7$log2FoldChange)>1,]#5416
write.table(res7f, "DEG_8cell_morula.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

res8 <- cbind(as.data.frame(res8),cpm2)
res8f <- res8[res8$padj<0.05 & abs(res8$log2FoldChange)>1,]#2503
write.table(res8f, "DEG_blastocyst_morula.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

res9 <- cbind(as.data.frame(res9),cpm2)
res9f <- res9[res9$padj<0.05 & abs(res9$log2FoldChange)>1,]#7568
write.table(res9f, "DEG_blastocyst_MII.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

res10 <- cbind(as.data.frame(res10),cpm2)
res10f <- res10[res10$padj<0.05 & abs(res10$log2FoldChange)>1,]#7115
write.table(res10f, "DEG_morula_MII.tsv", sep = "\t", row.names = TRUE, quote = FALSE)
dim(res10f)


res1f_UP <- res1f[res1f$log2FoldChange >1,]#2990
res1f_DOWN <- res1f[res1f$log2FoldChange <1,]#2145
dim(res1f_UP)
res2f_UP <- res2f[res2f$log2FoldChange >1,]#882
res2f_DOWN <- res2f[res2f$log2FoldChange <1,]#798
res3f_UP <- res3f[res3f$log2FoldChange >1,]#5407
res3f_DOWN <- res3f[res3f$log2FoldChange <1,]#2463
res4f_UP <- res4f[res4f$log2FoldChange >1,]#5267
res4f_DOWN <- res4f[res4f$log2FoldChange <1,]#2321
res5f_UP <- res5f[res5f$log2FoldChange >1,]#2399
res5f_DOWN <- res5f[res5f$log2FoldChange <1,]#2695
res6f_UP <- res6f[res6f$log2FoldChange >1,]#4713
res6f_DOWN <- res6f[res6f$log2FoldChange <1,]#1630
res7f_UP <- res7f[res7f$log2FoldChange >1,]#4178
res7f_DOWN <- res7f[res7f$log2FoldChange <1,]#1238
res8f_UP <- res8f[res8f$log2FoldChange >1,]#971
res8f_DOWN <- res8f[res8f$log2FoldChange <1,]#1532
res9f_UP <- res9f[res9f$log2FoldChange >1,]#2536
res9f_DOWN <- res9f[res9f$log2FoldChange <1,]#5032
res10f_UP <- res10f[res10f$log2FoldChange >1,]#2282
res10f_DOWN <- res10f[res10f$log2FoldChange <1,]#4833
dim(res10f_DOWN)


# volcano plot between 2Cell vs 8CEll####
category1=ifelse(res1$padj>0.05|abs(res1$log2FoldChange)<1,"Regularly expressed",ifelse(res1$log2FoldChange>1,"Highly expressed in 2Cell","Highly expressed in 8Cell"))
index=as.factor(category1)
if(length(levels(index))==3) { index=factor(category1,levels=c("Regularly expressed","Highly expressed in 2Cell","Highly expressed in 8Cell"),ordered=T) }
plot=data.frame(logFC=res1$log2FoldChange,log10p=-log10(res1$padj),index=index)
#pdf("volcano_plot_2Cellvs8Cell.pdf",width=8,height=5)
print(ggplot(plot %>% arrange(index),aes(x=logFC,y=log10p,col=index))+
        geom_point(aes(size=index))+
        ggtitle("Volcano plot of 2Cell vs 8Cell")+
        theme_bw()+theme_classic()+
        scale_colour_manual(values=c("#999999","#E69F00","#56B4E9"))+
        scale_size_manual(values=c(1,2,2))+
        scale_y_continuous(name=expression(bold(-log[10]("Adjusted p-value"))))+
        geom_vline(xintercept=1,colour="black",linetype=2)+
        geom_vline(xintercept=-1,colour="black",linetype=2)+
        geom_hline(yintercept=-log10(0.05),colour="black",linetype=2)+
        theme(plot.title=element_text(size=14,family="Helvetica",face="bold",hjust=0.2),
              text=element_text(size=12,family="Helvetica"),
              axis.title=element_text(face="bold"),
              axis.text.x=element_text(size=10,face="bold"),
              axis.text.y=element_text(size=10,face="bold"),
              legend.text=element_text(size=10,face="bold"),
              legend.title=element_text(size=10,face="bold")))
#dev.off()

# volcano plot between 2Cell vs MII####
category1=ifelse(res2$padj>0.05|abs(res2$log2FoldChange)<1,"Regularly expressed",ifelse(res2$log2FoldChange>1,"Highly expressed in 2Cell","Highly expressed in MII"))
index=as.factor(category1)
if(length(levels(index))==3) { index=factor(category1,levels=c("Regularly expressed","Highly expressed in 2Cell","Highly expressed in MII"),ordered=T) }
plot=data.frame(logFC=res2$log2FoldChange,log10p=-log10(res2$padj),index=index)
#pdf("volcano_plot_Anapc4_stl745_vs_WT.pdf",width=8,height=5)
print(ggplot(plot %>% arrange(index),aes(x=logFC,y=log10p,col=index))+
        geom_point(aes(size=index))+
        ggtitle("Volcano plot of 2Cell vs MII")+
        theme_bw()+theme_classic()+
        scale_colour_manual(values=c("#999999","#E69F00","#56B4E9"))+
        scale_size_manual(values=c(1,2,2))+
        scale_y_continuous(name=expression(bold(-log[10]("Adjusted p-value"))))+
        geom_vline(xintercept=1,colour="black",linetype=2)+
        geom_vline(xintercept=-1,colour="black",linetype=2)+
        geom_hline(yintercept=-log10(0.05),colour="black",linetype=2)+
        theme(plot.title=element_text(size=14,family="Helvetica",face="bold",hjust=0.2),
              text=element_text(size=12,family="Helvetica"),
              axis.title=element_text(face="bold"),
              axis.text.x=element_text(size=10,face="bold"),
              axis.text.y=element_text(size=10,face="bold"),
              legend.text=element_text(size=10,face="bold"),
              legend.title=element_text(size=10,face="bold")))
#dev.off()

# volcano plot between 2Cell vs Blastocyst####

category1=ifelse(res3$padj>0.05|abs(res3$log2FoldChange)<1,"Regularly expressed",ifelse(res3$log2FoldChange>1,"Highly expressed in 2Cell","Highly expressed in Blastocyst"))
index=as.factor(category1)
if(length(levels(index))==3) { index=factor(category1,levels=c("Regularly expressed","Highly expressed in 2Cell","Highly expressed in Blastocyst"),ordered=T) }
plot=data.frame(logFC=res3$log2FoldChange,log10p=-log10(res3$padj),index=index)
#pdf("volcano_plot_Dscaml1_stl740_vs_WT.pdf",width=8,height=5)
print(ggplot(plot %>% arrange(index),aes(x=logFC,y=log10p,col=index))+
        geom_point(aes(size=index))+
        ggtitle("Volcano plot of 2Cell vs Blastocyst")+
        theme_bw()+theme_classic()+
        scale_colour_manual(values=c("#999999","#E69F00","#56B4E9"))+
        scale_size_manual(values=c(1,2,2))+
        scale_y_continuous(name=expression(bold(-log[10]("Adjusted p-value"))))+
        geom_vline(xintercept=1,colour="black",linetype=2)+
        geom_vline(xintercept=-1,colour="black",linetype=2)+
        geom_hline(yintercept=-log10(0.05),colour="black",linetype=2)+
        theme(plot.title=element_text(size=14,family="Helvetica",face="bold",hjust=0.2),
              text=element_text(size=12,family="Helvetica"),
              axis.title=element_text(face="bold"),
              axis.text.x=element_text(size=10,face="bold"),
              axis.text.y=element_text(size=10,face="bold"),
              legend.text=element_text(size=10,face="bold"),
              legend.title=element_text(size=10,face="bold")))
#dev.off()

# volcano plot between 2Cell vs Morula####

category1=ifelse(res4$padj>0.05|abs(res4$log2FoldChange)<1,"Regularly expressed",ifelse(res4$log2FoldChange>1,"Highly expressed in 2Cell","Highly expressed in Morula"))
index=as.factor(category1)
if(length(levels(index))==3) { index=factor(category1,levels=c("Regularly expressed","Highly expressed in 2Cell","Highly expressed in Morula"),ordered=T) }
plot=data.frame(logFC=res4$log2FoldChange,log10p=-log10(res4$padj),index=index)
#pdf("volcano_plot_Dscaml1_stl740_vs_WT.pdf",width=8,height=5)
print(ggplot(plot %>% arrange(index),aes(x=logFC,y=log10p,col=index))+
        geom_point(aes(size=index))+
        ggtitle("Volcano plot of 2Cell vs Morula")+
        theme_bw()+theme_classic()+
        scale_colour_manual(values=c("#999999","#E69F00","#56B4E9"))+
        scale_size_manual(values=c(1,2,2))+
        scale_y_continuous(name=expression(bold(-log[10]("Adjusted p-value"))))+
        geom_vline(xintercept=1,colour="black",linetype=2)+
        geom_vline(xintercept=-1,colour="black",linetype=2)+
        geom_hline(yintercept=-log10(0.05),colour="black",linetype=2)+
        theme(plot.title=element_text(size=14,family="Helvetica",face="bold",hjust=0.2),
              text=element_text(size=12,family="Helvetica"),
              axis.title=element_text(face="bold"),
              axis.text.x=element_text(size=10,face="bold"),
              axis.text.y=element_text(size=10,face="bold"),
              legend.text=element_text(size=10,face="bold"),
              legend.title=element_text(size=10,face="bold")))
#dev.off()

# volcano plot between 8Cell vs MII####

category1=ifelse(res5$padj>0.05|abs(res5$log2FoldChange)<1,"Regularly expressed",ifelse(res5$log2FoldChange>1,"Highly expressed in 8Cell","Highly expressed in MII"))
index=as.factor(category1)
if(length(levels(index))==3) { index=factor(category1,levels=c("Regularly expressed","Highly expressed in 8Cell","Highly expressed in MII"),ordered=T) }
plot=data.frame(logFC=res5$log2FoldChange,log10p=-log10(res5$padj),index=index)
#pdf("volcano_plot_Dscaml1_stl740_vs_WT.pdf",width=8,height=5)
print(ggplot(plot %>% arrange(index),aes(x=logFC,y=log10p,col=index))+
        geom_point(aes(size=index))+
        ggtitle("Volcano plot of 8Cell vs MII")+
        theme_bw()+theme_classic()+
        scale_colour_manual(values=c("#999999","#E69F00","#56B4E9"))+
        scale_size_manual(values=c(1,2,2))+
        scale_y_continuous(name=expression(bold(-log[10]("Adjusted p-value"))))+
        geom_vline(xintercept=1,colour="black",linetype=2)+
        geom_vline(xintercept=-1,colour="black",linetype=2)+
        geom_hline(yintercept=-log10(0.05),colour="black",linetype=2)+
        theme(plot.title=element_text(size=14,family="Helvetica",face="bold",hjust=0.2),
              text=element_text(size=12,family="Helvetica"),
              axis.title=element_text(face="bold"),
              axis.text.x=element_text(size=10,face="bold"),
              axis.text.y=element_text(size=10,face="bold"),
              legend.text=element_text(size=10,face="bold"),
              legend.title=element_text(size=10,face="bold")))
#dev.off()

# volcano plot between 8Cell vs Blastocyst####

category1=ifelse(res6$padj>0.05|abs(res6$log2FoldChange)<1,"Regularly expressed",ifelse(res6$log2FoldChange>1,"Highly expressed in 8Cell","Highly expressed in Blastocyst"))
index=as.factor(category1)
if(length(levels(index))==3) { index=factor(category1,levels=c("Regularly expressed","Highly expressed in 8Cell","Highly expressed in Blastocyst"),ordered=T) }
plot=data.frame(logFC=res6$log2FoldChange,log10p=-log10(res6$padj),index=index)
#pdf("volcano_plot_Dscaml1_stl740_vs_WT.pdf",width=8,height=5)
print(ggplot(plot %>% arrange(index),aes(x=logFC,y=log10p,col=index))+
        geom_point(aes(size=index))+
        ggtitle("Volcano plot of 8Cell vs Blastocyst")+
        theme_bw()+theme_classic()+
        scale_colour_manual(values=c("#999999","#E69F00","#56B4E9"))+
        scale_size_manual(values=c(1,2,2))+
        scale_y_continuous(name=expression(bold(-log[10]("Adjusted p-value"))))+
        geom_vline(xintercept=1,colour="black",linetype=2)+
        geom_vline(xintercept=-1,colour="black",linetype=2)+
        geom_hline(yintercept=-log10(0.05),colour="black",linetype=2)+
        theme(plot.title=element_text(size=14,family="Helvetica",face="bold",hjust=0.2),
              text=element_text(size=12,family="Helvetica"),
              axis.title=element_text(face="bold"),
              axis.text.x=element_text(size=10,face="bold"),
              axis.text.y=element_text(size=10,face="bold"),
              legend.text=element_text(size=10,face="bold"),
              legend.title=element_text(size=10,face="bold")))
#dev.off()

# volcano plot between #8Cell vs Morula####

category1=ifelse(res7$padj>0.05|abs(res7$log2FoldChange)<1,"Regularly expressed",ifelse(res7$log2FoldChange>1,"Highly expressed in 8Cell","Highly expressed in Morula"))
index=as.factor(category1)
if(length(levels(index))==3) { index=factor(category1,levels=c("Regularly expressed","Highly expressed in 8Cell","Highly expressed in Morula"),ordered=T) }
plot=data.frame(logFC=res7$log2FoldChange,log10p=-log10(res7$padj),index=index)
#pdf("volcano_plot_Dscaml1_stl740_vs_WT.pdf",width=8,height=5)
print(ggplot(plot %>% arrange(index),aes(x=logFC,y=log10p,col=index))+
        geom_point(aes(size=index))+
        ggtitle("Volcano plot of 8Cell vs Morula")+
        theme_bw()+theme_classic()+
        scale_colour_manual(values=c("#999999","#E69F00","#56B4E9"))+
        scale_size_manual(values=c(1,2,2))+
        scale_y_continuous(name=expression(bold(-log[10]("Adjusted p-value"))))+
        geom_vline(xintercept=1,colour="black",linetype=2)+
        geom_vline(xintercept=-1,colour="black",linetype=2)+
        geom_hline(yintercept=-log10(0.05),colour="black",linetype=2)+
        theme(plot.title=element_text(size=14,family="Helvetica",face="bold",hjust=0.2),
              text=element_text(size=12,family="Helvetica"),
              axis.title=element_text(face="bold"),
              axis.text.x=element_text(size=10,face="bold"),
              axis.text.y=element_text(size=10,face="bold"),
              legend.text=element_text(size=10,face="bold"),
              legend.title=element_text(size=10,face="bold")))
#dev.off()

# volcano plot between Blastocyst vs Morula####

category1=ifelse(res8$padj>0.05|abs(res8$log2FoldChange)<1,"Regularly expressed",ifelse(res8$log2FoldChange>1,"Highly expressed in Blastocyst","Highly expressed in Morula"))
index=as.factor(category1)
if(length(levels(index))==3) { index=factor(category1,levels=c("Regularly expressed","Highly expressed in Blastocyst","Highly expressed in Morula"),ordered=T) }
plot=data.frame(logFC=res8$log2FoldChange,log10p=-log10(res8$padj),index=index)
#pdf("volcano_plot_Dscaml1_stl740_vs_WT.pdf",width=8,height=5)
print(ggplot(plot %>% arrange(index),aes(x=logFC,y=log10p,col=index))+
        geom_point(aes(size=index))+
        ggtitle("Volcano plot of Blastocyst vs Morula")+
        theme_bw()+theme_classic()+
        scale_colour_manual(values=c("#999999","#E69F00","#56B4E9"))+
        scale_size_manual(values=c(1,2,2))+
        scale_y_continuous(name=expression(bold(-log[10]("Adjusted p-value"))))+
        geom_vline(xintercept=1,colour="black",linetype=2)+
        geom_vline(xintercept=-1,colour="black",linetype=2)+
        geom_hline(yintercept=-log10(0.05),colour="black",linetype=2)+
        theme(plot.title=element_text(size=14,family="Helvetica",face="bold",hjust=0.2),
              text=element_text(size=12,family="Helvetica"),
              axis.title=element_text(face="bold"),
              axis.text.x=element_text(size=10,face="bold"),
              axis.text.y=element_text(size=10,face="bold"),
              legend.text=element_text(size=10,face="bold"),
              legend.title=element_text(size=10,face="bold")))
#dev.off()

# volcano plot between Blastocyst vs MII####

category1=ifelse(res9$padj>0.05|abs(res9$log2FoldChange)<1,"Regularly expressed",ifelse(res9$log2FoldChange>1,"Highly expressed in Blastocyst","Highly expressed in MII"))
index=as.factor(category1)
if(length(levels(index))==3) { index=factor(category1,levels=c("Regularly expressed","Highly expressed in Blastocyst","Highly expressed in MII"),ordered=T) }
plot=data.frame(logFC=res9$log2FoldChange,log10p=-log10(res9$padj),index=index)
#pdf("volcano_plot_Dscaml1_stl740_vs_WT.pdf",width=8,height=5)
print(ggplot(plot %>% arrange(index),aes(x=logFC,y=log10p,col=index))+
        geom_point(aes(size=index))+
        ggtitle("Volcano plot of Blastocyst vs MII")+
        theme_bw()+theme_classic()+
        scale_colour_manual(values=c("#999999","#E69F00","#56B4E9"))+
        scale_size_manual(values=c(1,2,2))+
        scale_y_continuous(name=expression(bold(-log[10]("Adjusted p-value"))))+
        geom_vline(xintercept=1,colour="black",linetype=2)+
        geom_vline(xintercept=-1,colour="black",linetype=2)+
        geom_hline(yintercept=-log10(0.05),colour="black",linetype=2)+
        theme(plot.title=element_text(size=14,family="Helvetica",face="bold",hjust=0.2),
              text=element_text(size=12,family="Helvetica"),
              axis.title=element_text(face="bold"),
              axis.text.x=element_text(size=10,face="bold"),
              axis.text.y=element_text(size=10,face="bold"),
              legend.text=element_text(size=10,face="bold"),
              legend.title=element_text(size=10,face="bold")))
#dev.off()

# volcano plot between Morula vs MII####

category1=ifelse(res10$padj>0.05|abs(res10$log2FoldChange)<1,"Regularly expressed",ifelse(res10$log2FoldChange>1,"Highly expressed in Morula","Highly expressed in MII"))
index=as.factor(category1)
if(length(levels(index))==3) { index=factor(category1,levels=c("Regularly expressed","Highly expressed in Morula","Highly expressed in MII"),ordered=T) }
plot=data.frame(logFC=res10$log2FoldChange,log10p=-log10(res10$padj),index=index)
#pdf("volcano_plot_Dscaml1_stl740_vs_WT.pdf",width=8,height=5)
print(ggplot(plot %>% arrange(index),aes(x=logFC,y=log10p,col=index))+
        geom_point(aes(size=index))+
        ggtitle("Volcano plot of Morula vs MII")+
        theme_bw()+theme_classic()+
        scale_colour_manual(values=c("#999999","#E69F00","#56B4E9"))+
        scale_size_manual(values=c(1,2,2))+
        scale_y_continuous(name=expression(bold(-log[10]("Adjusted p-value"))))+
        geom_vline(xintercept=1,colour="black",linetype=2)+
        geom_vline(xintercept=-1,colour="black",linetype=2)+
        geom_hline(yintercept=-log10(0.05),colour="black",linetype=2)+
        theme(plot.title=element_text(size=14,family="Helvetica",face="bold",hjust=0.2),
              text=element_text(size=12,family="Helvetica"),
              axis.title=element_text(face="bold"),
              axis.text.x=element_text(size=10,face="bold"),
              axis.text.y=element_text(size=10,face="bold"),
              legend.text=element_text(size=10,face="bold"),
              legend.title=element_text(size=10,face="bold")))
#dev.off()

