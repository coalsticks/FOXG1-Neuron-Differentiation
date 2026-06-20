library(tidyr)
library(tidyverse)
library(readr)
library(ggplot2)
library(omicsdata)
library(limma)
library(EDASeq)
library(edgeR)
library(dplyr)
library(plyr)
library(corrplot)
library(ggfortify)
library(genefilter)
library(ggrepel)
library(pheatmap)
library(cluster)
library(knitr)
library(DESeq2)
library(clusterProfiler)
library(pathview)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggnewscale)
library(ggupset)
library(ggridges)

data <- readRDS("/opt/omicsdata/datasets/SRP047194.rds")


coldata_list <- lapply(data$sra.sample_attributes, function(x) {
  kv <- strsplit(strsplit(x, "\\|")[[1]], ";;")
  keys <- sapply(kv, `[`, 1)
  values <- sapply(kv, `[`, 2)
  df <- as.data.frame(t(values), stringsAsFactors = FALSE)
  colnames(df) <- keys
  df
})
# coldata_list$sample <- colnames(data)
coldata_df <- do.call(rbind, coldata_list)
coldata_df <- separate_wider_delim(coldata_df, source_name, delim = ",", names = c("cell-type", "DaIV", "Control"))
coldata_df <- mutate(coldata_df, DIV = `number of days in vitro in terminal differentiation conditions`)
coldata_df <- coldata_df |>
  dplyr::select(-`cell type`, -`cell-type`, -`DaIV`, -`number of days in vitro in terminal differentiation conditions`)
glimpse(coldata_df)
coldata_df$sample <- colnames(data)
coldata_df <- mutate_at(coldata_df, "DIV", as.integer) |>
  filter(DIV != 0)


table(rowSums(assay(data, "raw_counts")) >= 20)
data <- data[rowSums(assay(data, "raw_counts")) > 20, ]
data <- data[, colnames(data) %in% coldata_df$sample]
  
data$sra.sample_attributes <- coldata_df
glimpse(data$sra.sample_attributes)

assay(data, "raw_counts")[1:5,1:5]
# glimpse(colnames(data))

raw_count <- assay(data, "raw_counts")
coldata_df <- mutate(coldata_df, factorDIV = factor(DIV))
col <- as.factor(coldata_df$Control)
levels(col) <- c("red","blue")
col <- as.character(col)

boxplot(log1p(raw_count), col=col, xlab = "Sample", ylab = "log1p expression", ylim = c(0, 20),  las = 2, cex.axis = 0.7)
limma::plotMA(log1p(raw_count), main = "MA_graph", xlab = "A", ylab = "M")
MDPlot(raw_count, c(1, 7))
plotRLE(raw_count, outline = FALSE, las = 2,
        ylab = "RLE", main = "RLE of raw counts", col = col)

pca <- prcomp(t(log1p(raw_count)))
summary(pca)$importance
autoplot(pca)
autoplot(pca, data = coldata_df, colour = "Control", main = "Condition colour")
autoplot(pca, data = coldata_df, colour = "DIV", main = "DIV colour")

t11 <- as.data.frame(pca$x[,c(1,2)])
t11$subtype <- col
ggplot(t11, aes(PC1,PC2, color=subtype)) + geom_point()
# dfpca <- as.data.frame(pca$x)
# ggplot(dfpca, aes(PC1,PC3, color=col)) + geom_point()

ensid <- as.character(rowRanges(data)$gene_id)
ensid <- substr(ensid, 1, 15)

# View(data$gene_id)

gc <- getGeneLengthAndGCContent(id = ensid, org = "hg38", mode = "org.db")
table(is.na(gc[, 2]))
gc <- gc[!is.na(gc[, 2]), ]
rowData(data)$ensid <- ensid
data <- data[rowData(data)$ensid  %in% rownames(gc), ]
raw_count <- assay(data, "raw_counts")
row.names(raw_count) <- row.names(gc)

raw_count

# table(duplicated(row.names(raw_count)))

set <- newSeqExpressionSet(raw_count,
                           phenoData = AnnotatedDataFrame(data.frame(conditions = factor(1:ncol(raw_count)),
                                                                     row.names = colnames(raw_count))),
                           featureData = AnnotatedDataFrame(data.frame(gc = gc[, 2], l = gc[, 1])))
set

before <- EDASeq::biasPlot(set, "gc", ylim = c(0, 10), log = TRUE)
lrt <- log(raw_count[, 1] + 1) - log(raw_count[, 13] + 1) # ASK WHY/WHAT WE NEED TO SEE CHANGE
biasBoxplot(lrt, gc[, 2], outline = FALSE, xlab = "GC-content", ylab = "Gene expression", las = 2, cex.axis = 0.5) 

# 

wt <- withinLaneNormalization(set, "gc" , which = "upper", offset = TRUE)
wt
after <- biasPlot(wt, "gc", ylim = c(0, 10), log = TRUE)

set_full <- betweenLaneNormalization(set, which = "full", offset = TRUE)
set_uq <- betweenLaneNormalization(set, which = "upper", offset = TRUE)
full <- set_full@assayData$normalizedCounts
uq_edaseq <- set_uq@assayData$normalizedCounts

# raw_count <- counts(set)

lbs <- calcNormFactors(raw_count, method = "RLE")
rle <- cpm(raw_count, lib.size = lbs * colSums(raw_count))

lbs <- calcNormFactors(raw_count, method = "TMM")
tmm <- cpm(raw_count, lib.size = lbs * colSums(raw_count))

lbs <- calcNormFactors(raw_count, method = "upperquartile")
uq <- cpm(raw_count, lib.size = lbs * colSums(raw_count))

tc <- cpm(raw_count)


boxplot(log2(raw_count + 1), outline = FALSE,  las=2,
        ylab = "log(count+1)", main = "Before normalization", xaxt = 'n')

boxplot(log2(tc + 1), outline = FALSE, las = 2,
        ylab = "log(count+1)", main = "Boxplot Total count", xaxt='n')

boxplot(log2(rle + 1), outline=FALSE,  las=2,
        ylab = "log(count+1)", main = "Boxplot RLE", xaxt = 'n')

boxplot(log2(tmm + 1), outline=FALSE,  las=2,
        ylab = "log(count+1)", main = "Boxplot TMM", xaxt = 'n')

boxplot(log2(uq + 1), outline=FALSE,  las=2,
        ylab = "log(count+1)", main = "Boxplot Upper quartile", xaxt = 'n')

boxplot(log2(full + 1), outline=FALSE,  las=2,
        ylab = "log(count+1)", main = "Boxplot full", xaxt = 'n')


EDASeq::plotRLE(raw_count, outline = FALSE,  las=2, col=col,
                ylab = "RLE", main = "RLE before normalization", xaxt='n')

EDASeq::plotRLE(tc, outline = FALSE,  las=2, col=col,
                ylab = "RLE", main = "RLE total count", xaxt='n')

EDASeq::plotRLE(tmm, outline = FALSE,  las=2, col=col,
                ylab = "RLE", main = "RLE TMM", xaxt='n')

EDASeq::plotRLE(rle, outline = FALSE,  las=2, col=col,
                ylab = "RLE", main = "RLE RLE", xaxt='n')

EDASeq::plotRLE(uq, outline = FALSE,  las=2, col=col,
                ylab = "RLE", main = "RLE upper quartile", xaxt='n')

EDASeq::plotRLE(full, outline = FALSE,  las=2, col=col,
                ylab = "RLE", main = "RLE full", xaxt='n')

colnames <- colnames(tmm)
colnames(tmm) <- substr(colnames(tmm), 9, 12)
m <- cor(log1p(tmm))

corrplot(m, method = "color", order = 'hclust', tl.cex=0.2)

data_pca <- prcomp(t(log1p(tmm)))
tmp <- summary(data_pca)
tmp$importance

autoplot(data_pca)
screeplot(data_pca, type = c("lines"))

# Review 'cause it doesn't show points
EDASeq::plotPCA(tc, k=3, labels=F, col=as.numeric(col), pch=20)
EDASeq::plotPCA(tmm, k=3, labels=F, col=as.numeric(col), pch=20)
EDASeq::plotPCA(rle, k=3, labels=F, col=as.numeric(col), pch=20)

table(data$sra.sample_attributes[,2:3])
data_analysis <- data[, data$sra.sample_attributes[,3] != 0]
# View(data_analysis$sra.sample_attributes)
data_analysis$sra.sample_attributes <- dplyr::select(data_analysis$sra.sample_attributes, DIV)
data_analysis$sra.sample_attributes <- as.character(as.factor(unlist(data_analysis$sra.sample_attributes)))
# data_analysis$sra.sample_attributes <- as.character(data_analysis$sra.sample_attributes)
# table(data_analysis$sra.sample_attributes)
# data$sra.sample_attributes[,3]
design <- model.matrix(~ sra.sample_attributes, data = colData(data_analysis))
# head(design)
dge <- DGEList(assay(data_analysis, "raw_counts"))
rowRanges(data_analysis)$gene_id

colSums(assay(data_analysis, "raw_counts"))
colData(data_analysis)
dge <- calcNormFactors(dge, method="upperquartile") 
dge <- estimateDisp(dge, design) 

plotMeanVar(dge, show.raw.vars = TRUE, show.ave.raw.vars = FALSE, NBline = TRUE, show.tagwise.vars = TRUE)

fit <- glmFit(dge, design)
fit
res <- glmLRT(fit, coef=2)
topTags(res, sort.by="logFC")

tab <- topTags(res, n=Inf)
tab <- tab$table
head(tab)

de <- tab[tab$FDR <0.05,]
dim(de)

plotMD(res)

names <- row.names(de)

tab$sign <- NA
tab$sign[tab$logFC < 0] <- "Down"
tab$sign[tab$logFC > 0] <- "Up"
tab$sign[tab$FDR >= 0.05] <- "NoSig"
tab$gene <- rownames(tab)

ggplot(tab, aes(x = logFC, y = -log10(PValue), color = sign)) + geom_point() + theme_classic() +
  scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NoSig" = "black")) +
  geom_text_repel(data = tab[1:20, ], aes(label = gene), size = 2.5, max.overlaps = 20)

hist(tab$PValue, main="p-values distribution")

raw_count <- assay(data_analysis, "raw_counts")
lbs <- calcNormFactors(raw_count, method = "TMM")
tmm <- cpm(raw_count, lib.size = lbs * colSums(raw_count))

luad_degs <- tmm[names, ]

Col <- data.frame(Subtype=data_analysis$sra.sample_attributes)
rownames(Col) <- colnames(data_analysis)

ann_colors = list(Subtype = c("11"="blue","31"="lightblue"))

g <- pheatmap(log1p(luad_degs), show_colnames = FALSE, scale = "row", show_rownames = FALSE,
              drop_levels = TRUE, annotation_colors = ann_colors, annotation_col = Col)

cut <- cutree(g$tree_col, k = 2)

table(cut)

dist_mat <- dist(t(log1p(luad_degs)))
plot(silhouette(cut, dist_mat), col=c("#2E86C1", "#EC7063"), border=NA)

top <- topTags(res, n = Inf)$table
deg <- top[top$FDR <= 0.05,]
dim(deg)

universo <- substr(row.names(top), 1, 15)
sign <- substr(row.names(deg), 1, 15)

universo
sign

universo.entrez <- mapIds(org.Hs.eg.db, # genome wide annotation for Human that contains Entrez Gene identifiers
                          keys = universo, #Column containing Ensembl gene ids
                          column = "ENTREZID",
                          keytype = "ENSEMBL")  

sign.entrez <- mapIds(org.Hs.eg.db,
                      keys = sign, #Column containing Ensembl gene ids
                      column = "ENTREZID",
                      keytype = "ENSEMBL")
universo.entrez <- as.character(universo.entrez)
sign.entrez <- as.character(sign.entrez)

kk <- enrichKEGG(gene = sign.entrez,
                 universe = universo.entrez,
                 organism = 'hsa',
                 keyType = 'ncbi-geneid', #ENTREZID
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
dim(kk)
View(kk@result[1:10,])

barplot(kk) #Quickly rank terms by the number of genes involved
dotplot(kk) #Useful when you want to distinguish between a term that is significant but contains few of your genes vs. one where your DEGs dominate the pathway.
cnetplot(kk) # Useful to see gene overlap between terms and identify multi-pathway genes
heatplot(kk, showCategory = 10) #Reveals which genes are shared across multiple terms and which are unique to one pathway

edo <- pairwise_termsim(kk)
emapplot(edo) #Enrichment Map, terms are nodes; edges connect terms that share genes
obj <- upsetplot(edo, 8)
obj

pathview(gene.data = sign.entrez, pathway.id = "hsa04110", species = "hsa")
knitr::include_graphics("hsa04110.pathview.png")

kk <- enrichPathway(gene = sign.entrez,
                    universe = universo.entrez,
                    organism = 'human',
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH")

View(kk@result[1:10, 1:6])

barplot(kk)
dotplot(kk)
cnetplot(kk)
heatplot(kk, showCategory = 10)

edo <- pairwise_termsim(kk)
emapplot(edo)
upsetplot(edo, 5)

ego <- enrichGO(gene = sign.entrez,
                universe = universo.entrez,
                OrgDb = 'org.Hs.eg.db',
                keyType = "ENTREZID",
                ont = "BP",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                readable = TRUE)

View(ego@result[1:10,])

barplot(ego)
upsetplot(ego, 5)

tt <- top$logFC # logFC relies on the mean but does not take into account the variability of data -> always better to use the statistical value of a test
names(tt) <- universo.entrez

table(duplicated(names(tt)))
tt <- tt[!duplicated(names(tt))]

sort_tt <- sort(tt, decreasing = TRUE)

kk2 <- gseKEGG(geneList = sort_tt,
               organism = "hsa",
               pvalueCutoff = 0.05)

View(kk2@result[, 1:7])

kkGO <- gseGO(geneList = sort_tt,
              ont = "BP",
              OrgDb='org.Hs.eg.db',
              keyType = "ENTREZID",
              pAdjustMethod = "BH")

View(kkGO@result[1:10,])
