library(tidyr)
library(dplyr)
library(omicsdata)
library(limma)
library(EDASeq)
library(edgeR)
library(plyr)
library(corrplot)
library(ggfortify)
library(SummarizedExperiment)
library(tidyverse)


data <- readRDS("/opt/omicsdata/datasets/SRP047194.rds")																								

# raw_count_before_filter <- 

class(data)


# extract raw strings into data frame
raw <- data.frame(
  sample = colnames(data), 
  raw_attr = as.character(data$sra.sample_attributes)
)

# split key-value pairs for attributes and attribute values
tidy_attr <- raw %>%
  separate_rows(raw_attr, sep = "\\|") %>%
  separate(raw_attr, into = c("key", "value"),
           sep = ";;", fill = "right") %>%
  mutate(key = str_trim(key),
         value = str_trim(value))

# pivot wider so each attribute becomes its own column
# modifying last attribute to specify ASD or Non-ASD and remove redundant data
sample_meta <- tidy_attr %>%
  pivot_wider(names_from = key, values_from = value) %>%
  mutate(condition = if_else(str_detect(source_name, "non-ASD"), "non-ASD", "ASD")) %>%
  dplyr::select(-source_name, -`cell type`) %>%
  dplyr::rename(days_in_vitro = `number of days in vitro in terminal differentiation conditions`) %>%
  filter(days_in_vitro != 0)  # ---> CHECKPOINT CORRECT

# exploratory data analysis

# sanity check on sample_meta data set
glimpse(sample_meta)

# how many samples per condition
sample_meta <- sample_meta %>%
  mutate(days_in_vitro = as.numeric(days_in_vitro))

dplyr::count(sample_meta, days_in_vitro)
dplyr::count(sample_meta, condition)
dplyr::count(sample_meta, condition, days_in_vitro)


# filter raw count inside summarized experiment

data <- data[rowSums(assay(data, "raw_counts")) > 20, ]  

# only samples present in sample_meta
data <- data[, colnames(data) %in% sample_meta$sample]  # --> checkpoint 2


# -------------------------------------------CHECKPOINT-1-----------------------------------------------> CONFIRMED

# EXPLORATORY DATA ANALYSIS

# sample distribution
library(ggplot2)

ggplot(sample_meta, aes(x = condition, fill = condition)) +
  geom_bar() + 
  facet_wrap(~ days_in_vitro, labeller = label_both) + 
  labs(
    title = "Sample distribution by condition and timepoint",
    x = "Condition",
    y = "Number of sample"
  ) +
  theme_minimal()

# BOXPLOT

# log-transformed
col <- factor(sample_meta$condition)
levels(col) <- c("red", "blue")
col <- as.character(col)

raw_count <- assay(data, "raw_counts")

par(mfrow = c(1, 1))

boxplot(log1p(raw_count),
        main = "Log-transformed counts",
        ylab = "log1p counts",
        las = 2,
        cex.axis = 0.7,
        col = col)

legend("topright",
       legend = c("ASD", "non-ASD"),
       fill = c("red", "blue"))

# RLE
library(EDASeq)

plotRLE(raw_count,
        outline = FALSE,
        las = 2,
        ylab = "RLE",
        main = "RLE of raw counts",
        col = col,
        cex.axis = 0.5)

legend("topright",
       legend = c("ASD", "non-ASD"),
       fill = c("red", "blue"))

# RLE - check outlier samples
rle_matrix <- log1p(raw_count) - rowMedians(log1p(raw_count))
rle_medians <- colMedians(rle_matrix)

# find samples that deviate most from zero
names(rle_medians) <- colnames(raw_count)
sort(abs(rle_medians), decreasing = TRUE)


# PCA
library(ggfortify)

pca <- prcomp(t(log1p(raw_count)))
summary(pca)$importance

autoplot(pca)

# by condition -> confirm paper findings
autoplot(pca, data = sample_meta, colour = "condition", 
         main = "PCA colored by condition")

autoplot(pca, data = sample_meta, colour = "days_in_vitro",
         main = "PCA colored by timepoint")


# MA plot
limma::plotMA(log1p(raw_count),
              main = "MA plot",
              xlab = "A",
              ylab = "M")

# MD plot
MDPlot(raw_count, c(1, 7))

# check how many genes have all zero counts
table(rowSums(raw_count) == 0)

# check how many genes have very low counts
table(rowSums(raw_count) < 10)

table(rowSums(raw_count) < 20)

# filter out lowly expressed genes
raw_count <- raw_count[rowSums(raw_count) >= 20, ]

dim(raw_count)


# RLE before normalization
plotRLE(raw_count, outline = FALSE, las = 2,
        ylab = "RLE", main = "RLE of raw counts", col = col)


# GC content
ensid <- as.character(rowRanges(data)$gene_id)
ensid <- substr(ensid, 1, 15)

gc <- getGeneLengthAndGCContent(id = ensid, org = "hg38", mode = "org.db")
table(is.na(gc[, 2]))

gc <- gc[!is.na(gc[, 2]), ]
rowData(data)$ensid <- ensid
data <- data[rowData(data)$ensid %in% rownames(gc), ]
raw_count <- assay(data, "raw_counts")
row.names(raw_count) <- row.names(gc)
raw_count


# Normalization

set <- newSeqExpressionSet(raw_count,
                           phenoData = AnnotatedDataFrame(data.frame(conditions = factor(1:ncol(raw_count)),
                                                                     row.names = colnames(raw_count))),
                           featureData = AnnotatedDataFrame(data.frame(gc = gc[, 2], l = gc[, 1])))
set

# Before normalization

before <- EDASeq::biasPlot(set, "gc", ylim = c(0, 10), log = TRUE)
lrt <- log(raw_count[, 1] + 1) - log(raw_count[, 13] + 1)
biasBoxplot(lrt, gc[, 2], outline = FALSE, xlab = "GC-content", ylab = "Gene expression", las = 2, cex.axis = 0.5)


# withinlanenormalization
wt <- withinLaneNormalization(set, "gc", which = "upper", offset = TRUE)
wt

after <- biasPlot(wt, "gc", ylim = c(0, 10), log = TRUE)



# EDASeq between-sample
set_full <- betweenLaneNormalization(set, which = "full",  offset = TRUE)
set_uq   <- betweenLaneNormalization(set, which = "upper", offset = TRUE)
full      <- set_full@assayData$normalizedCounts
uq_edaseq <- set_uq@assayData$normalizedCounts


raw_count <- counts(set)

# edgeR

lbs <- calcNormFactors(raw_count, method = "RLE")
rle <- cpm(raw_count, lib.size = lbs * colSums(raw_count))

lbs <- calcNormFactors(raw_count, method = "TMM")
tmm <- cpm(raw_count, lib.size = lbs * colSums(raw_count))

lbs <- calcNormFactors(raw_count, method = "upperquartile")
uq  <- cpm(raw_count, lib.size = lbs * colSums(raw_count))

tc  <- cpm(raw_count)

# PLOTS EVALUATION

col <- factor(sample_meta$condition)
levels(col) <- c("red", "blue")
col <- as.character(col)
col_numeric <- as.numeric(factor(sample_meta$condition))

# --- Boxplots ---
boxplot(log2(raw_count + 1), outline = FALSE, las = 2,
        ylab = "log2(count+1)", main = "Before normalization", xaxt = "n")
boxplot(log2(tc   + 1), outline = FALSE, las = 2,
        ylab = "log2(count+1)", main = "Total count", xaxt = "n")
boxplot(log2(rle  + 1), outline = FALSE, las = 2,
        ylab = "log2(count+1)", main = "RLE", xaxt = "n")
boxplot(log2(tmm  + 1), outline = FALSE, las = 2,
        ylab = "log2(count+1)", main = "TMM", xaxt = "n")
boxplot(log2(uq   + 1), outline = FALSE, las = 2,
        ylab = "log2(count+1)", main = "Upper quartile", xaxt = "n")
boxplot(log2(full + 1), outline = FALSE, las = 2,
        ylab = "log2(count+1)", main = "Full quantile", xaxt = "n")

# RLE Plots
EDASeq::plotRLE(raw_count, outline = FALSE, las = 2, col = col,
                ylab = "RLE", main = "Before normalization", xaxt = "n")
EDASeq::plotRLE(tc,   outline = FALSE, las = 2, col = col,
                ylab = "RLE", main = "Total count",          xaxt = "n")
EDASeq::plotRLE(rle,  outline = FALSE, las = 2, col = col,
                ylab = "RLE", main = "RLE",                  xaxt = "n")
EDASeq::plotRLE(tmm,  outline = FALSE, las = 2, col = col,
                ylab = "RLE", main = "TMM",                  xaxt = "n")
EDASeq::plotRLE(uq,   outline = FALSE, las = 2, col = col,
                ylab = "RLE", main = "Upper quartile",       xaxt = "n")
EDASeq::plotRLE(full, outline = FALSE, las = 2, col = col,
                ylab = "RLE", main = "Full quantile",        xaxt = "n")

# Correlation
tmm_colnames <- colnames(tmm)
colnames(tmm) <- substr(colnames(tmm), 1, 8)
m <- cor(log1p(tmm))
corrplot(m, method = "color", order = "hclust", tl.cex = 0.2,
         main = "Correlation - TMM")
colnames(tmm) <- tmm_colnames

# PCA - after normalization
asd_pca <- prcomp(t(log1p(tmm)))
summary(asd_pca)$importance[, 1:10]

autoplot(asd_pca, data = sample_meta, colour = "condition",
         main = "PCA - TMM, coloured by condition")
screeplot(asd_pca, type = "lines", main = "Screeplot - TMM")

EDASeq::plotPCA(tc,   k = 3, labels = FALSE, col = col_numeric, pch = 20,
                main = "PCA - Total count")
EDASeq::plotPCA(tmm,  k = 3, labels = FALSE, col = col_numeric, pch = 20,
                main = "PCA - TMM")
EDASeq::plotPCA(rle,  k = 3, labels = FALSE, col = col_numeric, pch = 20,
                main = "PCA - RLE")
EDASeq::plotPCA(full, k = 3, labels = FALSE, col = col_numeric, pch = 20,
                main = "PCA - Full quantile")



