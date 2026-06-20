# test R script
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

data <- readRDS("/opt/omicsdata/datasets/SRP047194.rds")


# Split each sample into fields
x <- strsplit(data$sra.sample_attributes, "\\|")

# Function to extract values in fixed order
extract <- function(v) {
  sapply(v, function(y) strsplit(y, ";;")[[1]][2])
}

# Apply to all samples
mat <- t(sapply(x, extract))

# Convert to dataframe
df <- as.data.frame(mat, stringsAsFactors = FALSE)

# Add column names manually
colnames(df) <- c("cell_type", "donor_id", "days", "source_name")

head(df)

df$days <- as.numeric(df$days)
df$ASD <- grepl("ASD", df$source_name)


