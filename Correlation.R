###############################################################################
## 計算 Sample 間相似度 (考慮 Cell_Type, 使用 log1p CPM)
###############################################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(pheatmap)
})

#----------------------------
# Step 1: 準備資料
#----------------------------
meta   <- seuratObject_Sample@meta.data
counts <- GetAssayData(seuratObject_Sample, slot = "counts")

#----------------------------
# Step 2: 建立 pseudo-bulk matrix
#   - 聚合到 (sample_id × Cell_Type)
#----------------------------
pseudo_bulk <- t(apply(counts, 1, function(gene_expr){
  tapply(gene_expr, paste(meta$sample_id, meta$Cell_Type, sep="__"), sum)
}))

# 移除全零基因
pseudo_bulk <- pseudo_bulk[rowSums(pseudo_bulk) > 0, ]

#----------------------------
# Step 3: 計算 CPM + log1p
#----------------------------
calc_cpm <- function(counts) {
  libsize <- colSums(counts)
  cpm <- t(t(counts) / libsize * 1e6)
  return(cpm)
}
log_expr <- log1p(calc_cpm(pseudo_bulk))

#----------------------------
# Step 4: 對每個 sample_id 建立 profile
#   - 把該 sample 的不同 Cell_Type profile 取平均
#   - 也可以改成挑特定 Cell_Type (只 subset column)
#----------------------------
sample_profiles <- sapply(unique(meta$sample_id), function(sid){
  cols <- grep(paste0("^", sid, "__"), colnames(log_expr))
  rowMeans(log_expr[, cols, drop=FALSE])
})

#----------------------------
# Step 5: 計算樣本間相似度 (Pearson correlation)
#----------------------------
sample_cor <- cor(sample_profiles, method = "pearson")

# 顯示結果矩陣
print(round(sample_cor, 3))

#----------------------------
# Step 6: 視覺化
#----------------------------
pheatmap(sample_cor,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Sample similarity (log1p CPM, considering Cell_Type)")
