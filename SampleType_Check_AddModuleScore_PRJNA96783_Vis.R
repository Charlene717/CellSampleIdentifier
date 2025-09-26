###############################################################################
## 0) 套件準備
###############################################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(pheatmap)
})

###############################################################################
## 1) 定義基因集合
###############################################################################
genesets <- list(
  KF_MF          = c("POSTN","COL11A1","COL12A1","COL5A2","COMP","ADAM12","TNFRSF12A","SDC1",
                     "ITGBL1","FN1","THBS2","ASPN","TAGLN","ACTA2"),
  IL6_STAT3      = c("IL6","IL6R","STAT3","SOCS3","CCL2","CXCL1","CXCL2","CXCL3","CXCL8",
                     "NFKBIA","WNT5A","PTGS2","OSMR"),
  SPF            = c("APCDD1","ID1","WIF1","COL13A1","COL18A1","COL23A1","PI16","CXCL14","DPP4"),
  SRF            = c("WISP2","SLPI","TSPAN8","COL18A1","SFRP2","PDGFRB"),
  FB_NormalPlus  = c("COL14A1","PCOLCE2","MFAP5","OGN","CFD"),
  KC_Homeo       = c("KRT14","KRT5","DST","DSG3","ITGA6","ITGB1","COL17A1","PLEC",
                     "LAMA3","LAMB3","LAMB2","EPCAM","CLDN1",
                     "KLF4","OVOL1","GRHL1","GRHL3"),
  Endo_Quiescent = c("KDR","FLT1","ESAM","PECAM1","CDH5","CLDN5","RBP7","GPIHBP1","PROX1","KLF2","KLF4"),
  TGFB_targets   = c("CTGF","SERPINE1","THBS1","COMP","TAGLN","ACTA2","COL1A1","COL1A2","COL3A1","FN1"),
  Tx_Up          = c("MMP1","MMP3","DUSP1","KLF2","KLF6","ZFP36","TSC22D3","FKBP5"),
  KF_KC_Up       = c("S100A8", "S100A9", "KRT6A", "KRT6B", "KRT16", "KRT17", "SPRR1A", "SPRR1B", "SPRR2A", 
                     "SPRR2D", "PI3", "ELAFIN", "LCN2", "DEFB4A", "CXCL1", "CXCL8", "IL8", "CCL20", "IL1B", 
                     "IL36G", "TNFAIP3", "NFKBIA", "AREG", "HBEGF", "TGFA", "EREG", "MMP9", "MMP10", "MMP13", 
                     "PLAU", "SERPINB2"),
  KF_KC_Down      =c("FLG", "FLG2", "LOR", "KRT1", "KRT10", "KRT2", "DSG1", "DSC1", "CASP14", 
                     "TGM1", "SPINK5", "KLK5", "KLK7", "SBSN", "ALOX12B", "ALOXE3", "PNPLA1", "CDSN")
)

###############################################################################
## 2) 計算 Module Scores
###############################################################################
seuratObject_Sample <- AddModuleScore(
  object   = seuratObject_Sample,
  features = genesets,
  name     = names(genesets)   # 每組基因集會產生一個 *_1 欄位
)

###############################################################################
## 3) 小提琴圖（每個 module score 依 sample_id）
###############################################################################
plot_list <- list()
for (nm in names(genesets)) {
  # 找出 AddModuleScore 生成的對應欄位（可能是 nm1 或 nm_1）
  colname <- grep(paste0("^", nm), colnames(seuratObject_Sample@meta.data), value = TRUE)
  if (length(colname) == 0) next  # 如果沒找到就跳過
  
  p <- VlnPlot(
    seuratObject_Sample,
    features = colname,
    group.by = "sample_id",
    pt.size  = 0
  ) +
    labs(title = nm, y = "Module score") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  plot_list[[nm]] <- p
}

# 拼接成大圖
combined_plot <- wrap_plots(plot_list, ncol = 4)
print(combined_plot)

###############################################################################
## 4) Heatmap：顯示每個 sample_id 的平均 module score
###############################################################################
# 抽取 meta.data 的 module score 欄位
score_cols <- grep("_1$", colnames(seuratObject_Sample@meta.data), value = TRUE)

# 計算每個 sample_id 的平均值
df_scores <- seuratObject_Sample@meta.data %>%
  select(all_of(c("sample_id", score_cols))) %>%
  group_by(sample_id) %>%
  summarise(across(all_of(score_cols), mean, na.rm = TRUE)) %>%
  as.data.frame()

rownames(df_scores) <- df_scores$sample_id
df_scores <- df_scores[ , -1, drop = FALSE]


summary(df_scores)

# 移除全為 NA 的欄位
df_scores <- df_scores[, colSums(is.na(df_scores)) < nrow(df_scores), drop = FALSE]

# 把剩餘欄位的 NA 用 0 補上（避免 scale 出錯）
df_scores[is.na(df_scores)] <- 0
# 畫 heatmap
pheatmap(df_scores,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         scale        = "row",
         main         = "Average Module Scores per sample_id")
