#### =========================================================
#### Keloid / Normal / Post-Tx 分類（Sample-level）+ UMAP 輸出
#### =========================================================

suppressPackageStartupMessages({
  if (!require("Seurat")) install.packages("Seurat"); library(Seurat)
  if (!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
  if (!require("patchwork")) install.packages("patchwork"); library(patchwork)
})

# ---------- 0) 基本設定（依需求可調整） ----------
# 你的 Seurat 物件名稱請為 seuratObject_Sample，且 meta.data 內有 sample_id 欄位
# seuratObject_Sample@meta.data$sample_id <- seuratObject_Sample@meta.data$orig.ident

stopifnot(exists("seuratObject_Sample"))
if (!"sample_id" %in% colnames(seuratObject_Sample@meta.data)) {
  stop("seuratObject_Sample@meta.data 缺少 'sample_id' 欄位")
}

# [保留你原本的前處理為註解，不更動]
# DefaultAssay(seuratObject_Sample) <- if ("RNA" %in% names(seuratObject_Sample@assays)) "RNA" else DefaultAssay(seuratObject_Sample)
# if (!"data" %in% slotNames(GetAssay(seuratObject_Sample))) seuratObject_Sample <- NormalizeData(seuratObject_Sample, verbose = FALSE)
# if (!"scale.data" %in% slotNames(GetAssay(seuratObject_Sample))) {
#   seuratObject_Sample <- ScaleData(seuratObject_Sample, features = rownames(seuratObject_Sample), verbose = FALSE)
# }

# 期望數量的參數（可調整）
enforce_expected <- TRUE
expected_counts  <- c("Keloid" = 2L, "Normal" = 2L, "Post-Tx" = 2L)

# 輸出資料夾
out_dir <- "ModuleScore_Output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------- 1) 定義基因集合（微調 + 加入 Tx_Up） ----------
# Keloid / Fibrosis（含 mesenchymal CAF 標記）
genes_KF_MF <- c(
  "POSTN","COL11A1","COL12A1","COL5A2","COMP","ADAM12","TNFRSF12A","SDC1",
  "ITGBL1","FN1","THBS2","ASPN","TAGLN","ACTA2"
)

# IL6/JAK/STAT3 發炎（偏 PIF / 發炎軸）
genes_IL6_STAT3 <- c(
  "IL6","IL6R","STAT3","SOCS3","CCL2","CXCL1","CXCL2","CXCL3","CXCL8",
  "NFKBIA","WNT5A","PTGS2","OSMR"
)

# —— Normal：Dermal fibroblast（Papillary/Reticular/抗纖維化特徵）——
genes_SPF <- c("APCDD1","ID1","WIF1","COL13A1","COL18A1","COL23A1","PI16","CXCL14","DPP4")
genes_SRF <- c("WISP2","SLPI","TSPAN8","COL18A1","SFRP2","PDGFRB")
# Papillary/抗纖維化補強（文獻常見於正常真皮）：COL14A1、PCOLCE2、MFAP5、OGN、CFD
genes_FB_NormalPlus <- c("COL14A1","PCOLCE2","MFAP5","OGN","CFD")

# —— Normal：Keratinocyte homeostasis（基底層/屏障/分化）——
genes_KC_Homeo <- c(
  "KRT14","KRT5","DST","DSG3","ITGA6","ITGB1","COL17A1","PLEC",
  "LAMA3","LAMB3","LAMB2","EPCAM","CLDN1",
  "KLF4","OVOL1","GRHL1","GRHL3"
)

# —— Normal：Endothelial quiescence（靜息/流動應力應答）——
genes_Endo_Quiescent <- c(
  "KDR","FLT1","ESAM","PECAM1","CDH5","CLDN5","RBP7","GPIHBP1","PROX1","KLF2","KLF4"
)

# TGF-β 目標 / ECM（治療後通常下降）
genes_TGFB_targets <- c("CTGF","SERPINE1","THBS1","COMP","TAGLN","ACTA2","COL1A1","COL1A2","COL3A1","FN1")

# Tx_Up：治療後上升（降發炎/即時早期反應/ECM 降解）
genes_Tx_Up <- c("MMP1","MMP3","DUSP1","KLF2","KLF6","ZFP36","TSC22D3","FKBP5")

DefaultAssay(seuratObject_Sample) <- if ("RNA" %in% names(seuratObject_Sample@assays)) "RNA" else DefaultAssay(seuratObject_Sample)
.keep_present <- function(g) unique(intersect(g, rownames(seuratObject_Sample)))

signatures <- list(
  KF_MF            = .keep_present(genes_KF_MF),
  IL6_STAT3        = .keep_present(genes_IL6_STAT3),
  SPF              = .keep_present(genes_SPF),
  SRF              = .keep_present(genes_SRF),
  FB_NormalPlus    = .keep_present(genes_FB_NormalPlus),
  KC_Homeo         = .keep_present(genes_KC_Homeo),
  Endo_Quiescent   = .keep_present(genes_Endo_Quiescent),
  TGFB_targets     = .keep_present(genes_TGFB_targets),
  Tx_Up            = .keep_present(genes_Tx_Up)
)

# 友善提醒：列出每個 signature 實際可用基因數
sig_sizes <- vapply(signatures, length, integer(1))
message("Signature sizes (present genes only): ", paste(names(sig_sizes), sig_sizes, sep="=", collapse="; "))

# 若未 Normalize/Scale，至少先 Normalize（避免 AddModuleScore 走 raw）
if (!"data" %in% slotNames(GetAssay(seuratObject_Sample))) {
  seuratObject_Sample <- NormalizeData(seuratObject_Sample, verbose = FALSE)
}


# ---------- 2) AddModuleScore（相容版：不傳 k/seed/search） ----------
add_one_score <- function(obj, genes, name){
  genes <- intersect(genes, rownames(obj))
  if (length(genes) < 2) {
    warning(sprintf("Signature '%s' 只有 %d 個可用基因，略過。", name, length(genes)))
    return(obj)
  }
  Seurat::AddModuleScore(
    object   = obj,
    features = list(genes),
    name     = name,
    assay    = DefaultAssay(obj),
    nbin     = 24,
    ctrl     = 50
  )
}

score_names <- character(0)
for (nm in names(signatures)) {
  seuratObject_Sample <- add_one_score(seuratObject_Sample, signatures[[nm]], paste0("MS_", nm))
  added_col <- paste0("MS_", nm, "1") # Seurat 會自動加 "1"
  if (added_col %in% colnames(seuratObject_Sample[[]])) {
    score_names <- c(score_names, added_col)
  }
}
if (length(score_names) == 0) stop("沒有成功產生任何 ModuleScore 欄位，請檢查 signatures 與基因名稱大小寫。")

# ---------- 3) 匯總到 sample 層級（median） + Z-score ----------
score_df <- FetchData(seuratObject_Sample, vars = c("sample_id", score_names))
score_sample <- score_df %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(across(all_of(score_names), median, na.rm = TRUE), .groups = "drop")

# Z-score（跨 signature 可比）
score_z <- scale(score_sample[,-1]) %>% as.data.frame()
colnames(score_z) <- gsub("^MS_", "", gsub("1$", "", colnames(score_sample)[-1]))
score_sample_z <- cbind(sample_id = score_sample$sample_id, score_z)

# ---------- 4) 建立加權分類分數 ----------
# ---------- 4) 建立加權分類分數（強化 Normal 權重） ----------
# 權重可依資料微調；此組合通常能拉高 Normal 的辨識力，同時抑制假陽性
w_kf    <- 1.00   # 纖維化
w_il6   <- 0.70   # 發炎
w_tgf   <- 1.00   # TGFβ 目標
w_tx    <- 1.00   # 治療後上升
w_fbN   <- 0.60   # FB Normal（SPF/SRF/Plus）
w_kcN   <- 0.40   # KC Homeostasis
w_endoN <- 0.30   # Endothelial quiescence

score_sample_z <- score_sample_z %>%
  mutate(
    # Keloid：纖維化 + 發炎 + TGFβ，上升的 Tx_Up 略作扣分（治療相反）
    Score_Keloid =  w_kf * KF_MF + w_il6 * IL6_STAT3 + w_tgf * TGFB_targets - 0.3 * w_tx * Tx_Up,
    
    # Normal：Dermal FB（SPF、SRF、FB_NormalPlus）＋ KC_Homeo ＋ Endo_Quiescent
    #         再扣掉纖維化/發炎（避免「發炎很低就被當 Normal」）
    Score_Normal =  w_fbN * ( (SPF + SRF + FB_NormalPlus) / 3 ) +
      w_kcN * KC_Homeo +
      w_endoN * Endo_Quiescent -
      0.6 * (KF_MF + IL6_STAT3),
    
    # Post-Tx：Tx_Up 抬升；並對纖維化/發炎/TGFβ 予以負向權重
    Score_PostTx =  1.2 * w_tx * Tx_Up - 0.6 * (KF_MF + IL6_STAT3 + TGFB_targets)
  ) %>%
  mutate(
    Pred_Label = c("Keloid","Normal","Post-Tx")[
      max.col(cbind(Score_Keloid, Score_Normal, Score_PostTx), ties.method = "first")
    ]
  )

# ---------- 4b) 平衡指派（可開關 enforce_expected） ----------
score_sample_z$Balanced_Label <- NA_character_

if (isTRUE(enforce_expected)) {
  long_scores <- score_sample_z %>%
    dplyr::select(sample_id, Score_Keloid, Score_Normal, Score_PostTx) %>%
    tidyr::pivot_longer(-sample_id, names_to = "Class", values_to = "Score") %>%
    mutate(Class = dplyr::recode(Class,
                                 "Score_Keloid" = "Keloid", "Score_Normal" = "Normal", "Score_PostTx" = "Post-Tx"
    ))
  
  caps <- expected_counts
  assignments <- setNames(rep(NA_character_, nrow(score_sample_z)), score_sample_z$sample_id)
  ord <- order(-long_scores$Score)
  for (idx in ord) {
    sid <- long_scores$sample_id[idx]
    cls <- long_scores$Class[idx]
    if (is.na(assignments[sid]) && !is.na(caps[cls]) && caps[cls] > 0L) {
      assignments[sid] <- cls
      caps[cls] <- caps[cls] - 1L
    }
    if (all(!is.na(caps) & caps == 0L)) break
  }
  fallback <- is.na(assignments)
  if (any(fallback)) {
    assignments[fallback] <- score_sample_z$Pred_Label[match(names(assignments)[fallback], score_sample_z$sample_id)]
  }
  score_sample_z$Balanced_Label <- assignments[score_sample_z$sample_id]
}

# 最終標籤（依參數決定是否用平衡）
score_sample_z$Final_Label <- if (isTRUE(enforce_expected)) score_sample_z$Balanced_Label else score_sample_z$Pred_Label
print(score_sample_z)

# ---------- 4c) 寫回 Seurat 物件（sample_type） ----------
sample_to_type <- setNames(score_sample_z$Final_Label, score_sample_z$sample_id)
seuratObject_Sample$sample_type <- unname(sample_to_type[seuratObject_Sample$sample_id])
DimPlot(seuratObject_Sample, group.by = "sample_type", label = TRUE, reduction = "umap")
DimPlot(seuratObject_Sample, group.by = "Cell_Type", label = TRUE, reduction = "umap")

# ---------- 5) 輸出表格 ----------
write.csv(score_sample_z, file.path(out_dir, "SampleLevel_ModuleScore_Prediction.csv"), row.names = FALSE)
message("已輸出：", file.path(out_dir, "SampleLevel_ModuleScore_Prediction.csv"))
message("Pred counts: "); print(table(score_sample_z$Pred_Label))
if (isTRUE(enforce_expected)) { message("Balanced counts: "); print(table(score_sample_z$Balanced_Label)) }
message("Final counts: "); print(table(score_sample_z$Final_Label))
message("Obj sample_type counts: "); print(table(seuratObject_Sample$sample_type, useNA = "ifany"))

# ---------- 6) 畫圖與存檔：UMAP + ModuleScore 投影 ----------
# 若沒有 UMAP，就快速計算一個
if (!"umap" %in% names(seuratObject_Sample@reductions)) {
  message("No UMAP found — computing PCA/UMAP for visualization...")
  if (!"pca" %in% names(seuratObject_Sample@reductions)) {
    seuratObject_Sample <- NormalizeData(seuratObject_Sample, verbose = FALSE)
    seuratObject_Sample <- FindVariableFeatures(seuratObject_Sample, verbose = FALSE)
    seuratObject_Sample <- ScaleData(seuratObject_Sample, verbose = FALSE)
    seuratObject_Sample <- RunPCA(seuratObject_Sample, npcs = 30, verbose = FALSE)
  }
  seuratObject_Sample <- RunUMAP(seuratObject_Sample, dims = 1:min(30, ncol(Embeddings(seuratObject_Sample[["pca"]]))), verbose = FALSE)
}

# 6a) Final_Label 的 UMAP
p_label <- DimPlot(seuratObject_Sample, group.by = "sample_type", label = TRUE) +
  ggtitle("Final Sample Type (per cell)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path(out_dir, "UMAP_Final_SampleType.png"), plot = p_label, width = 7, height = 6, dpi = 300)
ggsave(filename = file.path(out_dir, "UMAP_Final_SampleType.pdf"), plot = p_label, width = 7, height = 6)

# 6b) 各 ModuleScore 的 UMAP（逐一輸出）
# 把剛剛的 score_names（例如 MS_KF_MF1, MS_IL6_STAT31, ...）拿來投影
for (sc in score_names) {
  if (!sc %in% colnames(seuratObject_Sample@meta.data)) next
  p_sc <- FeaturePlot(seuratObject_Sample, features = sc, reduction = "umap") +
    ggtitle(paste0("UMAP - ", sc))
  fn_base <- file.path(out_dir, paste0("UMAP_", sc))
  ggsave(filename = paste0(fn_base, ".png"), plot = p_sc, width = 7, height = 6, dpi = 300)
  ggsave(filename = paste0(fn_base, ".pdf"), plot = p_sc, width = 7, height = 6)
}

# 6c) （可選）把主要幾個分數合併成一頁快速檢視
top_show <- intersect(c("MS_KF_MF1","MS_IL6_STAT31","MS_TGFB_targets1","MS_Tx_Up1","MS_SPF1","MS_SRF1"), score_names)
if (length(top_show) > 0) {
  plist <- lapply(top_show, function(sc) FeaturePlot(seuratObject_Sample, features = sc, reduction = "umap") + ggtitle(sc))
  p_combined <- wrap_plots(plist, ncol = min(3, length(plist)))
  ggsave(filename = file.path(out_dir, "UMAP_ModuleScores_combined.png"), plot = p_combined, width = 6*min(3, length(plist)), height = 5*ceiling(length(plist)/min(3, length(plist))), dpi = 300)
  ggsave(filename = file.path(out_dir, "UMAP_ModuleScores_combined.pdf"), plot = p_combined, width = 6*min(3, length(plist)), height = 5*ceiling(length(plist)/min(3, length(plist))))
}

message("所有圖已輸出到資料夾：", normalizePath(out_dir))
