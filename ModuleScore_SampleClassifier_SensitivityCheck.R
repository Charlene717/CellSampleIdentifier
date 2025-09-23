## =========================================================
## ModuleScore 分組「權重敏感度檢驗」完整腳本（無 desc() 依賴版）
## 需求：已載入 Seurat 物件 seuratObject_Sample，且 meta.data 有 sample_id、Cell_Type
## =========================================================

suppressPackageStartupMessages({
  if (!require("Seurat")) install.packages("Seurat"); library(Seurat)
  if (!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
  if (!require("patchwork")) install.packages("patchwork"); library(patchwork)
})

## -------------------------------
## 0) 參數設定
## -------------------------------
stopifnot(exists("seuratObject_Sample"))
if (!"sample_id" %in% colnames(seuratObject_Sample@meta.data)) {
  stop("seuratObject_Sample@meta.data 缺少 'sample_id' 欄位")
}

# 預期每組樣本數（若不想強制平衡，設 enforce_expected <- FALSE）
enforce_expected <- TRUE
expected_counts  <- c("Keloid" = 2L, "Normal" = 2L, "Post-Tx" = 2L)

# 敏感度迭代次數與擾動強度
n_iter       <- 500        # 推薦 200~1000
sdlog_jitter <- 0.25       # 權重乘上 log-normal 噪音，數值越大擾動越強

# 輸出資料夾
out_dir <- "Sensitivity"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## -------------------------------
## 1) 基因 signature
## -------------------------------
.keep_present <- function(g) unique(intersect(g, rownames(seuratObject_Sample)))

genes_KF_MF <- c("POSTN","COL11A1","COL12A1","COL5A2","COMP","ADAM12","TNFRSF12A","SDC1",
                 "ITGBL1","FN1","THBS2","ASPN","TAGLN","ACTA2")
genes_IL6_STAT3 <- c("IL6","IL6R","STAT3","SOCS3","CCL2","CXCL1","CXCL2","CXCL3","CXCL8",
                     "NFKBIA","WNT5A","PTGS2","OSMR")
genes_SPF <- c("APCDD1","ID1","WIF1","COL13A1","COL18A1","COL23A1","PI16","CXCL14","DPP4")
genes_SRF <- c("WISP2","SLPI","TSPAN8","COL18A1","SFRP2","PDGFRB")
genes_FB_NormalPlus <- c("COL14A1","PCOLCE2","MFAP5","OGN","CFD")
genes_KC_Homeo <- c("KRT14","KRT5","DST","DSG3","ITGA6","ITGB1","COL17A1","PLEC",
                    "LAMA3","LAMB3","LAMB2","EPCAM","CLDN1","KLF4","OVOL1","GRHL1","GRHL3")
genes_Endo_Quiescent <- c("KDR","FLT1","ESAM","PECAM1","CDH5","CLDN5","RBP7","GPIHBP1","PROX1","KLF2","KLF4")
genes_TGFB_targets <- c("CTGF","SERPINE1","THBS1","COMP","TAGLN","ACTA2","COL1A1","COL1A2","COL3A1","FN1")
genes_Tx_Up <- c("MMP1","MMP3","DUSP1","KLF2","KLF6","ZFP36","TSC22D3","FKBP5")

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

message("Signature sizes: ",
        paste(names(signatures), vapply(signatures, length, integer(1)), sep="=", collapse="; "))

## -------------------------------
## 2) 取得（或計算）sample-level Z-score 表 (score_sample_z)
## -------------------------------
add_one_score <- function(obj, genes, name){
  genes <- intersect(genes, rownames(obj))
  if (length(genes) < 2) return(obj)
  Seurat::AddModuleScore(
    object   = obj,
    features = list(genes),
    name     = name,
    assay    = DefaultAssay(obj),
    nbin     = 24,
    ctrl     = 50
  )
}

compute_score_sample_z <- function(obj, signatures){
  DefaultAssay(obj) <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)
  if (!"data" %in% slotNames(GetAssay(obj))) {
    obj <- NormalizeData(obj, verbose = FALSE)
  }
  score_names <- character(0)
  for (nm in names(signatures)) {
    obj <- add_one_score(obj, signatures[[nm]], paste0("MS_", nm))
    added_col <- paste0("MS_", nm, "1")
    if (added_col %in% colnames(obj[[]])) score_names <- c(score_names, added_col)
  }
  if (length(score_names) == 0) stop("ModuleScore 皆未產生，請檢查基因名稱。")
  
  score_df <- FetchData(obj, vars = c("sample_id", score_names))
  score_sample <- score_df |>
    dplyr::group_by(sample_id) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(score_names), median, na.rm = TRUE), .groups = "drop")
  score_z <- scale(score_sample[,-1]) |> as.data.frame()
  colnames(score_z) <- gsub("^MS_", "", gsub("1$", "", colnames(score_sample)[-1]))
  score_sample_z <- cbind(sample_id = score_sample$sample_id, score_z)
  list(obj=obj, score_sample_z=score_sample_z, score_names=score_names)
}

if (!exists("score_sample_z")) {
  res0 <- compute_score_sample_z(seuratObject_Sample, signatures)
  seuratObject_Sample <- res0$obj
  score_sample_z      <- res0$score_sample_z
}

## -------------------------------
## 3) 定義打分與（可選）平衡分配
## -------------------------------
# 基準權重
base_weights <- c(
  w_kf    = 1.00,  # 纖維化
  w_il6   = 1.00,  # 發炎
  w_tgf   = 1.00,  # TGFβ
  w_tx    = 1.00,  # 治療上升群
  w_fbN   = 1.00,  # Normal: FB（SPF/SRF/Plus）
  w_kcN   = 1.00,  # Normal: Keratinocyte
  w_endoN = 1.00   # Normal: Endothelial quiescence
)

compute_labels_once <- function(score_sample_z, weights, enforce_expected=FALSE, expected_counts=NULL){
  z <- score_sample_z
  # 三類打分（沿用你的公式與係數）
  Score_Keloid <-  weights["w_kf"]  * z$KF_MF +
    weights["w_il6"] * z$IL6_STAT3 +
    weights["w_tgf"] * z$TGFB_targets -
    0.3 * weights["w_tx"] * z$Tx_Up
  
  Score_Normal <-  weights["w_fbN"]  * ((z$SPF + z$SRF + z$FB_NormalPlus)/3) +
    weights["w_kcN"]  * z$KC_Homeo +
    weights["w_endoN"]* z$Endo_Quiescent -
    0.6 * (z$KF_MF + z$IL6_STAT3)
  
  Score_PostTx <-  1.2 * weights["w_tx"] * z$Tx_Up -
    0.6 * (z$KF_MF + z$IL6_STAT3 + z$TGFB_targets)
  
  S <- cbind(Score_Keloid, Score_Normal, Score_PostTx)
  pred <- c("Keloid","Normal","Post-Tx")[max.col(S, ties.method="first")]
  
  if (!enforce_expected) return(pred)
  
  # 平衡分配：從最高分往下填滿配額（用 base R 排序避免 desc 衝突）
  long_scores <- tibble::tibble(sample_id = z$sample_id,
                                Keloid = Score_Keloid,
                                Normal = Score_Normal,
                                `Post-Tx` = Score_PostTx) |>
    tidyr::pivot_longer(-sample_id, names_to = "Class", values_to = "Score") |>
    as.data.frame()
  long_scores <- long_scores[order(-long_scores$Score), ]
  
  caps <- expected_counts
  assignments <- setNames(rep(NA_character_, nrow(z)), z$sample_id)
  for (i in seq_len(nrow(long_scores))){
    sid <- long_scores$sample_id[i]; cls <- long_scores$Class[i]
    if (is.na(assignments[sid]) && !is.na(caps[cls]) && caps[cls] > 0L){
      assignments[sid] <- cls
      caps[cls] <- caps[cls] - 1L
    }
    if (all(!is.na(caps) & caps==0L)) break
  }
  # 補上未分配者：用原始 pred
  fallback <- is.na(assignments)
  if (any(fallback)) assignments[fallback] <- pred[match(names(assignments)[fallback], z$sample_id)]
  unname(assignments[z$sample_id])
}

## 若尚未有基準 Final_Label，先以 base_weights 計一次
if (!"Final_Label" %in% colnames(score_sample_z)) {
  score_sample_z$Final_Label <- compute_labels_once(
    score_sample_z, base_weights, enforce_expected, expected_counts
  )
}

## -------------------------------
## 4) 敏感度模擬：擾動權重多次
## -------------------------------
set.seed(123)
sample_ids <- score_sample_z$sample_id
res_mat <- matrix(NA_character_, nrow=length(sample_ids), ncol=n_iter,
                  dimnames=list(sample_ids, paste0("iter", seq_len(n_iter))))

for (it in seq_len(n_iter)){
  jitter_mult <- rlnorm(length(base_weights), meanlog=0, sdlog=sdlog_jitter)
  names(jitter_mult) <- names(base_weights)
  w_perturbed <- base_weights * jitter_mult
  res_mat[,it] <- compute_labels_once(
    score_sample_z, w_perturbed, enforce_expected, expected_counts
  )
}

## -------------------------------
## 5) 穩定度統計與輸出（修正版）
## -------------------------------
# 各 sample 的類別次數 -> 直接產生 matrix（避免 list/rbind 問題）
counts_mat <- t(apply(res_mat, 1, function(v) {
  as.integer(table(factor(v, levels = c("Keloid","Normal","Post-Tx"))))
}))
colnames(counts_mat) <- c("Keloid","Normal","Post-Tx")
rownames(counts_mat) <- rownames(res_mat)

counts_df <- as.data.frame(counts_mat)
counts_df$sample_id <- rownames(counts_df)
counts_df$total <- counts_df$Keloid + counts_df$Normal + counts_df$`Post-Tx`
counts_df$prop_Keloid <- counts_df$Keloid / counts_df$total
counts_df$prop_Normal <- counts_df$Normal / counts_df$total
counts_df$prop_PostTx <- counts_df$`Post-Tx` / counts_df$total

# 穩定度（majority proportion）
majority_label <- apply(res_mat, 1, function(v) names(sort(table(v), decreasing=TRUE))[1])
majority_prop  <- apply(res_mat, 1, function(v) max(prop.table(table(v))))

summary_df <- tibble::tibble(
  sample_id       = sample_ids,
  baseline_label  = score_sample_z$Final_Label,
  majority_label  = majority_label,
  stability       = round(majority_prop, 3),
  agree_with_baseline = (baseline_label == majority_label)
)
# 用 base R 排序，避免 desc 衝突
summary_df <- summary_df[order(-summary_df$stability), ]

# 輸出 CSV
write.csv(summary_df, file.path(out_dir, "sensitivity_summary.csv"), row.names=FALSE)
write.csv(counts_df,  file.path(out_dir, "sensitivity_counts.csv"),  row.names=FALSE)

message("已輸出：", normalizePath(file.path(out_dir, "sensitivity_summary.csv")))
message("已輸出：", normalizePath(file.path(out_dir, "sensitivity_counts.csv")))

## -------------------------------
## 6) Visualization（不使用 desc 排序）
## -------------------------------
# (a) Stability barplot（依 stability 降冪的 factor level）
summary_df$sample_id <- factor(summary_df$sample_id, levels = summary_df$sample_id)
p_stab <- ggplot(summary_df, aes(x=sample_id, y=stability, fill=agree_with_baseline)) +
  geom_col(width=0.7) +
  coord_flip() +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  scale_fill_manual(values=c("TRUE"="#1F78B4","FALSE"="#E31A1C"),
                    labels=c("FALSE"="≠ Baseline","TRUE"="= Baseline")) +
  labs(x=NULL, y="Stability (majority proportion)", fill=NULL,
       title="Sample-level Stability under Perturbation") +
  theme_minimal(base_size = 13)
ggsave(file.path(out_dir, "stability_barplot.png"), p_stab, width=7, height=6, dpi=300)

# (b) Stacked barplot of class proportions（依 summary_df 的順序顯示）
counts_long <- counts_df |>
  dplyr::select(sample_id, Keloid, Normal, `Post-Tx`) |>
  tidyr::pivot_longer(-sample_id, names_to="Class", values_to="Count") |>
  dplyr::group_by(sample_id) |>
  dplyr::mutate(Prop = Count / sum(Count)) |>
  dplyr::ungroup()

counts_long$sample_id <- factor(counts_long$sample_id, levels = summary_df$sample_id)

p_stack <- ggplot(counts_long, aes(x=sample_id, y=Prop, fill=Class)) +
  geom_col(width=0.7) +
  coord_flip() +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values=c("Keloid"="#F8766D","Normal"="#00BA38","Post-Tx"="#619CFF")) +
  labs(x=NULL, y="Proportion across iterations", fill=NULL,
       title="Class Proportions under Perturbation") +
  theme_minimal(base_size = 13)
ggsave(file.path(out_dir, "stacked_props.png"), p_stack, width=7, height=6, dpi=300)

## -------------------------------
## 7) 螢幕輸出重點
## -------------------------------
print(head(summary_df, 10))
message("完成！請查看資料夾：", normalizePath(out_dir))
