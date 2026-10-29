## ============================================================================
## Tumor microenvironment: tumor-infiltrating immune cells, TAMs and CAFs
## 12 mouse tongue samples (6 tumor, 6 normal), integrated object from
## integrate_12_samples.R
##
##   1. Helpers (subset + re-integrate + cluster + score marker sets)
##   2. Immune compartment: T/NK/B, monocytes, macrophages (TAM states), DCs,
##      neutrophils, mast  -> subtype per cluster
##   3. Fibroblast compartment: normal fibroblasts vs CAF states
##      (myCAF / iCAF / apCAF)  -> subtype per cluster
##   4. Proportions tumor vs normal (speckle::propeller, samples = replicates)
##   5. TAM and CAF signature scores per sample
##   6. Pseudobulk DE tumor vs normal within each subtype (DESeq2)
## Review the DotPlots after step 2/3 and fix labels in `immune_fix` / `fib_fix`.
## ============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
## once: BiocManager::install(c("speckle", "limma", "DESeq2"))
library(speckle)
library(DESeq2)
library(dplyr, warn.conflicts = FALSE)       # load last so dplyr::count/filter win
count <- dplyr::count

int_dir <- "D:/cross_species_paper/integrated_12"
tme_dir <- file.path(int_dir, "TME")
dir.create(tme_dir, showWarnings = FALSE)
obj <- readRDS(file.path(int_dir, "tongue12_integrated_annotated.rds"))
cond_cols <- c(Normal = "#2E86C1", Tumor = "#C0392B")
set.seed(1)

## =============================================================================
## 1. Helpers
## =============================================================================
reprocess <- function(x, res = 0.6, dims = 1:20) {
  x[["RNA"]] <- split(x[["RNA"]], f = x$sample)
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, nfeatures = 2000, verbose = FALSE)
  x <- ScaleData(x, verbose = FALSE)
  x <- RunPCA(x, npcs = 30, verbose = FALSE)
  x <- IntegrateLayers(x, method = RPCAIntegration, orig.reduction = "pca",
                       new.reduction = "rpca_sub", dims = dims, k.weight = 50, verbose = FALSE)
  x[["RNA"]] <- JoinLayers(x[["RNA"]])
  x <- FindNeighbors(x, reduction = "rpca_sub", dims = dims, verbose = FALSE)
  x <- FindClusters(x, resolution = res, verbose = FALSE)
  x <- RunUMAP(x, reduction = "rpca_sub", dims = dims, verbose = FALSE)
  x
}

label_clusters <- function(x, sets, fix = NULL) {
  sets <- lapply(sets, intersect, rownames(x)); sets <- sets[lengths(sets) > 0]
  x <- AddModuleScore(x, features = sets, name = "SET_", seed = 1)
  cols <- paste0("SET_", seq_along(sets))
  cl <- x@meta.data %>% dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(n = dplyr::n(), pct_tumor = round(100 * mean(condition == "Tumor"), 1),
                     dplyr::across(dplyr::all_of(cols), median), .groups = "drop")
  sc <- as.matrix(cl[, cols]); ord <- t(apply(sc, 1, order, decreasing = TRUE))
  cl$subtype <- names(sets)[ord[, 1]]
  cl$second  <- names(sets)[ord[, 2]]
  cl$gap     <- round(sc[cbind(seq_len(nrow(sc)), ord[, 1])] - sc[cbind(seq_len(nrow(sc)), ord[, 2])], 2)
  if (!is.null(fix)) for (k in names(fix)) cl$subtype[cl$seurat_clusters == k] <- fix[[k]]
  x$subtype <- cl$subtype[match(x$seurat_clusters, cl$seurat_clusters)]
  x@meta.data[, cols] <- NULL
  list(obj = x, table = cl %>% dplyr::select(seurat_clusters, n, pct_tumor, subtype, second, gap))
}

save_plots <- function(x, sets, prefix) {
  p1 <- DimPlot(x, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle(prefix)
  p2 <- DimPlot(x, group.by = "condition", cols = cond_cols, shuffle = TRUE)
  p3 <- DimPlot(x, group.by = "seurat_clusters", label = TRUE) + NoLegend()
  ggsave(file.path(tme_dir, paste0(prefix, "_UMAP.png")), p1 | p2 | p3, width = 20, height = 6, dpi = 180)
  feats <- lapply(sets, intersect, rownames(x)); feats <- feats[lengths(feats) > 0]
  feats <- lapply(seq_along(feats), function(i) setdiff(feats[[i]], unlist(feats[seq_len(i - 1)])))
  names(feats) <- names(lapply(sets, intersect, rownames(x))[lengths(lapply(sets, intersect, rownames(x))) > 0])
  feats <- feats[lengths(feats) > 0]
  p4 <- DotPlot(x, features = feats, group.by = "seurat_clusters") + RotatedAxis() +
    theme(strip.text.x = element_text(angle = 90, size = 7))
  ggsave(file.path(tme_dir, paste0(prefix, "_DotPlot_clusters.png")), p4,
         width = max(10, 0.28 * length(unlist(feats)) + 4), height = 7, dpi = 180)
}

## =============================================================================
## 2. Immune compartment
## =============================================================================
immune_sets <- list(
  CD4_T          = c("Cd3e", "Cd4", "Cd40lg", "Il7r"),
  CD8_T          = c("Cd3e", "Cd8a", "Cd8b1", "Gzmk"),
  Treg           = c("Foxp3", "Il2ra", "Ctla4", "Ikzf2"),
  Exhausted_T    = c("Pdcd1", "Havcr2", "Lag3", "Tox", "Tigit"),
  gd_T           = c("Trdc", "Tcrg-C1", "Il17a", "Rorc"),
  NK             = c("Nkg7", "Klrb1c", "Ncr1", "Gzma", "Prf1"),
  ILC2           = c("Gata3", "Il1rl1", "Areg", "Klrg1"),
  B              = c("Cd79a", "Ms4a1", "Cd19", "Pax5"),
  Plasma         = c("Jchain", "Mzb1", "Xbp1", "Sdc1"),
  Monocytes      = c("Ly6c2", "Ccr2", "Plac8", "Chil3"),
  Resident_Mac   = c("Lyve1", "Folr2", "F13a1", "Cd163", "Mrc1"),
  TAM_Spp1_Arg1  = c("Spp1", "Arg1", "Trem2", "Apoe", "Fn1", "Vegfa"),
  TAM_C1q_MHCII  = c("C1qa", "C1qb", "C1qc", "H2-Ab1", "Cd74"),
  Inflam_Mac     = c("Il1b", "Cxcl2", "Tnf", "Nfkbia", "Ccl4"),
  IFN_Mac        = c("Isg15", "Ifit3", "Cxcl10", "Rsad2"),
  cDC1           = c("Xcr1", "Clec9a", "Cd24a", "Irf8"),
  cDC2           = c("Cd209a", "Clec10a", "Mgl2", "Irf4"),
  mregDC         = c("Ccr7", "Fscn1", "Ccl22", "Cd200"),
  pDC            = c("Siglech", "Bst2", "Ccr9", "Klk1"),
  Neutrophils    = c("S100a8", "S100a9", "Retnlg", "Csf3r", "Cxcr2"),
  Mast           = c("Cpa3", "Mcpt4", "Tpsb2", "Ms4a2"),
  Basophils      = c("Mcpt8", "Prss34", "Cd200r3"),
  Proliferating  = c("Mki67", "Top2a", "Stmn1"))

imm <- subset(obj, subset = cell_type == "Immune")
message("Immune cells: ", ncol(imm))
imm <- reprocess(imm, res = 0.8)
immune_fix <- list()                          # after review: list("7" = "Exhausted_T", ...)
r <- label_clusters(imm, immune_sets, immune_fix); imm <- r$obj
print(r$table, n = 80)
write.csv(r$table, file.path(tme_dir, "immune_cluster_labels.csv"), row.names = FALSE)
save_plots(imm, immune_sets, "Immune")
saveRDS(imm, file.path(tme_dir, "immune_subset.rds"))

## =============================================================================
## 3. Fibroblast compartment (normal fibroblasts vs CAF states)
## =============================================================================
fib_sets <- list(
  Universal_Pi16   = c("Pi16", "Dpt", "Ly6c1", "Cd34", "Clec3b"),
  Matrix_Col15a1   = c("Col15a1", "Col4a1", "Lama2", "Col3a1"),
  myCAF            = c("Acta2", "Tagln", "Postn", "Col11a1", "Lrrc15", "Tnc", "Fap"),
  iCAF             = c("Il6", "Cxcl1", "Cxcl12", "Has1", "Ccl2", "Lif"),
  apCAF            = c("H2-Ab1", "H2-Aa", "Cd74", "Saa3", "Slpi"),
  ECM_remodel_CAF  = c("Mmp3", "Mmp13", "Timp1", "Sfrp2", "Inhba"),
  Proliferating    = c("Mki67", "Top2a", "Stmn1"))

fibs <- subset(obj, subset = cell_type == "Fibroblast")
message("Fibroblasts: ", ncol(fibs))
fibs <- reprocess(fibs, res = 0.6)
fib_fix <- list()
r <- label_clusters(fibs, fib_sets, fib_fix); fibs <- r$obj
print(r$table, n = 60)
write.csv(r$table, file.path(tme_dir, "fibroblast_cluster_labels.csv"), row.names = FALSE)
save_plots(fibs, fib_sets, "Fibroblast")
saveRDS(fibs, file.path(tme_dir, "fibroblast_subset.rds"))
## Tip: a cluster with pct_tumor near 100 is tumor-specific (CAF / TAM candidate);
##      near 50 = shared state; near 0 = lost in tumors.

## =============================================================================
## 4. Proportions tumor vs normal (propeller; samples are the replicates)
## =============================================================================
run_propeller <- function(x, prefix) {
  res <- propeller(clusters = x$subtype, sample = x$sample, group = x$condition, transform = "logit")
  res$cell_type <- rownames(res)
  write.csv(res, file.path(tme_dir, paste0(prefix, "_propeller.csv")), row.names = FALSE)
  print(res[, c("cell_type", "PropMean.Normal", "PropMean.Tumor", "PropRatio", "FDR")])
  pr <- x@meta.data %>% dplyr::count(sample, condition, subtype) %>%
    dplyr::group_by(sample) %>% dplyr::mutate(pct = 100 * n / sum(n)) %>% dplyr::ungroup()
  p <- ggplot(pr, aes(condition, pct, colour = condition)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, size = 2) +
    facet_wrap(~ subtype, scales = "free_y") + scale_colour_manual(values = cond_cols) +
    labs(x = NULL, y = paste("% of", prefix, "cells")) + theme_classic() + NoLegend()
  ggsave(file.path(tme_dir, paste0(prefix, "_proportions.png")), p, width = 13, height = 9, dpi = 180)
  res
}
prop_imm <- run_propeller(imm, "Immune")
prop_fib <- run_propeller(fibs, "Fibroblast")
## immune cells as % of ALL cells per sample (infiltration)
infil <- obj@meta.data %>% dplyr::group_by(sample, condition) %>%
  dplyr::summarise(pct_immune = 100 * mean(cell_type == "Immune"),
                   pct_fibroblast = 100 * mean(cell_type == "Fibroblast"), .groups = "drop")
print(infil)
write.csv(infil, file.path(tme_dir, "immune_fibroblast_fraction_per_sample.csv"), row.names = FALSE)

## =============================================================================
## 5. TAM and CAF signature scores per sample
## =============================================================================
sig <- list(
  TAM_M2like   = c("Mrc1", "Cd163", "Arg1", "Retnla", "Chil3", "Folr2"),
  TAM_M1like   = c("Nos2", "Il1b", "Tnf", "Cxcl9", "Cxcl10", "Cd86"),
  TAM_SPP1     = c("Spp1", "Trem2", "Apoe", "Gpnmb", "Lgals3"),
  T_exhaustion = c("Pdcd1", "Havcr2", "Lag3", "Tox", "Tigit", "Entpd1"),
  T_cytotoxic  = c("Gzmb", "Gzma", "Prf1", "Ifng", "Nkg7"),
  CAF_myo      = c("Acta2", "Tagln", "Postn", "Col11a1", "Lrrc15"),
  CAF_inflam   = c("Il6", "Cxcl1", "Cxcl12", "Has1", "Ccl2"))
score_compartment <- function(x, sets, keep_subtypes, prefix) {
  x <- subset(x, subset = subtype %in% keep_subtypes)
  sets <- lapply(sets, intersect, rownames(x)); sets <- sets[lengths(sets) >= 3]
  x <- AddModuleScore(x, features = sets, name = "SIG_", seed = 1)
  cols <- paste0("SIG_", seq_along(sets))
  per_sample <- x@meta.data %>% dplyr::group_by(sample, condition) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(cols), mean), .groups = "drop")
  colnames(per_sample)[match(cols, colnames(per_sample))] <- names(sets)
  write.csv(per_sample, file.path(tme_dir, paste0(prefix, "_signature_scores_per_sample.csv")), row.names = FALSE)
  long <- tidyr::pivot_longer(per_sample, names(sets), names_to = "signature", values_to = "score")
  tests <- long %>% dplyr::group_by(signature) %>%
    dplyr::summarise(p_wilcox = wilcox.test(score[condition == "Tumor"], score[condition == "Normal"])$p.value)
  print(tests)
  p <- ggplot(long, aes(condition, score, colour = condition)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2) + facet_wrap(~ signature, scales = "free_y") +
    scale_colour_manual(values = cond_cols) + theme_classic() + NoLegend() +
    labs(x = NULL, y = "Mean module score per sample")
  ggsave(file.path(tme_dir, paste0(prefix, "_signature_scores.png")), p, width = 10, height = 6, dpi = 180)
}
mac_types <- c("Resident_Mac", "TAM_Spp1_Arg1", "TAM_C1q_MHCII", "Inflam_Mac", "IFN_Mac", "Monocytes")
t_types   <- c("CD4_T", "CD8_T", "Treg", "Exhausted_T")
score_compartment(imm, sig[c("TAM_M2like", "TAM_M1like", "TAM_SPP1")], mac_types, "Macrophage")
score_compartment(imm, sig[c("T_exhaustion", "T_cytotoxic")], t_types, "Tcell")
score_compartment(fibs, sig[c("CAF_myo", "CAF_inflam")], unique(fibs$subtype), "Fibroblast")

## =============================================================================
## 6. Pseudobulk DE tumor vs normal within each subtype (needs >= 3 samples per group)
## =============================================================================
pb_de <- function(x, prefix, min_cells = 20) {
  out <- list()
  for (st in sort(unique(x$subtype))) {
    xs <- subset(x, subset = subtype == st)
    n_per <- table(xs$sample)
    ok <- names(n_per)[n_per >= min_cells]
    xs <- subset(xs, subset = sample %in% ok)
    grp <- tapply(xs$condition, xs$sample, unique)
    if (sum(grp == "Tumor") < 3 || sum(grp == "Normal") < 3) { message(st, ": not enough samples"); next }
    pb <- AggregateExpression(xs, assays = "RNA", group.by = "sample", return.seurat = FALSE)$RNA
    colnames(pb) <- sub("^g", "", colnames(pb))         # Seurat may prefix names starting with a letter
    cd <- data.frame(row.names = colnames(pb),
                     condition = factor(grp[colnames(pb)], levels = c("Normal", "Tumor")),
                     sex = factor(ifelse(substr(colnames(pb), 1, 1) == "F", "F", "M")))
    pb <- pb[rowSums(pb >= 10) >= 3, ]
    dds <- DESeqDataSetFromMatrix(round(as.matrix(pb)), cd, design = ~ sex + condition)
    dds <- DESeq(dds, quiet = TRUE)
    rr <- as.data.frame(results(dds, contrast = c("condition", "Tumor", "Normal")))
    rr$gene <- rownames(rr); rr$subtype <- st
    rr <- rr[order(rr$padj), ]
    out[[st]] <- rr
    message(st, ": ", sum(rr$padj < 0.05, na.rm = TRUE), " DE genes (padj < 0.05)")
  }
  res <- dplyr::bind_rows(out)
  write.csv(res, file.path(tme_dir, paste0(prefix, "_pseudobulk_DE_tumor_vs_normal.csv")), row.names = FALSE)
  res
}
de_imm <- pb_de(imm, "Immune")
de_fib <- pb_de(fibs, "Fibroblast")

message("Done. Outputs in ", tme_dir)