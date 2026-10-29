## ============================================================================
## Compare RPCA, CCA and Harmony integration on the 12 mouse tongue samples
##
## Labels used to judge biology come from the PER-SAMPLE annotation
## (<sample>_seurat_annotated.rds: main_type), so they do not depend on any
## integration. Malignant epithelium is excluded from batch-mixing metrics:
## it exists only in tumors and should NOT mix with normal cells.
##
## Metrics per method
##   Batch mixing   iLISI (sample) on non-malignant cells  -> higher = better mixed
##                  % of clusters containing cells from >= 10 of 12 samples
##   Biology kept   cLISI (cell type)                       -> close to 1 = types kept apart
##                  ARI / NMI between clusters and per-sample labels
##                  silhouette of cell types in the embedding
##                  purity of rare types (Schwann, satellite, lymphatic, mesenchymal-like...)
##                  malignant vs normal epithelium separation
## ============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
## once: install.packages(c("harmony", "aricode", "cluster"))
##       remotes::install_github("immunogenomics/lisi")
library(harmony)
library(lisi)
library(aricode)
library(cluster)

int_dir  <- "D:/cross_species_paper/integrated_12"
out_root <- "D:/cross_species_paper/copykat_out"
cmp_dir  <- file.path(int_dir, "integration_comparison")
dir.create(cmp_dir, showWarnings = FALSE)
dims <- 1:30
set.seed(1)

## =============================================================================
## 1. Object: the annotated RPCA object from integrate_12_samples.R
## =============================================================================
obj <- readRDS(file.path(int_dir, "tongue12_integrated_annotated.rds"))

## integration-independent labels from the per-sample annotation
lab <- bind_rows(lapply(unique(obj$sample), function(s) {
  f <- file.path(out_root, s, paste0(s, "_seurat_annotated.rds"))
  if (!file.exists(f)) return(NULL)
  o <- readRDS(f)
  data.frame(cell = paste0(s, "_", colnames(o)), sample_label = o$main_type)
}))
obj$sample_label <- lab$sample_label[match(colnames(obj), lab$cell)]
## split epithelium into malignant / normal (from the CopyKAT-based status)
obj$sample_label <- ifelse(obj$sample_label %in% "Epithelial" & obj$epi_status == "Malignant_epithelium",
                           "Malignant_epithelium",
                           ifelse(obj$sample_label %in% "Epithelial", "Normal_epithelium", obj$sample_label))
obj$sample_label[is.na(obj$sample_label)] <- "Unlabeled"
print(table(obj$sample_label))

## =============================================================================
## 2. Run the three integrations on the same PCA
## =============================================================================
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)      # IntegrateLayers needs split layers
if (!"pca" %in% Reductions(obj)) {
  obj <- NormalizeData(obj); obj <- FindVariableFeatures(obj, nfeatures = 3000)
  obj <- ScaleData(obj); obj <- RunPCA(obj, npcs = 50)
}

timing <- list()
if (!"integrated.rpca" %in% Reductions(obj)) {
  t0 <- Sys.time()
  obj <- IntegrateLayers(obj, method = RPCAIntegration, orig.reduction = "pca",
                         new.reduction = "integrated.rpca", dims = dims, verbose = FALSE)
  timing$RPCA <- difftime(Sys.time(), t0, units = "mins")
}
t0 <- Sys.time()
obj <- IntegrateLayers(obj, method = CCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.cca", dims = dims, verbose = FALSE)
timing$CCA <- difftime(Sys.time(), t0, units = "mins")
t0 <- Sys.time()
obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca",
                       new.reduction = "harmony", verbose = FALSE)
timing$Harmony <- difftime(Sys.time(), t0, units = "mins")
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

methods <- c(Unintegrated = "pca", RPCA = "integrated.rpca", CCA = "integrated.cca", Harmony = "harmony")
for (m in names(methods)) {
  red <- methods[[m]]
  obj <- FindNeighbors(obj, reduction = red, dims = dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, cluster.name = paste0("clu_", m), verbose = FALSE)
  obj <- RunUMAP(obj, reduction = red, dims = dims, reduction.name = paste0("umap_", m), verbose = FALSE)
}
saveRDS(obj, file.path(cmp_dir, "tongue12_three_integrations.rds"))

## =============================================================================
## 3. Metrics
## =============================================================================
md <- obj@meta.data
keep_lab <- md$sample_label != "Unlabeled"
nonmal   <- keep_lab & md$sample_label != "Malignant_epithelium"
rare <- c("Schwann", "Satellite", "Lymphatic_EC", "Mesenchymal_like", "Skeletal_muscle", "Mural")
sub_sil <- sample(which(keep_lab), min(15000, sum(keep_lab)))

res <- list(); purity_tab <- list()
for (m in names(methods)) {
  emb <- Embeddings(obj, methods[[m]])[, dims]
  clu <- md[[paste0("clu_", m)]]

  ## batch mixing (non-malignant cells only)
  il <- compute_lisi(emb[nonmal, ], md[nonmal, , drop = FALSE], "sample")$sample
  cl_samples <- tapply(md$sample, clu, function(x) n_distinct(x))
  ## biology
  cl_ <- compute_lisi(emb[keep_lab, ], md[keep_lab, , drop = FALSE], "sample_label")$sample_label
  ari <- ARI(clu[keep_lab], md$sample_label[keep_lab])
  nmi <- NMI(clu[keep_lab], md$sample_label[keep_lab])
  sil <- mean(silhouette(as.integer(factor(md$sample_label[sub_sil])), dist(emb[sub_sil, ]))[, 3])

  ## malignant vs normal epithelium: % of malignant cells whose cluster is mostly malignant
  epi <- md$sample_label %in% c("Malignant_epithelium", "Normal_epithelium")
  maj <- tapply(md$sample_label[epi], clu[epi], function(x) names(which.max(table(x))))
  mal_sep <- mean(maj[as.character(clu[epi & md$sample_label == "Malignant_epithelium"])] == "Malignant_epithelium")

  ## rare types: % of their cells sitting in a cluster where they are the majority
  for (r in intersect(rare, unique(md$sample_label))) {
    idx <- md$sample_label == r
    major <- tapply(md$sample_label[keep_lab], clu[keep_lab], function(x) names(which.max(table(x))))
    purity_tab[[paste(m, r)]] <- data.frame(method = m, cell_type = r, n = sum(idx),
      pct_in_own_cluster = round(100 * mean(major[as.character(clu[idx])] == r), 1))
  }

  res[[m]] <- data.frame(
    method = m,
    n_clusters = length(unique(clu)),
    iLISI_sample_median = round(median(il), 2),                       # max = 12
    pct_clusters_ge10_samples = round(100 * mean(cl_samples >= 10), 1),
    cLISI_celltype_median = round(median(cl_), 3),                     # 1 = perfectly separated
    ARI = round(ari, 3), NMI = round(nmi, 3),
    silhouette_celltype = round(sil, 3),
    malignant_epi_separated_pct = round(100 * mal_sep, 1),
    minutes = if (!is.null(timing[[m]])) round(as.numeric(timing[[m]]), 1) else NA)
}
res <- bind_rows(res); purity_tab <- bind_rows(purity_tab)
print(res); print(tidyr::pivot_wider(purity_tab, id_cols = cell_type, names_from = method,
                                      values_from = pct_in_own_cluster))
write.csv(res, file.path(cmp_dir, "integration_metrics.csv"), row.names = FALSE)
write.csv(purity_tab, file.path(cmp_dir, "rare_type_purity.csv"), row.names = FALSE)

## =============================================================================
## 4. Plots
## =============================================================================
cols_cond <- c(Normal = "#2E86C1", Tumor = "#C0392B")
panel <- lapply(names(methods), function(m) {
  r <- paste0("umap_", m)
  (DimPlot(obj, reduction = r, group.by = "sample", shuffle = TRUE, pt.size = 0.1) + ggtitle(paste(m, "- sample")) + NoLegend()) |
  (DimPlot(obj, reduction = r, group.by = "condition", cols = cols_cond, shuffle = TRUE, pt.size = 0.1) + ggtitle("condition") + NoLegend()) |
  (DimPlot(obj, reduction = r, group.by = "sample_label", label = TRUE, repel = TRUE, label.size = 3, pt.size = 0.1) +
     ggtitle("per-sample cell type") + NoLegend())
})
ggsave(file.path(cmp_dir, "UMAP_comparison_grid.png"), wrap_plots(panel, ncol = 1),
       width = 18, height = 22, dpi = 150)

## metric summary plot
long <- tidyr::pivot_longer(res, c(iLISI_sample_median, pct_clusters_ge10_samples, ARI, NMI,
                                   silhouette_celltype, malignant_epi_separated_pct),
                            names_to = "metric", values_to = "value")
p_m <- ggplot(long, aes(method, value, fill = method)) + geom_col() +
  facet_wrap(~ metric, scales = "free_y", nrow = 2) + theme_classic() + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = NULL, y = NULL)
ggsave(file.path(cmp_dir, "integration_metrics.png"), p_m, width = 11, height = 6, dpi = 200)

p_r <- ggplot(purity_tab, aes(cell_type, pct_in_own_cluster, fill = method)) +
  geom_col(position = "dodge") + theme_classic() + labs(x = NULL, y = "% cells in own cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(cmp_dir, "rare_type_purity.png"), p_r, width = 9, height = 4.5, dpi = 200)
message("Done. Outputs in ", cmp_dir)
