## ============================================================================
## Mouse tongue - integration of all 12 samples (6 tumor, 6 normal)
## Adds: sample, sex, condition, cell type, epithelial status (malignant / normal)
## and the CopyKAT call where a cell was tested.
##
##   1. Read + QC per sample (+ optional scDblFinder)
##   2. Merge, RPCA integration, clustering, UMAP
##   3. Main cell types on integrated clusters (same marker sets as before)
##   4. Epithelial status from the per-sample CopyKAT work
##   5. Plots + save
## ============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

cr_root  <- "D:/cross_species_paper/mouse_cellranger"
out_root <- "D:/cross_species_paper/copykat_out"
int_dir  <- "D:/cross_species_paper/integrated_12"
dir.create(int_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c(FT  = "FT_Count_outs",  F1 = "F1_count_outs",  F2 = "F2_count_outs",
             MTT = "MTT_Count_outs", M1 = "M1_count_outs",  M2 = "M2_count_outs",
             FN1 = "FN1_count_outs", FN2 = "FN2_count_outs", FC = "FC_Count_outs",
             MN1 = "MN1_count_outs", MN2 = "MN2_count_outs", MC = "MC_Count_outs")
tumor_samples <- c("FT", "F1", "F2", "MTT", "M1", "M2")

## Epithelial clusters with a normal (Krt5+, diploid) phenotype inside tumor samples
normal_looking <- list(FT = "13", F1 = "17", F2 = "14")

run_doublets <- TRUE    # needs: BiocManager::install("scDblFinder")

## =============================================================================
## 1. Read + QC
## =============================================================================
objs <- list()
qc_log <- list()
for (s in names(samples)) {
  counts <- Read10X(file.path(cr_root, samples[[s]], "filtered_feature_bc_matrix"))
  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  o <- CreateSeuratObject(counts, project = s, min.cells = 3, min.features = 200)
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^mt-")
  n0 <- ncol(o)
  
  ## upper cut per sample (median + 3 MAD) instead of a fixed 7,000
  max_feat <- median(o$nFeature_RNA) + 3 * mad(o$nFeature_RNA)
  o <- subset(o, subset = nFeature_RNA > 200 & nFeature_RNA < max_feat & percent.mt < 20)
  
  n_dbl <- NA
  if (run_doublets) {
    library(scDblFinder); library(SingleCellExperiment)
    sce <- scDblFinder(SingleCellExperiment(list(counts = GetAssayData(o, layer = "counts"))))
    o$dbl_class <- sce$scDblFinder.class
    o$dbl_score <- sce$scDblFinder.score
    n_dbl <- sum(o$dbl_class == "doublet")
    o <- subset(o, subset = dbl_class == "singlet")
  }
  
  o$sample    <- s
  o$sex       <- ifelse(substr(s, 1, 1) == "F", "Female", "Male")
  o$condition <- ifelse(s %in% tumor_samples, "Tumor", "Normal")
  o <- RenameCells(o, add.cell.id = s)                 # barcodes become "<sample>_<barcode>"
  objs[[s]] <- o
  qc_log[[s]] <- data.frame(sample = s, cells_raw = n0, max_nFeature = round(max_feat),
                            doublets_removed = n_dbl, cells_kept = ncol(o))
  message(s, ": ", n0, " -> ", ncol(o), " cells")
}
qc_log <- bind_rows(qc_log); print(qc_log)
write.csv(qc_log, file.path(int_dir, "QC_per_sample.csv"), row.names = FALSE)

## =============================================================================
## 2. Merge + RPCA integration (Seurat v5)
## =============================================================================
obj <- merge(objs[[1]], objs[-1])
rm(objs); gc()
## after merge the RNA assay has one counts layer per sample (split), as IntegrateLayers needs

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 50)
ElbowPlot(obj, ndims = 50)                          # check that 30 PCs is reasonable

obj <- IntegrateLayers(obj, method = RPCAIntegration, orig.reduction = "pca",
                       new.reduction = "integrated.rpca", dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:30)
obj <- FindClusters(obj, resolution = c(0.2, 0.5), cluster.name = c("clusters_0.2", "clusters_0.5"))
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap")
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])           # for markers / DotPlots
Idents(obj) <- "clusters_0.5"
saveRDS(obj, file.path(int_dir, "tongue12_integrated_raw.rds"))

## =============================================================================
## 3. Main cell types on integrated clusters
## =============================================================================
main_markers <- list(
  Epithelial       = c("Krt14", "Krt5", "Krt13", "Krt4", "Krt17", "Epcam", "Krt15"),
  Endothelial      = c("Pecam1", "Cdh5", "Cldn5", "Emcn", "Kdr"),
  Lymphatic_EC     = c("Prox1", "Flt4", "Ccl21a", "Lyve1", "Mmrn1"),
  Immune           = c("Ptprc", "Cd52", "Coro1a", "Tyrobp", "Lcp1"),
  Fibroblast       = c("Pdgfra", "Col1a1", "Col1a2", "Dcn", "Lum"),
  Mural            = c("Rgs5", "Pdgfrb", "Abcc9", "Kcnj8", "Acta2", "Myh11", "Tagln"),
  Skeletal_muscle  = c("Acta1", "Tnnc2", "Tnnt3", "Tnni2", "Myl1", "Ckm"),
  Satellite        = c("Pax7", "Myf5", "Calcr", "Cdh15"),
  Schwann          = c("Sox10", "Plp1", "Mpz", "Cdh19", "S100b"),
  Mesenchymal_like = c("Igfbp6", "Tenm2", "Sfrp5", "Bnc2", "Ptch1", "Crlf1", "Itgb4"))
mk <- lapply(main_markers, intersect, rownames(obj))
obj <- AddModuleScore(obj, features = mk, name = "MAIN_", seed = 1)
sc_cols <- paste0("MAIN_", seq_along(mk))

cl <- obj@meta.data %>% group_by(clusters_0.5) %>%
  summarise(n = n(),
            pct_tumor = round(100 * mean(condition == "Tumor"), 1),
            n_samples = n_distinct(sample),
            across(all_of(sc_cols), median), .groups = "drop")
sc <- as.matrix(cl[, sc_cols])
ord <- t(apply(sc, 1, order, decreasing = TRUE))
cl$main_type <- names(mk)[ord[, 1]]
cl$second    <- names(mk)[ord[, 2]]
cl$gap       <- round(sc[cbind(seq_len(nrow(sc)), ord[, 1])] - sc[cbind(seq_len(nrow(sc)), ord[, 2])], 2)
colnames(cl)[match(sc_cols, colnames(cl))] <- paste0("score_", names(mk))
write.csv(cl %>% mutate(across(where(is.numeric), ~ round(.x, 2))),
          file.path(int_dir, "integrated_cluster_annotation.csv"), row.names = FALSE)
print(cl %>% select(clusters_0.5, n, n_samples, pct_tumor, main_type, second, gap), n = 80)

## manual corrections after looking at the DotPlot (cluster = type; use "Doublet"/"LowQC" to drop)
overrides <- list()       # e.g. list("27" = "Doublet")
for (k in names(overrides)) cl$main_type[cl$clusters_0.5 == k] <- overrides[[k]]
obj$cell_type <- cl$main_type[match(obj$clusters_0.5, cl$clusters_0.5)]
obj@meta.data[, sc_cols] <- NULL

pan <- lapply(main_markers, function(g) intersect(g[1:min(4, length(g))], rownames(obj)))
p_dot <- DotPlot(obj, features = pan, group.by = "clusters_0.5") + RotatedAxis() +
  theme(strip.text.x = element_text(angle = 90, size = 7)) + ggtitle("Integrated clusters - pan markers")
ggsave(file.path(int_dir, "DotPlot_pan_markers_integrated_clusters.png"), p_dot, width = 18, height = 9, dpi = 200)

## =============================================================================
## 4. Epithelial status (malignant / normal) from the per-sample work
## =============================================================================
## 4a. per-sample epithelial clusters chosen before (per-sample objects, same barcodes)
status_tab <- list()
for (s in names(samples)) {
  f <- file.path(out_root, s, paste0(s, "_seurat_lineage.rds"))
  if (!file.exists(f)) { message(s, ": no per-sample object yet"); next }
  o <- readRDS(f)
  epi <- colnames(o)[o$lineage == "Epithelial"]
  if (!length(epi)) next
  cl_epi <- as.character(o$seurat_clusters[match(epi, colnames(o))])
  st <- if (s %in% tumor_samples) {
    ifelse(cl_epi %in% normal_looking[[s]], "Normal_epithelium", "Malignant_epithelium")
  } else "Normal_epithelium"
  status_tab[[s]] <- data.frame(cell = paste0(s, "_", epi), epi_status = st)
}
status_tab <- bind_rows(status_tab)
obj$epi_status <- status_tab$epi_status[match(colnames(obj), status_tab$cell)]

## 4b. cells without a per-sample label (new QC let them in, e.g. >7,000 genes):
##     take the majority status of labeled cells from the same sample + integrated cluster
md <- obj@meta.data
md$cell <- rownames(md)
maj <- md %>% filter(!is.na(epi_status)) %>%
  count(sample, clusters_0.5, epi_status) %>%
  group_by(sample, clusters_0.5) %>% slice_max(n, n = 1, with_ties = FALSE) %>%
  select(sample, clusters_0.5, maj_status = epi_status)
md <- left_join(md, maj, by = c("sample", "clusters_0.5"))
fill <- is.na(md$epi_status) & md$cell_type == "Epithelial"
md$epi_status[fill] <- md$maj_status[fill]
md$epi_status[is.na(md$epi_status)] <- ifelse(md$cell_type[is.na(md$epi_status)] == "Epithelial",
                                              "Epithelial_unassigned", "Non_epithelial")
obj$epi_status <- md$epi_status[match(colnames(obj), md$cell)]

## 4c. the CopyKAT call itself, where a cell was tested (normal-epithelium reference runs)
pred_files <- c(
  file.path(out_root, "TUMORS_vs_NORMAL_EPI", "TumVsNormEpi_copykat_prediction.txt"),
  list.files(out_root, pattern = "_vsNormEpi_copykat_prediction.txt$", recursive = TRUE, full.names = TRUE))
pred <- bind_rows(lapply(pred_files[file.exists(pred_files)], read.delim))
pred <- pred[!duplicated(pred$cell.names), ]
key  <- c(setNames(pred$copykat.pred, pred$cell.names),
          setNames(pred$copykat.pred, sub("\\.1$", "-1", pred$cell.names)))   # "-" may have become "."
obj$copykat_call <- unname(key[colnames(obj)])
obj$copykat_call[is.na(obj$copykat_call)] <- "not_tested"

## combined label: cell type, with epithelium split into malignant / normal
obj$cell_type_status <- ifelse(obj$cell_type == "Epithelial", obj$epi_status, obj$cell_type)

## checks
print(table(obj$sample, obj$epi_status))
print(table(obj$epi_status, obj$copykat_call))   # malignant should be ~all aneuploid where tested

## =============================================================================
## 5. Plots + save
## =============================================================================
p1 <- DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE) + ggtitle("Cell type")
p2 <- DimPlot(obj, group.by = "condition", cols = c(Normal = "#2E86C1", Tumor = "#C0392B"), shuffle = TRUE)
p3 <- DimPlot(obj, group.by = "sex", shuffle = TRUE)
p4 <- DimPlot(obj, group.by = "epi_status",
              cols = c(Malignant_epithelium = "#C0392B", Normal_epithelium = "#2E86C1",
                       Epithelial_unassigned = "orange", Non_epithelial = "grey90"),
              order = c("Malignant_epithelium", "Normal_epithelium"))
ggsave(file.path(int_dir, "UMAP_overview.png"), (p1 | p2) / (p3 | p4), width = 16, height = 13, dpi = 200)

p5 <- DimPlot(obj, group.by = "sample", shuffle = TRUE) + ggtitle("Sample mixing")
ggsave(file.path(int_dir, "UMAP_by_sample.png"), p5, width = 8, height = 6, dpi = 200)

prop <- obj@meta.data %>% count(sample, condition, cell_type_status) %>%
  group_by(sample) %>% mutate(pct = 100 * n / sum(n))
prop$sample <- factor(prop$sample, levels = names(samples))
p6 <- ggplot(prop, aes(sample, pct, fill = cell_type_status)) + geom_col() +
  facet_grid(~ condition, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "% of cells", fill = NULL) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(int_dir, "Barplot_composition_by_sample.png"), p6, width = 11, height = 5, dpi = 200)
write.csv(prop, file.path(int_dir, "composition_by_sample.csv"), row.names = FALSE)

saveRDS(obj, file.path(int_dir, "tongue12_integrated_annotated.rds"))
message("Done: ", ncol(obj), " cells. Files in ", int_dir)