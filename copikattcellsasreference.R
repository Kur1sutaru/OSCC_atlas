## ============================================================================
## Mouse tongue - per-sample annotation with the slide marker sets,
## then CopyKAT with T cells as the diploid reference
##
##   PART 1  Marker sets (main types + subtypes, from the compilation slides)
##   PART 2  Per-sample annotation: main type per cluster, subtype within type,
##           T cells gated per cell. Saves tables + dotplots for review.
##           -> REVIEW, then fix anything in `overrides` and re-run Part 2b
##   PART 3  CopyKAT per sample: test = epithelium, reference = T cells
##           (run a normal sample first as the control)
##   PART 4  Summary + comparison with the tumor signature
## ============================================================================

library(Seurat)
library(copykat)
library(dplyr)
library(ggplot2)

cr_root  <- "D:/cross_species_paper/mouse_cellranger"
out_root <- "D:/cross_species_paper/copykat_out"

samples <- c(FN1 = "FN1_count_outs", MN1 = "MN1_count_outs",       # presumed normal: run first
             FN2 = "FN2_count_outs", MN2 = "MN2_count_outs",
             FC  = "FC_Count_outs",  MC  = "MC_Count_outs",
             FT  = "FT_Count_outs",  MTT = "MTT_Count_outs",
             F1  = "F1_count_outs",  F2  = "F2_count_outs",
             M1  = "M1_count_outs",  M2  = "M2_count_outs")

## =============================================================================
## PART 1. Marker sets
## =============================================================================
## Main types (pan markers). "Mesenchymal-like" is kept separate from epithelium:
## its genes (Igfbp6, Sfrp5, Bnc2, Ptch1, Tenm2) are perineurial-fibroblast genes.
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
  Mesenchymal_like = c("Igfbp6", "Tenm2", "Sfrp5", "Bnc2", "Ptch1", "Crlf1", "Itgb4")
)

## Subtypes (genes exactly as in the slide dotplots; Neutrophils added)
sub_markers <- list(
  Epithelial = list(
    Suprabasal_keratinocytes    = c("Krt4", "Mt4", "Krt13", "Lgals7", "Dsp", "Pkp1"),
    Basal_keratinocytes         = c("Krt14", "Krt5", "Krt15", "Trp63", "Col17a1"),
    Proliferating_keratinocytes = c("Hmgb2", "Top2a", "Mki67", "Cenpf", "Stmn1", "Cks2")),
  Immune = list(
    T_cells     = c("Ptprc", "Cd3d", "Cd3g", "Cd3e"),
    Macrophages = c("Cd163", "Cd68", "Ccr2", "Cd14", "Cd74", "Tyrobp"),
    Mast        = c("Ms4a2", "Tpsb2", "Mcpt4", "Mrgprb2"),
    Neutrophils = c("S100a8", "S100a9", "Csf3r", "Retnlg")),
  Endothelial = list(
    Metab_active_capillary_EC = c("Gpihbp1", "Fabp4", "Ly6c1", "Cd36"),
    Capillary_EC              = c("Dach1", "Ptprm", "Flt1"),
    Stress_response_EC        = c("Rflnb", "Ybx1", "Unc45b", "Uba52", "Ccnd1", "Eef1a1"),
    IFN_responsive_EC         = c("Ifit1", "Isg15", "Rsad2", "Stat1", "Irf7"),
    Arterial_EC               = c("Hey1", "Sox17", "Gja4", "Efnb2"),
    Activated_arterial_EC     = c("Ier3", "Junb", "Fos", "Hes1"),
    Venous_EC                 = c("Ackr1", "Selp", "Sele", "Vcam1"),
    Angiogenic_EC             = c("Esm1", "Plvap")),
  Lymphatic_EC = list(
    Lymphatic_EC = c("Prox1", "Flt4", "Ccl21a")),
  Fibroblast = list(
    Inflammatory_fib        = c("Ccl2", "Ccl7", "Cxcl1", "Cxcl13", "Tnfaip6", "Fos"),
    Neural_associated_fib   = c("Nova1", "Pcolce2", "Ppp2r2b", "Sema3c"),
    ECM_remodeling_fib      = c("Abca8a", "Col3a1", "Dcn", "Lama2", "Htra3"),
    ECM_producing_fib       = c("Mfap4", "Ccdc80", "Col1a2", "Serpinf1", "Lum", "Pdgfra", "Fbln1"),
    Metabolically_active_fib= c("Gsn", "Inmt", "Pid1", "Clec3b", "Serping1"),
    Pi16_perivascular_fib   = c("Pi16", "Igfbp4", "Igfbp5", "Mfap5"),
    Vascular_associated_fib = c("Antxr2", "Pon1", "Adamts5", "Cd248"),
    Tendon_like_fib         = c("Cilp", "Thbs4", "Eln", "Aspn"),
    Wnt_inhibitory_fib      = c("Wif1", "Prg4", "Apcdd1", "S100a4", "Col1a1", "Cdh11"),
    Sdk1_adventitial_fib    = c("Sdk1", "Creb5", "Slit3", "Abi3bp", "Celf2", "Arhgap6", "Naaladl2"),
    Progenitor_fib          = c("Apod", "Myoc", "Abca6", "Gpc3", "Smoc2", "Il33", "Spp1", "Sox9")),
  Mural = list(
    Pericytes      = c("Pdgfrb", "Rgs5", "Abcc9", "Kcnj8"),
    Smooth_muscle  = c("Tagln", "Myh11", "Acta2", "Mylk", "Myl9")),
  Skeletal_muscle = list(
    Skeletal_muscle = c("Acta1", "Tnnc2", "Tnnt3", "Tnni2", "Myh4", "Myl1")),
  Satellite = list(
    Muscle_satellite = c("Pax7", "Myf5", "Calcr", "Sdc4", "Ncam1")),
  Schwann = list(
    Schwann_Lgi4_Prx = c("Lgi4", "Prx"),        # slide cluster 32
    Schwann_Plp1_Mpz = c("Plp1", "Mpz")),       # slide cluster 39
  Mesenchymal_like = list(
    Mesenchymal_like = c("Igfbp6", "Tenm2", "Sfrp5", "Bnc2", "Ptch1", "Dcn", "Mgp", "Crlf1", "Itgb4"))
)

t_genes <- c("Cd3e", "Cd3d", "Cd3g")      # per-cell T cell gate

## =============================================================================
## PART 2. Per-sample annotation
## =============================================================================
get_object <- function(s) {
  f <- file.path(out_root, s, paste0(s, "_seurat_lineage.rds"))
  if (file.exists(f)) return(readRDS(f))           # reuse clustering already done
  message(s, ": building object from Cell Ranger")
  counts <- Read10X(file.path(cr_root, samples[[s]], "filtered_feature_bc_matrix"))
  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  so <- CreateSeuratObject(counts, project = s, min.cells = 3, min.features = 200)
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
  so <- subset(so, subset = nFeature_RNA < 7000 & percent.mt < 20)
  so <- NormalizeData(so, verbose = FALSE)
  so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
  so <- ScaleData(so, verbose = FALSE)
  so <- RunPCA(so, npcs = 30, verbose = FALSE)
  so <- FindNeighbors(so, dims = 1:20, verbose = FALSE)
  so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
  so <- RunUMAP(so, dims = 1:20, verbose = FALSE)
  dir.create(file.path(out_root, s), recursive = TRUE, showWarnings = FALSE)
  saveRDS(so, f)
  so
}

add_scores <- function(so, sets, prefix) {
  sets <- lapply(sets, intersect, rownames(so))
  sets <- sets[lengths(sets) > 0]
  so <- AddModuleScore(so, features = sets, name = paste0(prefix, "_TMP_"), seed = 1)
  for (i in seq_along(sets)) {
    so[[paste0(prefix, "_", names(sets)[i])]] <- so[[paste0(prefix, "_TMP_", i)]]
    so[[paste0(prefix, "_TMP_", i)]] <- NULL
  }
  list(so = so, cols = paste0(prefix, "_", names(sets)))
}

annotate_sample <- function(s, overrides = NULL) {
  so <- get_object(s)
  sdir <- file.path(out_root, s, "annotation"); dir.create(sdir, recursive = TRUE, showWarnings = FALSE)

  ## main type per cluster (best median score; gap to 2nd best flags ambiguous clusters)
  r <- add_scores(so, main_markers, "MAIN"); so <- r$so; main_cols <- r$cols
  cl <- so@meta.data %>% group_by(seurat_clusters) %>%
    summarise(n = n(), nFeature = median(nFeature_RNA), pct_mt = round(median(percent.mt), 1),
              across(all_of(main_cols), median), .groups = "drop")
  sc <- as.matrix(cl[, main_cols])
  ord <- t(apply(sc, 1, order, decreasing = TRUE))
  cl$main_type <- sub("MAIN_", "", main_cols[ord[, 1]])
  cl$second    <- sub("MAIN_", "", main_cols[ord[, 2]])
  cl$gap       <- round(sc[cbind(seq_len(nrow(sc)), ord[, 1])] - sc[cbind(seq_len(nrow(sc)), ord[, 2])], 2)
  cl$best      <- round(sc[cbind(seq_len(nrow(sc)), ord[, 1])], 2)
  cl$flag <- case_when(cl$best < 0.1                      ~ "low_score",
                       cl$gap < 0.15                       ~ "ambiguous/doublet?",
                       cl$pct_mt > 10 | cl$nFeature < 500  ~ "low_QC?",
                       TRUE                                ~ "")
  ## manual overrides: list(cluster = "Main_type"), e.g. list("9" = "LowQC")
  if (!is.null(overrides)) for (k in names(overrides)) cl$main_type[cl$seurat_clusters == k] <- overrides[[k]]
  so$main_type <- cl$main_type[match(so$seurat_clusters, cl$seurat_clusters)]

  ## subtype within main type (per cluster)
  so$subtype <- so$main_type
  for (mt in intersect(names(sub_markers), unique(so$main_type))) {
    r <- add_scores(so, sub_markers[[mt]], "SUB"); so <- r$so
    if (length(r$cols) < 2) { so$subtype[so$main_type == mt] <- sub("SUB_", "", r$cols[1]); next }
    sc2 <- so@meta.data[so$main_type == mt, ] %>% group_by(seurat_clusters) %>%
      summarise(across(all_of(r$cols), median), .groups = "drop")
    best2 <- sub("SUB_", "", r$cols[apply(as.matrix(sc2[, r$cols]), 1, which.max)])
    idx <- so$main_type == mt
    so$subtype[idx] <- best2[match(so$seurat_clusters[idx], sc2$seurat_clusters)]
    so@meta.data[, r$cols] <- NULL
  }

  ## T cells gated per cell: immune + >= 2 of Cd3e/Cd3d/Cd3g detected
  tg <- intersect(t_genes, rownames(so))
  cd3 <- colSums(GetAssayData(so, layer = "counts")[tg, , drop = FALSE] > 0)
  so$is_T <- so$main_type == "Immune" & cd3 >= 2
  message(s, ": T cells gated = ", sum(so$is_T))

  ## outputs for review
  cl$subtype <- sapply(cl$seurat_clusters, function(k)
    names(sort(table(so$subtype[so$seurat_clusters == k]), decreasing = TRUE))[1])
  write.csv(cl %>% mutate(across(where(is.numeric), ~ round(.x, 2))),
            file.path(sdir, paste0(s, "_cluster_annotation.csv")), row.names = FALSE)

  p1 <- DimPlot(so, group.by = "main_type", label = TRUE, repel = TRUE) + ggtitle(paste(s, "- main type"))
  p2 <- DimPlot(so, group.by = "seurat_clusters", label = TRUE) + NoLegend() + ggtitle("clusters")
  ggsave(file.path(sdir, paste0(s, "_UMAP_main_type.png")), p1 | p2, width = 15, height = 6, dpi = 200)

  pan <- lapply(main_markers, function(g) intersect(g[1:min(4, length(g))], rownames(so)))
  p3 <- DotPlot(so, features = pan, group.by = "seurat_clusters") + RotatedAxis() +
    ggtitle(paste(s, "- pan markers by cluster")) + theme(strip.text.x = element_text(angle = 90, size = 7))
  ggsave(file.path(sdir, paste0(s, "_DotPlot_pan_by_cluster.png")), p3, width = 18, height = 6, dpi = 200)

  for (mt in intersect(names(sub_markers), unique(so$main_type))) {
    sub_o <- subset(so, subset = main_type == mt)
    if (ncol(sub_o) < 20) next
    feats <- lapply(sub_markers[[mt]], intersect, rownames(sub_o))
    p <- DotPlot(sub_o, features = feats, group.by = "seurat_clusters") + RotatedAxis() +
      ggtitle(paste(s, "-", mt)) + theme(strip.text.x = element_text(angle = 90, size = 7))
    ggsave(file.path(sdir, paste0(s, "_DotPlot_", mt, ".png")), p,
           width = max(6, 0.35 * length(unlist(feats)) + 3), height = 5, dpi = 200)
  }
  saveRDS(so, file.path(out_root, s, paste0(s, "_seurat_annotated.rds")))
  print(cl %>% select(seurat_clusters, n, main_type, subtype, best, second, gap, flag, nFeature, pct_mt), n = 60)
  invisible(so)
}

## --- 2a. First pass on every sample (~2-5 min each) ---
for (s in names(samples)) {
  message("\n==================== ANNOTATE ", s, " ====================")
  tryCatch(annotate_sample(s), error = function(e) message("ERROR ", s, ": ", e$message))
}
## >>> REVIEW <sample>/annotation/*.csv and *.png. Fix clusters below, then run 2b.

## --- 2b. Corrections (cluster = main type), then re-annotate only those samples ---
## Use "LowQC" or "Doublet" for clusters to exclude. Examples from earlier checks:
overrides <- list(
  FT  = list("9" = "LowQC"),
  M2  = list("9" = "Doublet"),
  F2  = list("3" = "Doublet"),
  MC  = list("18" = "LowQC"),
  FN1 = list("4" = "LowQC", "11" = "LowQC", "13" = "LowQC"),
  MN1 = list("2" = "LowQC")
)
for (s in names(overrides)) annotate_sample(s, overrides[[s]])

## =============================================================================
## PART 3. CopyKAT: test = epithelium, reference = T cells of the same sample
## =============================================================================
min_ref  <- 50       # minimum T cells; below this, fall back to all non-neutrophil immune cells
max_test <- 2000
max_ref  <- 1000
run_order <- c("FN1", "MN1", "FT", "MTT", "F1", "F2", "M1", "M2", "FC", "MC", "FN2", "MN2")

for (s in run_order) {
  f <- file.path(out_root, s, paste0(s, "_seurat_annotated.rds"))
  if (!file.exists(f)) { message(s, ": not annotated, skipping"); next }
  sdir <- file.path(out_root, s, "vs_Tcells"); dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  nm <- paste0(s, "_vsT")
  if (file.exists(file.path(sdir, paste0(nm, "_copykat_prediction.txt")))) { message(s, ": done already"); next }

  so <- readRDS(f)
  set.seed(1)
  test_cells <- colnames(so)[so$main_type == "Epithelial"]
  ref_cells  <- colnames(so)[so$is_T]
  ref_label  <- "T cells"
  if (length(ref_cells) < min_ref) {
    ref_cells <- colnames(so)[so$main_type == "Immune" & so$subtype != "Neutrophils"]
    ref_label <- "immune (T + macrophages + mast; T cells too few)"
  }
  if (length(test_cells) > max_test) test_cells <- sample(test_cells, max_test)
  if (length(ref_cells)  > max_ref)  ref_cells  <- sample(ref_cells,  max_ref)
  message("\n==== CopyKAT ", s, ": test = ", length(test_cells), " epithelial | reference = ",
          length(ref_cells), " ", ref_label, "  (", format(Sys.time(), "%H:%M"), ")")
  if (length(test_cells) < 30 || length(ref_cells) < 30) { message("Too few cells, skipping"); next }

  write.csv(data.frame(cell = c(test_cells, ref_cells),
                       role = rep(c("test", "reference"), c(length(test_cells), length(ref_cells))),
                       subtype = so$subtype[match(c(test_cells, ref_cells), colnames(so))]),
            file.path(sdir, paste0(nm, "_cells_used.csv")), row.names = FALSE)

  rawmat <- as.matrix(GetAssayData(so, layer = "counts")[, c(test_cells, ref_cells)])
  old_wd <- setwd(sdir); t0 <- Sys.time()
  try(copykat(rawmat = rawmat, id.type = "S", genome = "mm10",
              ngene.chr = 5, win.size = 25, KS.cut = 0.1, distance = "euclidean",
              norm.cell.names = ref_cells, output.seg = "FALSE",
              plot.genes = "FALSE", sam.name = nm, n.cores = 1))
  setwd(old_wd)
  message(s, " time: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  rm(rawmat, so); gc(); graphics.off()
}

## =============================================================================
## PART 4. Summary
## =============================================================================
summ <- list()
for (s in run_order) {
  sdir <- file.path(out_root, s, "vs_Tcells"); nm <- paste0(s, "_vsT")
  pf <- file.path(sdir, paste0(nm, "_copykat_prediction.txt"))
  uf <- file.path(sdir, paste0(nm, "_cells_used.csv"))
  if (!file.exists(pf)) next
  pred <- read.delim(pf); used <- read.csv(uf)
  pred$role    <- used$role[match(pred$cell.names, used$cell)]
  pred$subtype <- used$subtype[match(pred$cell.names, used$cell)]
  miss <- is.na(pred$role)                      # barcodes with "-" changed to "."
  pred$role[miss]    <- used$role[match(pred$cell.names[miss], make.names(used$cell))]
  pred$subtype[miss] <- used$subtype[match(pred$cell.names[miss], make.names(used$cell))]
  write.csv(as.data.frame.matrix(table(pred$subtype, pred$copykat.pred)),
            file.path(sdir, paste0(nm, "_calls_by_subtype.csv")))
  summ[[s]] <- data.frame(
    sample = s,
    n_epithelial = sum(pred$role == "test", na.rm = TRUE),
    pct_aneuploid_epithelial = round(100 * mean(pred$copykat.pred[pred$role %in% "test"] == "aneuploid"), 1),
    pct_diploid_reference    = round(100 * mean(pred$copykat.pred[pred$role %in% "reference"] == "diploid"), 1))
}
summ <- bind_rows(summ); print(summ)
write.csv(summ, file.path(out_root, "summary_vs_Tcells_all_samples.csv"), row.names = FALSE)
## Read the normal samples (FN1, MN1) FIRST: if their epithelium is mostly "aneuploid",
## the T-cell reference produces an expression artifact (as the immune+EC reference did before).