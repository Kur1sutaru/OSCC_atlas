## ============================================================================
## CopyKAT: tumor epithelium (FT, MTT) vs normal epithelium (FN1, MN1)
## Mouse tongue - full workflow
##   Part A. Fix lineage labels (already done once; safe to re-run)
##   Part B. Build the combined matrix and run CopyKAT (skipped if results exist)
##   Part C. Summary per sample
##   Part D. Heatmap + chromosome profile + FT vs MTT correlation
## Requires per-sample objects: <out_root>/<sample>/<sample>_seurat_lineage.rds
## ============================================================================

library(Seurat)
library(copykat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)    # BiocManager::install("ComplexHeatmap")
library(circlize)

out_root <- "D:/cross_species_paper/copykat_out"
sdir     <- file.path(out_root, "TUMORS_vs_NORMAL_EPI")
sam_name <- "TumVsNormEpi"
dir.create(sdir, showWarnings = FALSE)

tumors  <- c("FT", "MTT")
normals <- c("FN1", "MN1")

## ---------------------------------------------------------------------------
## Part A. Lineage fixes (manual, from the DotPlots + QC)
## ---------------------------------------------------------------------------
fix <- list(
  FT  = list(epi = c("0", "1", "2", "3", "8", "13"), lowqc = "9",
             other = list(Fibroblast = c("6", "10"), Schwann = "14")),
  MTT = list(epi = c("1", "3", "6", "8"), lowqc = character(0),
             other = list(Immune = c("0", "2", "9", "10", "11"), Endothelial = c("5", "12"),
                          Fibroblast = c("4", "7"), Schwann = "13")),
  FN1 = list(epi = c("9", "14"), supra = "10", lowqc = c("4", "11", "13")),
  MN1 = list(epi = c("6", "10"), supra = "11", lowqc = "2")
)

for (s in names(fix)) {
  f <- file.path(out_root, s, paste0(s, "_seurat_lineage.rds"))
  if (!file.exists(f)) { message("Missing object: ", f); next }
  o  <- readRDS(f)
  fx <- fix[[s]]
  if (!is.null(fx$other)) for (lab in names(fx$other))
    o$lineage[o$seurat_clusters %in% fx$other[[lab]]] <- lab
  o$lineage[o$seurat_clusters %in% fx$epi] <- "Epithelial"
  if (!is.null(fx$supra)) o$lineage[o$seurat_clusters %in% fx$supra] <- "Epithelial_suprabasal"
  o$lineage[o$seurat_clusters %in% fx$lowqc] <- "LowQC"
  message(s, ": ", paste(names(table(o$lineage)), table(o$lineage), sep = "=", collapse = "  "))
  saveRDS(o, f)
}

## ---------------------------------------------------------------------------
## Part B. Combined matrix + CopyKAT
## ---------------------------------------------------------------------------
pred_file <- file.path(sdir, paste0(sam_name, "_copykat_prediction.txt"))

if (!file.exists(pred_file)) {
  load_epi <- function(s, n_max = 1500) {
    o <- readRDS(file.path(out_root, s, paste0(s, "_seurat_lineage.rds")))
    o <- subset(o, subset = lineage == "Epithelial")
    set.seed(1)
    if (ncol(o) > n_max) o <- o[, sample(colnames(o), n_max)]
    m <- GetAssayData(o, layer = "counts")   # slot = "counts" on Seurat v4
    colnames(m) <- paste0(s, "_", colnames(m))   # prefix avoids barcode clashes
    m
  }
  mats <- lapply(setNames(c(tumors, normals), c(tumors, normals)), load_epi)
  print(sapply(mats, ncol))
  
  genes     <- Reduce(intersect, lapply(mats, rownames))
  rawmat    <- as.matrix(do.call(cbind, lapply(mats, function(m) m[genes, ])))
  ref_cells <- unlist(lapply(mats[normals], colnames))
  message("Matrix: ", nrow(rawmat), " genes x ", ncol(rawmat), " cells; reference = ", length(ref_cells))
  
  old_wd <- setwd(sdir)
  t0 <- Sys.time()
  try(copykat(rawmat = rawmat, id.type = "S", genome = "mm10",
              ngene.chr = 5, win.size = 25, KS.cut = 0.1, distance = "euclidean",
              norm.cell.names = ref_cells, output.seg = "FALSE",
              plot.genes = "FALSE", sam.name = sam_name, n.cores = 1))
  setwd(old_wd)
  message("CopyKAT time: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  rm(rawmat); gc(); graphics.off()
} else {
  message("CopyKAT results already exist - skipping the run")
}

## ---------------------------------------------------------------------------
## Part C. Summary per sample
## ---------------------------------------------------------------------------
pred <- read.delim(pred_file)
pred$sample <- sub("_.*", "", pred$cell.names)
tab <- table(pred$sample, pred$copykat.pred)
print(tab)
print(round(100 * prop.table(tab, 1), 1))
write.csv(as.data.frame.matrix(tab), file.path(sdir, "summary_by_sample_counts.csv"))
write.csv(round(100 * as.data.frame.matrix(prop.table(tab, 1)), 1),
          file.path(sdir, "summary_by_sample_percent.csv"))

## ---------------------------------------------------------------------------
## Part D. Heatmap + chromosome profile
## ---------------------------------------------------------------------------
## CopyKAT writes a bin-level file (chrom, chrompos, abspos) and a gene-level file
## (abspos, chromosome_name, ..., mgi_symbol, band). Use whichever is there.
print(list.files(sdir))
cna_file <- file.path(sdir, paste0(sam_name, "_copykat_CNA_results.txt"))
cna <- read.delim(cna_file, check.names = FALSE)

meta_cols <- intersect(colnames(cna),
                       c("chrom", "chrompos", "abspos", "chromosome_name", "start_position",
                         "end_position", "ensembl_gene_id", "mgi_symbol", "band"))
cell_cols <- setdiff(colnames(cna), meta_cols)
message("Metadata columns: ", paste(meta_cols, collapse = ", "), " | cells: ", length(cell_cols))

chr_col <- if ("chrom" %in% meta_cols) "chrom" else "chromosome_name"
cna <- cna[order(cna$abspos), ]
chr_chr <- as.character(cna[[chr_col]])
chr_chr[chr_chr == "X"] <- "23"
cna$chr_num <- as.numeric(chr_chr)

mat <- t(as.matrix(cna[, cell_cols]))
storage.mode(mat) <- "numeric"
print(summary(as.vector(mat)))

smp   <- sub("_.*", "", rownames(mat))
print(table(smp))
chrom <- factor(cna$chr_num, levels = sort(unique(cna$chr_num)))

## D1. Heatmap (cells x genome), split by sample
lim <- quantile(abs(mat), 0.99, na.rm = TRUE)
smp_cols <- c(FN1 = "#2E86C1", MN1 = "#85C1E9", FT = "#C0392B", MTT = "#F1948A")
ann <- rowAnnotation(sample = smp, col = list(sample = smp_cols), show_annotation_name = FALSE)
ht <- Heatmap(mat, name = "CNA",
              col = colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B")),
              cluster_columns = FALSE, column_split = chrom, column_gap = unit(0.5, "mm"),
              column_title_gp = gpar(fontsize = 8), show_column_names = FALSE,
              row_split = factor(smp, levels = c(normals, tumors)),
              cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE,
              left_annotation = ann, use_raster = TRUE)
png(file.path(sdir, paste0(sam_name, "_heatmap.png")), width = 3200, height = 2000, res = 200)
draw(ht)
dev.off()

## D2. Mean CNA per chromosome and sample (chrX = 23 excluded: different sexes)
prof <- sapply(split(seq_len(nrow(cna)), cna$chr_num), function(r)
  tapply(colMeans(as.matrix(cna[r, cell_cols])), smp, mean))
prof <- prof[, colnames(prof) != "23", drop = FALSE]
print(round(prof, 3))
write.csv(round(prof, 4), file.path(sdir, paste0(sam_name, "_mean_CNA_per_chromosome.csv")))

df <- data.frame(sample = rep(rownames(prof), times = ncol(prof)),
                 chrom  = rep(as.numeric(colnames(prof)), each = nrow(prof)),
                 mean_cna = as.vector(prof))
p <- ggplot(df, aes(factor(chrom), mean_cna, colour = sample, group = sample)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_line() + geom_point() +
  scale_colour_manual(values = smp_cols) +
  labs(x = "Chromosome", y = "Mean CNA (epithelial cells)",
       title = "Tumor vs normal epithelium - CopyKAT") + theme_classic()
ggsave(file.path(sdir, paste0(sam_name, "_chrom_profile.png")), p, width = 9, height = 4, dpi = 200)

## D3. Do the two tumors share the same pattern?
message("Correlation FT vs MTT (per-chromosome profile): ",
        round(cor(prof["FT", ], prof["MTT", ]), 2))

## D4. Chromosomes with the largest shift in BOTH tumors
shift <- data.frame(chrom = colnames(prof),
                    FT = prof["FT", ], MTT = prof["MTT", ],
                    normal_mean = colMeans(prof[normals, , drop = FALSE]))
shift$tumor_mean <- (shift$FT + shift$MTT) / 2
shift$same_direction <- sign(shift$FT) == sign(shift$MTT)
shift <- shift[order(-abs(shift$tumor_mean)), ]
print(head(shift, 10))
write.csv(shift, file.path(sdir, paste0(sam_name, "_chromosome_shifts.csv")), row.names = FALSE)


## 1. FT's 16 diploid cells: are they the normal-like basal cluster 13 (Krt5+)?
ft   <- readRDS("D:/cross_species_paper/copykat_out/FT/FT_seurat_lineage.rds")
pred <- read.delim("D:/cross_species_paper/copykat_out/TUMORS_vs_NORMAL_EPI/TumVsNormEpi_copykat_prediction.txt")
ft_pred <- pred[grepl("^FT_", pred$cell.names), ]
ft_pred$barcode <- sub("^FT_", "", ft_pred$cell.names)
ft_pred$cluster <- ft$seurat_clusters[match(ft_pred$barcode, colnames(ft))]
table(ft_pred$cluster, ft_pred$copykat.pred)

## 2. Does the chr11 "gain" survive without the keratin cluster?
cna <- read.delim("D:/cross_species_paper/copykat_out/TUMORS_vs_NORMAL_EPI/TumVsNormEpi_copykat_CNA_results.txt",
                  check.names = FALSE)
meta <- intersect(colnames(cna), c("chrom","chrompos","abspos","chromosome_name","start_position",
                                   "end_position","ensembl_gene_id","mgi_symbol","band"))
cells <- setdiff(colnames(cna), meta)
c11 <- cna[cna$chromosome_name == 11, ]
krt <- grepl("^Krt", c11$mgi_symbol)
smp <- sub("_.*", "", cells)
sapply(list(all_chr11 = rep(TRUE, nrow(c11)), without_keratins = !krt), function(k)
  tapply(colMeans(as.matrix(c11[k, cells])), smp, mean))

message("Done. Files in ", normalizePath(sdir))