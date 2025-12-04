## analysis for USASP abstract Jaclyn
setwd("C:/Users/crist/Desktop/sex_diff_analysis")

library(readxl)
df <- read_excel("RPKM over 5 AllComparisonTable_with_RPKM_NormReadCounts_exclude2Samples.xls")
colnames(df)

# Male-specific tumor genes = (MTvsMN) − (FTvsFN)
# Female-specific tumor genes = (FTvsFN) − (MTvsMN)
# Shared tumor genes = intersection
# This is AMAZING for:
# Biological sex divergence
# Evolutionary/functional interpretation
# CellChat filtering
# Visualizations (UPSET, Venn, radar)


library(dplyr)
library(readr)
library(stringr)

# 1) Read the four tables -----------------------------------------
ft_05 <- read_csv("FTvsFNFC_0_5.csv")    # female tumor vs female normal
mt_05 <- read_csv("MTvsMNFC_0_5.csv")    # male tumor vs male normal

# 2) Try to standardize column names for logFC ---------------------
# Adapt this if your column name is different (e.g., "LOG FC", "log2FoldChange", "FC")
guess_logfc_col <- function(df) {
  cn <- colnames(df)
  # look for something with "log" and "fc" or just "LOG FC"
  log_like <- cn[str_detect(tolower(cn), "log") & str_detect(tolower(cn), "fc")]
  if (length(log_like) == 0) {
    # maybe there is "FC" only and all FC > 0.5, so we just treat as up
    return(NULL)
  } else {
    return(log_like[1])
  }
}

logfc_ft <- guess_logfc_col(ft_05)
logfc_mt <- guess_logfc_col(mt_05)

# 3) Label direction: up vs down genes -----------------------------
# If you *only* have FC>0.5 in these tables (no downregulated), you can skip
# the down sets and just treat these as "tumor-up" signatures.

if (!is.null(logfc_ft)) {
  ft_05 <- ft_05 %>%
    rename(logFC = !!logfc_ft) %>%
    mutate(direction = ifelse(logFC > 0, "up", "down"))
} else {
  ft_05 <- ft_05 %>% mutate(direction = "up")
}

if (!is.null(logfc_mt)) {
  mt_05 <- mt_05 %>%
    rename(logFC = !!logfc_mt) %>%
    mutate(direction = ifelse(logFC > 0, "up", "down"))
} else {
  mt_05 <- mt_05 %>% mutate(direction = "up")
}

# 4) Extract SYMBOL-based lists ------------------------------------
# Use YOUR gene symbol column name here ("Symbol" in your screenshot)
sym_col <- "Symbol"

female_up   <- ft_05 %>% filter(direction == "up")   %>% pull(!!sym_col) %>% unique()
female_down <- ft_05 %>% filter(direction == "down") %>% pull(!!sym_col) %>% unique()

male_up     <- mt_05 %>% filter(direction == "up")   %>% pull(!!sym_col) %>% unique()
male_down   <- mt_05 %>% filter(direction == "down") %>% pull(!!sym_col) %>% unique()

# 5) Sex-specific & shared signatures ------------------------------

shared_up         <- intersect(female_up, male_up)
shared_down       <- intersect(female_down, male_down)

male_specific_up   <- setdiff(male_up, female_up)
female_specific_up <- setdiff(female_up, male_up)

male_specific_down   <- setdiff(male_down, female_down)
female_specific_down <- setdiff(female_down, male_down)

# Save for record / enrichment analysis
write_csv(data.frame(Symbol = shared_up),         "shared_up_genes.csv")
write_csv(data.frame(Symbol = shared_down),       "shared_down_genes.csv")
write_csv(data.frame(Symbol = male_specific_up),  "male_specific_up_genes.csv")
write_csv(data.frame(Symbol = female_specific_up),"female_specific_up_genes.csv")
write_csv(data.frame(Symbol = male_specific_down),"male_specific_down_genes.csv")
write_csv(data.frame(Symbol = female_specific_down),"female_specific_down_genes.csv")

#FN1 cell type annotation
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
FN1.markers <- FindAllMarkers(FN1, only.pos = TRUE)
FN1.markers %>%
  group_by(seurat_cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(FN1.markers, 
          file = "FN1_markers_log2FC_gt1.csv",
          row.names = FALSE)


new.cluster.ids <- c(
  "0"  = "Endothelial",
  "1"  = "Fibro",
  "2"  = "Epithelial",
  "3"  = "Endothelial",
  "4"  = "Endothelial",
  "5"  = "Endothelial",
  "6"  = "Pericytes",
  "7"  = "Endothelial",
  "8"  = "Epithelial",
  "9"  = "Endothelial",
  "10" = "Immune",
  "11" = "Fibro",
  "12" = "Keratinocytes",
  "13" = "Schwann",
  "14" = "Endothelial",
  "15" = "Skeletal Muscle",
  "16" = "Immune",
  "17" = "Immune",
  "18" = "Skeletal Muscle"
)
FN1$seurat_clusters <- as.character(FN1$seurat_clusters)
FN1$celltype_annot <- new.cluster.ids[ FN1$seurat_clusters ]
DimPlot(FN1, group.by = "celltype_annot", label = TRUE, repel = TRUE)

#FT
#FT cell type annotation
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
FT.markers <- FindAllMarkers(FT, only.pos = TRUE)
FT.markers %>%
  group_by(seurat_cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(FT.markers, 
          file = "FT_markers_log2FC_gt1.csv",
          row.names = FALSE)


new.cluster.ids <- c(
  "0"  = "Fibroblast",
  "1"  = "Neutrophils",
  "2"  = "Epithelial",
  "3"  = "Epithelial",
  "4"  = "Epithelial",
  "5"  = "Keratinocytes",
  "6"  = "Endothelial",
  "7"  = "Fibroblast",
  "8"  = "Macrophages",
  "9"  = "Epithelial",
  "10" = "Pericytes",
  "11" = "T cell",
  "12" = "Endothelial",
  "13" = "Keratinocytes",
  "14" = "Schwann"
)
FT$seurat_clusters <- as.character(FT$seurat_clusters)
FT$celltype_annot <- new.cluster.ids[ FT$seurat_clusters ]
DimPlot(FT, group.by = "celltype_annot", label = TRUE, repel = TRUE)
saveRDS(FT, "FT_annot.rds")



###MN1
MN1.markers <- FindAllMarkers(MN1, only.pos = TRUE)
MN1.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(MN1.markers, 
          file = "MN1_markers_log2FC_gt1.csv",
          row.names = FALSE)


new.cluster.ids <- c(
  "0"  = "Endothelial",
  "1"  = "Endothelial",
  "2"  = "Endothelial",
  "3"  = "Endothelial",
  "4"  = "Endothelial",
  "5"  = "Keratinocytes",
  "6"  = "Mural - SMC - Pericytes",
  "7"  = "Fibroblast",
  "8"  = "Endothelial",
  "9"  = "Keratinocytes",
  "10" = "Endothelial",
  "11" = "Endothelial",
  "12" = "Fibroblast",
  "13" = "Skeletal Muscle",
  "14" = "Immune",
  "15" = "Endothelial",
  "16" = "Skeletal Muscle"
)

MN1$seurat_clusters <- as.character(MN1$seurat_clusters)

MN1$celltype_annot <- unname(
  new.cluster.ids[ MN1$seurat_clusters ]
)
setdiff(unique(MN1$seurat_clusters), names(new.cluster.ids))


DimPlot(MN1, group.by = "celltype_annot", label = TRUE, repel = TRUE)
saveRDS(MN1, "MN1_annot.rds")



###MT
MT.markers <- FindAllMarkers(MT, only.pos = TRUE)
MT.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv(MT.markers, 
          file = "MT_markers_log2FC_gt1.csv",
          row.names = FALSE)


new.cluster.ids <- c(
  "0"  = "Neutrophil",
  "1"  = "Epithelial",
  "2"  = "Macrophages",
  "3"  = "Dendritic cells",
  "4"  = "Epithelial",
  "5"  = "Endothelial",
  "6"  = "Epithelial",
  "7"  = "Fibroblast",
  "8"  = "Mural cells",
  "9"  = "Epithelial",
  "10" = "T cells",
  "11" = "Macrophages",
  "12" = "Endothelial",
  "13" = "Schwann cells",
  "14" = "Dendritic cells"
)

MT$seurat_clusters <- as.character(MT$seurat_clusters)

MT$celltype_annot <- unname(
  new.cluster.ids[ MT$seurat_clusters ]
)
setdiff(unique(MT$seurat_clusters), names(new.cluster.ids))


DimPlot(MT, group.by = "celltype_annot", label = TRUE, repel = TRUE)
saveRDS(MT, "MT_annot.rds")


# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Load annotated datasets
FN <- readRDS("FN1_annot.rds")  # Female Normal
FT <- readRDS("FT_annot.rds")  # Female Tumor
MN <- readRDS("MN1_annot.rds") # Male Normal
MT <- readRDS("MT_annot.rds")  # Male Tumor

# Add sample metadata
FN$sample <- "Female_Normal"
FT$sample <- "Female_Tumor"
MN$sample <- "Male_Normal"
MT$sample <- "Male_Tumor"



cluster_annotations <- c(
  "0"  = "Endothelial",
  "1"  = "Fibroblast",
  "2"  = "Epithelial",
  "3"  = "Endothelial",
  "4"  = "Endothelial",
  "5"  = "Endothelial",
  "6"  = "Pericytes",
  "7"  = "Endothelial",
  "8"  = "Epithelial",
  "9"  = "Endothelial",
  "10" = "Immune",
  "11" = "Fibroblast",
  "12" = "Keratinocytes",
  "13" = "Schwann",
  "14" = "Endothelial",
  "15" = "Skeletal Muscle",
  "16" = "Immune",
  "17" = "Immune",
  "18" = "Skeletal Muscle"
)

# Apply mapping
FN <- RenameIdents(FN, cluster_annotations)

levels(FN)   # should now show "Endothelial", "Fibroblast", etc.



cluster_annotations <- c(
  "0"  = "Fibroblast",
  "1"  = "Neutrophils",
  "2"  = "Epithelial",
  "3"  = "Epithelial",
  "4"  = "Epithelial",
  "5"  = "Keratinocytes",
  "6"  = "Endothelial",
  "7"  = "Fibroblast",
  "8"  = "Macrophages",
  "9"  = "Epithelial",
  "10" = "Pericytes",
  "11" = "T cell",
  "12" = "Endothelial",
  "13" = "Keratinocytes",
  "14" = "Schwann"
)

# Apply mapping
FT <- RenameIdents(FT, cluster_annotations)

levels(FT)

cluster_annotations <- c(
  "0"  = "Endothelial",
  "1"  = "Endothelial",
  "2"  = "Endothelial",
  "3"  = "Endothelial",
  "4"  = "Endothelial",
  "5"  = "Keratinocytes",
  "6"  = "Mural - SMC - Pericytes",
  "7"  = "Fibroblast",
  "8"  = "Endothelial",
  "9"  = "Keratinocytes",
  "10" = "Endothelial",
  "11" = "Endothelial",
  "12" = "Fibroblast",
  "13" = "Skeletal Muscle",
  "14" = "Immune",
  "15" = "Endothelial",
  "16" = "Skeletal Muscle"
)

# Apply mapping
MN <- RenameIdents(MN, cluster_annotations)

levels(MN)

cluster_annotations <- c(
  "0"  = "Neutrophil",
  "1"  = "Epithelial",
  "2"  = "Macrophages",
  "3"  = "Dendritic cells",
  "4"  = "Epithelial",
  "5"  = "Endothelial",
  "6"  = "Epithelial",
  "7"  = "Fibroblast",
  "8"  = "Mural cells",
  "9"  = "Epithelial",
  "10" = "T cells",
  "11" = "Macrophages",
  "12" = "Endothelial",
  "13" = "Schwann cells",
  "14" = "Dendritic cells"
)


# Apply mapping
MT <- RenameIdents(MT, cluster_annotations)

levels(MT)

# Combine all
combined <- merge(FN, y = list(FT, MN, MT))



standard_labels <- c(
  "Endothelial" = "Endothelial",
  "Mural cells" = "Mural cells",
  "Mural - SMC - Pericytes" = "Mural cells",
  "Pericytes" = "Mural cells",
  "Fibro" = "Fibroblast",
  "Fibroblast" = "Fibroblast",
  "Cancer-associated fibroblasts" = "Fibroblast",
  "Epithelial" = "Epithelial",
  "Keratinocytes" = "Keratinocytes",
  "Schwann" = "Schwann cells",
  "Schwann cells" = "Schwann cells",
  "Skeletal Muscle" = "Skeletal muscle",
  "Immune" = "Immune",
  "T cell" = "T cells",
  "T cells" = "T cells",
  "Macrophages" = "Macrophages",
  "Dendritic cells" = "Dendritic cells",
  "Neutrophil" = "Neutrophils",
  "Neutrophils" = "Neutrophils"
)

# Apply to each object
FN <- RenameIdents(FN, standard_labels)
FT <- RenameIdents(FT, standard_labels)
MN <- RenameIdents(MN, standard_labels)
MT <- RenameIdents(MT, standard_labels)



# Extract identities into meta.data
FN$cell_type <- Idents(FN)
FT$cell_type <- Idents(FT)
MN$cell_type <- Idents(MN)
MT$cell_type <- Idents(MT)
# Combine all
combined <- merge(FN, y = list(FT, MN, MT))

# Compute proportions
prop_df <- combined@meta.data %>%
  group_by(sample, cell_type) %>%
  summarise(count = n()) %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

ggplot(prop_df, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal(base_size = 14) +
  labs(title = "Cell Type Proportions by Sample",
       x = "Sample",
       y = "Proportion",
       fill = "Cell Type") +
  scale_fill_brewer(palette = "Set3")


  # Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Read tables (replace with your actual file paths)
MTvsMN <- read.csv("MTvsMNFC_1_5.csv")
FTvsFN <- read.csv("FTvsFNFC_0.5.csv.csv")
FT_FN  <- read.csv("FT_FNFC_1_5.csv")
MTvsMN2 <- read.csv("MTvsMNFC_0_5.csv.csv")

# Extract top 20 genes from each table
top_MTvsMN  <- head(MTvsMN$Symbol, 5)
top_FTvsFN  <- head(FTvsFN$Symbol, 5)
top_FT_FN   <- head(FT_FN$Symbol, 5)
top_MTvsMN2 <- head(MTvsMN2$Symbol, 5)

# Combine into one gene list
top_genes <- unique(c(top_MTvsMN, top_FTvsFN, top_FT_FN, top_MTvsMN2))

DotPlot(combined, features = top_genes, group.by = "sample") +
  RotatedAxis() +
  theme_minimal(base_size = 12) +
  labs(title = "Top 5 DE Genes DotPlot",
       x = "Genes",
       y = "Sample Groups")
