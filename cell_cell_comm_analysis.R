# ============================================================
# 1) Add sex column to OSCC (scRNA tumor object)
#    FT/FT2/FT3 = Female; MT2/MT3/MTT = Male
# ============================================================
OSCC$sex <- NA_character_
OSCC$sex[OSCC$orig.ident %in% c("FT","FT2","FT3")] <- "Female"
OSCC$sex[OSCC$orig.ident %in% c("MT2","MT3","MTT")] <- "Male"
OSCC$sex <- factor(OSCC$sex, levels = c("Female","Male"))

table(OSCC$sex, useNA="ifany")
table(OSCC$orig.ident, OSCC$sex, useNA="ifany")

# ============================================================
# 2) Create neur_tumor by excluding FN*/MN* normals from neur
#    and add sex column (FT* Female; MT* Male)
# ============================================================
stopifnot("sample_id" %in% colnames(neur@meta.data))

neur$normal_bulk <- grepl("_Sample_FN|_Sample_MN", neur$sample_id)
table(neur$normal_bulk)

neur_tumor <- subset(neur, subset = normal_bulk == FALSE)

# Add sex for tumor-only bulk neurons
neur_tumor$sex <- ifelse(
  grepl("_Sample_FT", neur_tumor$sample_id), "Female",
  ifelse(grepl("_Sample_MT", neur_tumor$sample_id), "Male", NA)
)
neur_tumor$sex <- factor(neur_tumor$sex, levels = c("Female","Male"))

table(neur_tumor$sample_id)
table(neur_tumor$sex, useNA="ifany")

# Ensure consistent grouping columns for CellChat later
neur_tumor$celltype <- "Neurons"

# Make cell names unique before merge
neur_tumor <- RenameCells(neur_tumor, add.cell.id = "NEURONBULK")

# ============================================================
# 3) Harmonize metadata columns (optional but recommended)
#    Make sure BOTH objects have: celltype, sex
# ============================================================
stopifnot("celltype" %in% colnames(OSCC@meta.data))
stopifnot("celltype" %in% colnames(neur_tumor@meta.data))
stopifnot("sex" %in% colnames(OSCC@meta.data))
stopifnot("sex" %in% colnames(neur_tumor@meta.data))

# If OSCC has extra levels, unify factor levels
sex_levels <- c("Female","Male")
OSCC$sex <- factor(as.character(OSCC$sex), levels = sex_levels)
neur_tumor$sex <- factor(as.character(neur_tumor$sex), levels = sex_levels)

# ============================================================
# 4) Merge: OSCC + tumor-associated bulk neurons
# ============================================================
combined_tumor <- merge(OSCC, neur_tumor)

# Sanity checks
table(combined_tumor$sex, useNA="ifany")
table(combined_tumor$sex, combined_tumor$celltype, useNA="ifany")

# ============================================================
# UPDATED FULL SNIPPET (fixes CellChat factor-level error)
# Female vs Male CellChat with Seurat v5 layered RNA union
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(CellChat)
  library(dplyr)
})

stopifnot(inherits(combined_tumor, "Seurat"))
stopifnot(all(c("sex","celltype") %in% colnames(combined_tumor@meta.data)))

combined_tumor$sex <- factor(as.character(combined_tumor$sex), levels = c("Female","Male"))
combined_tumor$celltype <- as.character(combined_tumor$celltype)  # make clean before split

# -----------------------------
# 1) Helper: merge layered matrices by UNION of genes
# -----------------------------
get_layer_matrix_union <- function(seu, assay="RNA", prefer_prefix="data") {
  DefaultAssay(seu) <- assay
  lay <- Layers(seu[[assay]])
  layers_use <- lay[grepl(paste0("^", prefer_prefix), lay)]
  if (length(layers_use) == 0) stop("No layers starting with '", prefer_prefix, "' in assay ", assay)
  
  mats <- lapply(layers_use, function(L) {
    m <- LayerData(seu, assay = assay, layer = L)
    if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
    m
  })
  
  all_genes <- Reduce(union, lapply(mats, rownames))
  
  mats2 <- lapply(mats, function(m) {
    missing <- setdiff(all_genes, rownames(m))
    if (length(missing) > 0) {
      zeros <- Matrix(0, nrow = length(missing), ncol = ncol(m), sparse = TRUE)
      rownames(zeros) <- missing
      colnames(zeros) <- colnames(m)
      m <- rbind(m, zeros)
    }
    m[all_genes, , drop = FALSE]
  })
  
  mat <- do.call(Matrix::cbind2, mats2)
  mat <- mat[, colnames(seu), drop = FALSE]
  mat
}

# -----------------------------
# 2) Helper: make CellChat-safe grouping factor
# -----------------------------
sanitize_group_factor <- function(vec) {
  # vec is a character or factor
  v <- as.character(vec)
  v[is.na(v) | v == ""] <- "Unknown"
  f <- factor(v)
  f <- droplevels(f)
  f
}

# -----------------------------
# 3) Helper: run CellChat (Seurat v5 safe + factor fix)
# -----------------------------
run_cellchat_v5 <- function(seu,
                            group_col = "celltype",
                            assay = "RNA",
                            species = c("mouse","human"),
                            use_raw = FALSE,
                            min_cells = 1) {
  species <- match.arg(species)
  stopifnot(group_col %in% colnames(seu@meta.data))
  
  # IMPORTANT: drop unused levels / sanitize group labels
  seu[[group_col]] <- sanitize_group_factor(seu[[group_col]][,1])
  
  # Extract expression matrix
  mat <- get_layer_matrix_union(
    seu,
    assay = assay,
    prefer_prefix = ifelse(use_raw, "counts", "data")
  )
  
  # Align metadata to matrix columns
  meta <- seu@meta.data[colnames(mat), , drop = FALSE]
  stopifnot(identical(rownames(meta), colnames(mat)))
  
  # Ensure meta group col is also dropped (CellChat uses meta too)
  meta[[group_col]] <- sanitize_group_factor(meta[[group_col]])
  
  # Create CellChat
  cc <- createCellChat(object = mat, meta = meta, group.by = group_col)
  
  # Double-safety: drop idents levels
  cc@idents <- droplevels(cc@idents)
  
  # DB
  cc@DB <- if (species == "mouse") CellChatDB.mouse else CellChatDB.human
  
  # Robust signaling subset
  lr.genes <- unique(c(cc@DB$interaction$ligand, cc@DB$interaction$receptor))
  gene.use <- intersect(rownames(cc@data), lr.genes)
  if (length(gene.use) < 10) stop("Too few LR genes overlap (", length(gene.use), "). Check gene symbols/species.")
  cc@data.signaling <- cc@data[gene.use, , drop = FALSE]
  
  # Standard pipeline
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  
  if (!use_raw) {
    cc <- projectData(cc, if (species == "mouse") PPI.mouse else PPI.human)
  }
  
  cc <- computeCommunProb(
    cc,
    raw.use = use_raw,
    type = "truncatedMean",
    trim = 0.1,
    population.size = FALSE
  )
  
  cc <- filterCommunication(cc, min.cells = min_cells)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)
  
  cc
}

# -----------------------------
# 4) Split by sex
# -----------------------------
female <- subset(combined_tumor, subset = sex == "Female")
male   <- subset(combined_tumor, subset = sex == "Male")

# Make sure sex subsets have clean, dropped group labels
female$celltype <- sanitize_group_factor(female$celltype)
male$celltype   <- sanitize_group_factor(male$celltype)

# Optional check
cat("Female celltype levels:", length(levels(female$celltype)), "\n")
cat("Male celltype levels:", length(levels(male$celltype)), "\n")

# -----------------------------
# 5) Run CellChat per sex
# -----------------------------
cellchat_female <- run_cellchat_v5(
  female,
  group_col = "celltype",
  assay = "RNA",
  species = "mouse",
  use_raw = FALSE,
  min_cells = 1
)

cellchat_male <- run_cellchat_v5(
  male,
  group_col = "celltype",
  assay = "RNA",
  species = "mouse",
  use_raw = FALSE,
  min_cells = 1
)

# -----------------------------
# 6) Merge & compare
# -----------------------------
cellchat_merged <- mergeCellChat(
  list(Female = cellchat_female, Male = cellchat_male),
  add.names = c("Female","Male")
)

compareInteractions(cellchat_merged, group = c(1,2), show.legend = FALSE)
netVisual_diffInteraction(cellchat_merged, measure = "count",  weight.scale = TRUE)
netVisual_diffInteraction(cellchat_merged, measure = "weight", weight.scale = TRUE)
rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)

# -----------------------------
# 7) Export CSV tables
# -----------------------------
outdir <- "CellChat_Female_vs_Male_CSV"
dir.create(outdir, showWarnings = FALSE)

df_female_all <- subsetCommunication(cellchat_female)
df_male_all   <- subsetCommunication(cellchat_male)
write.csv(df_female_all, file.path(outdir, "Female_all_LR_interactions.csv"), row.names = FALSE)
write.csv(df_male_all,   file.path(outdir, "Male_all_LR_interactions.csv"),   row.names = FALSE)

df_female_path <- subsetCommunication(cellchat_female, slot.name = "netP")
df_male_path   <- subsetCommunication(cellchat_male,   slot.name = "netP")
write.csv(df_female_path, file.path(outdir, "Female_pathway_interactions.csv"), row.names = FALSE)
write.csv(df_male_path,   file.path(outdir, "Male_pathway_interactions.csv"),   row.names = FALSE)

df_female_oscc_neur <- subsetCommunication(cellchat_female, sources.use = "OSCC", targets.use = "Neurons")
df_male_oscc_neur   <- subsetCommunication(cellchat_male,   sources.use = "OSCC", targets.use = "Neurons")
write.csv(df_female_oscc_neur, file.path(outdir, "Female_OSCC_to_Neurons_LR.csv"), row.names = FALSE)
write.csv(df_male_oscc_neur,   file.path(outdir, "Male_OSCC_to_Neurons_LR.csv"),   row.names = FALSE)

df_diff <- subsetCommunication(cellchat_merged, slot.name = "net", comparison = c(1,2))
write.csv(df_diff, file.path(outdir, "Female_vs_Male_differential_LR.csv"), row.names = FALSE)

rank_df <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE, return.data = TRUE)
write.csv(rank_df, file.path(outdir, "Female_vs_Male_pathway_ranking.csv"), row.names = FALSE)

write.csv(cellchat_female@net$centrality, file.path(outdir, "Female_signaling_centrality.csv"))
write.csv(cellchat_male@net$centrality,   file.path(outdir, "Male_signaling_centrality.csv"))

message("Done. CSVs written to: ", normalizePath(outdir))



# ============================================================
# Female vs Male OSCC <-> Neurons: delta LR, pathway comparison,
# bubble plots, heatmaps, and alluvial (incoming/outgoing) patterns
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(CellChat)
  library(ggalluvial)   # for netAnalysis_river
  library(NMF)          # for identifyCommunicationPatterns
})

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(patchwork)
  library(CellChat)
  library(NMF)
  library(ggalluvial)
})

# ============================================================
# A) Top OSCC -> Neurons barplots (Female + Male)
# ============================================================

top_n <- 20

# Make sure prob is numeric (very common CSV issue)
Female_OSCC_to_Neurons_LR$prob <- as.numeric(Female_OSCC_to_Neurons_LR$prob)
Male_OSCC_to_Neurons_LR$prob   <- as.numeric(Male_OSCC_to_Neurons_LR$prob)

df_female_top <- Female_OSCC_to_Neurons_LR %>%
  filter(source == "OSCC", target == "Neurons") %>%
  arrange(desc(prob)) %>%
  slice(1:top_n) %>%
  mutate(LR = str_trunc(paste0(ligand, " → ", receptor), 45)) %>%
  arrange(prob)

p_female <- ggplot(df_female_top, aes(x = prob, y = reorder(LR, prob))) +
  geom_col(fill = "#E64B35FF") +
  labs(
    title = "Female: Top OSCC → Neurons LR pairs",
    x = "CellChat communication probability",
    y = "Ligand → Receptor"
  ) +
  theme_classic(base_size = 12)

df_male_top <- Male_OSCC_to_Neurons_LR %>%
  filter(source == "OSCC", target == "Neurons") %>%
  arrange(desc(prob)) %>%
  slice(1:top_n) %>%
  mutate(LR = str_trunc(paste0(ligand, " → ", receptor), 45)) %>%
  arrange(prob)

p_male <- ggplot(df_male_top, aes(x = prob, y = reorder(LR, prob))) +
  geom_col(fill = "#4DBBD5FF") +
  labs(
    title = "Male: Top OSCC → Neurons LR pairs",
    x = "CellChat communication probability",
    y = "Ligand → Receptor"
  ) +
  theme_classic(base_size = 12)

# Save
ggsave("Female_OSCC_to_Neurons_topLR_barplot.pdf", p_female, width=11, height=8.5, units="in")
ggsave("Male_OSCC_to_Neurons_topLR_barplot.pdf",   p_male,   width=11, height=8.5, units="in")


# ============================================================
# B) Delta table + delta barplot (Female − Male) for same LR pairs
#    (No F/M objects needed)
# ============================================================

# Ensure pathway_name exists in both; if not, keep it anyway (join will handle)
df_delta <- full_join(
  Female_OSCC_to_Neurons_LR %>%
    filter(source=="OSCC", target=="Neurons") %>%
    transmute(source, target, ligand, receptor,
              pathway_name,
              prob_F = as.numeric(prob)),
  Male_OSCC_to_Neurons_LR %>%
    filter(source=="OSCC", target=="Neurons") %>%
    transmute(source, target, ligand, receptor,
              pathway_name,
              prob_M = as.numeric(prob)),
  by = c("source","target","ligand","receptor","pathway_name")
) %>%
  mutate(
    prob_F = ifelse(is.na(prob_F), 0, prob_F),
    prob_M = ifelse(is.na(prob_M), 0, prob_M),
    delta  = prob_F - prob_M,
    LR = str_trunc(paste0(ligand, " → ", receptor), 55),
    direction = ifelse(delta > 0, "Female higher", "Male higher")
  )

write.csv(df_delta, "Delta_OSCC_to_Neurons_all_LR.csv", row.names = FALSE)

# Plot top by absolute delta
top_delta <- 25
df_delta_top <- df_delta %>%
  arrange(desc(abs(delta))) %>%
  slice(1:top_delta) %>%
  arrange(delta)

p_delta <- ggplot(df_delta_top, aes(x = delta, y = reorder(LR, delta), fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(
    title = paste0("OSCC → Neurons: Top ", top_delta, " sex-differential LR pairs"),
    x = "Delta probability (Female − Male)",
    y = "Ligand → Receptor",
    fill = ""
  ) +
  theme_classic(base_size = 12)

print(p_delta)
ggsave("Delta_OSCC_to_Neurons_TopLR.pdf", p_delta, width=12, height=9, units="in")
ggsave("Delta_OSCC_to_Neurons_TopLR.png", p_delta, width=12, height=9, units="in", dpi=300)


# ============================================================
# C) Side-by-side dotplots: Female vs Male (Outgoing & Incoming)
#    Requires cellchat_female and cellchat_male objects
# ============================================================

# Ensure deterministic NMF
set.seed(123)
k_out <- 5
k_in  <- 5

# Identify patterns (outgoing + incoming) for each sex
cellchat_female <- identifyCommunicationPatterns(cellchat_female, pattern="outgoing", k=k_out)
cellchat_female <- identifyCommunicationPatterns(cellchat_female, pattern="incoming", k=k_in)

set.seed(123)
cellchat_male <- identifyCommunicationPatterns(cellchat_male, pattern="outgoing", k=k_out)
cellchat_male <- identifyCommunicationPatterns(cellchat_male, pattern="incoming", k=k_in)

# Create dotplots (these return ggplot objects in most CellChat versions)
pF_out <- netAnalysis_dot(cellchat_female, pattern="outgoing") + ggtitle("Female – Outgoing patterns")
pM_out <- netAnalysis_dot(cellchat_male,   pattern="outgoing") + ggtitle("Male – Outgoing patterns")

pF_in  <- netAnalysis_dot(cellchat_female, pattern="incoming") + ggtitle("Female – Incoming patterns")
pM_in  <- netAnalysis_dot(cellchat_male,   pattern="incoming") + ggtitle("Male – Incoming patterns")

# Side-by-side panels
p_out_side <- pF_out | pM_out
p_in_side  <- pF_in  | pM_in

# Save BIG
ggsave("DotPatterns_Outgoing_Female_vs_Male.pdf", p_out_side, width=18, height=7, units="in")
ggsave("DotPatterns_Incoming_Female_vs_Male.pdf", p_in_side,  width=18, height=7, units="in")

ggsave("DotPatterns_Outgoing_Female_vs_Male.png", p_out_side, width=18, height=7, units="in", dpi=300)
ggsave("DotPatterns_Incoming_Female_vs_Male.png", p_in_side,  width=18, height=7, units="in", dpi=300)


# ============================================================
# D) OPTIONAL: If netAnalysis_dot() does NOT return ggplot in your version
#    (rare) — save via PDF device instead
# ============================================================

# pdf("DotPatterns_Outgoing_Female.pdf", width=16, height=7)
# netAnalysis_dot(cellchat_female, pattern="outgoing")
# dev.off()
#
# pdf("DotPatterns_Outgoing_Male.pdf", width=16, height=7)
# netAnalysis_dot(cellchat_male, pattern="outgoing")
# dev.off()
#
# (Repeat for incoming)



suppressPackageStartupMessages({
  library(CellChat)
  library(dplyr)
  library(ggplot2)
  library(ggalluvial)
})

dir.create("CellChat_OSCC_Neurons", showWarnings = FALSE)

# -----------------------------
# Helper: pull ALL interactions (disable pval filtering)
# -----------------------------
get_comm_all <- function(cc, pval_thresh = 1) {
  subsetCommunication(cc, slot.name = "net", thresh = pval_thresh)
}

# -----------------------------
# Helper: restrict to OSCC <-> Neurons
# mode = "outgoing": OSCC -> Neurons only
# mode = "incoming": Neurons <- OSCC (same edges, just receiver framing)
# -----------------------------
restrict_oscc_neurons <- function(df, mode = c("outgoing","incoming")) {
  mode <- match.arg(mode)
  if (mode == "outgoing") {
    df %>% filter(source == "OSCC", target == "Neurons")
  } else {
    df %>% filter(target == "Neurons", source == "OSCC")
  }
}

# ============================================================
# 1) BUBBLE PLOTS: OSCC -> Neurons (Female + Male)
# ============================================================
pdf("CellChat_OSCC_Neurons/Bubble_Female_OSCC_to_Neurons.pdf", width=16, height=10)
netVisual_bubble(cellchat_female, sources.use = "OSCC", targets.use = "Neurons", remove.isolate = TRUE)
dev.off()

pdf("CellChat_OSCC_Neurons/Bubble_Male_OSCC_to_Neurons.pdf", width=16, height=10)
netVisual_bubble(cellchat_male, sources.use = "OSCC", targets.use = "Neurons", remove.isolate = TRUE)
dev.off()

# ============================================================
# 2) HEATMAPS: OSCC -> Neurons (Pathway x LR summary)
#    We'll build a clean heatmap-like tile plot from the table:
#    y = ligand→receptor, x = pathway, fill = prob
# ============================================================
make_lr_heatmap <- function(df, title, outfile, top_lr = 40) {
  df2 <- df %>%
    mutate(
      prob = as.numeric(prob),
      LR = paste0(ligand, "→", receptor)
    ) %>%
    arrange(desc(prob)) %>%
    slice(1:top_lr) %>%
    mutate(
      LR = factor(LR, levels = rev(unique(LR))),
      pathway_name = factor(pathway_name, levels = unique(pathway_name))
    )
  
  p <- ggplot(df2, aes(x = pathway_name, y = LR, fill = prob)) +
    geom_tile() +
    labs(title = title, x = "Pathway", y = "Top LR pairs", fill = "prob") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(outfile, p, width = 12, height = 10, units = "in")
  p
}

dfF_all <- get_comm_all(cellchat_female, pval_thresh = 1)
dfM_all <- get_comm_all(cellchat_male,   pval_thresh = 1)

dfF_osccN <- restrict_oscc_neurons(dfF_all, "outgoing")
dfM_osccN <- restrict_oscc_neurons(dfM_all, "outgoing")

make_lr_heatmap(dfF_osccN,
                "Female: OSCC → Neurons (top LR, grouped by pathway)",
                "CellChat_OSCC_Neurons/Heatmap_Female_OSCC_to_Neurons_topLR.pdf",
                top_lr = 50)

make_lr_heatmap(dfM_osccN,
                "Male: OSCC → Neurons (top LR, grouped by pathway)",
                "CellChat_OSCC_Neurons/Heatmap_Male_OSCC_to_Neurons_topLR.pdf",
                top_lr = 50)

# ============================================================
# 3) ALLUVIAL: OSCC -> Neurons flows by pathway
#    Outgoing: OSCC -> pathway -> Neurons
#    Incoming: Neurons -> pathway -> OSCC (same info, receiver framing)
# ============================================================
make_alluvial_oscc_neurons <- function(df, prefix, mode = c("outgoing","incoming"), top_pathways = 20) {
  mode <- match.arg(mode)
  
  df2 <- df %>%
    mutate(prob = as.numeric(prob)) %>%
    filter(source == "OSCC", target == "Neurons") %>%
    group_by(pathway_name, ligand, receptor) %>%
    summarise(flow = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    filter(flow > 0)
  
  # limit to top pathways for readability
  top_pw <- df2 %>%
    group_by(pathway_name) %>%
    summarise(total = sum(flow), .groups = "drop") %>%
    arrange(desc(total)) %>%
    slice(1:top_pathways) %>%
    pull(pathway_name)
  
  df2 <- df2 %>% filter(pathway_name %in% top_pw) %>%
    mutate(LR = paste0(ligand, "→", receptor))
  
  # Build a 3-axis alluvial
  # Outgoing framing: OSCC -> Pathway -> Neurons
  # Incoming framing:  Neurons -> Pathway -> OSCC
  if (mode == "outgoing") {
    df_plot <- df2 %>% transmute(A="OSCC", B=pathway_name, C="Neurons", flow=flow)
    ttl <- paste0(prefix, " OSCC → Neurons alluvial (Outgoing): OSCC → pathway → Neurons")
  } else {
    df_plot <- df2 %>% transmute(A="Neurons", B=pathway_name, C="OSCC", flow=flow)
    ttl <- paste0(prefix, " OSCC → Neurons alluvial (Incoming): Neurons → pathway → OSCC")
  }
  
  p <- ggplot(df_plot, aes(axis1=A, axis2=B, axis3=C, y=flow)) +
    geom_alluvium(aes(fill=B), alpha=0.85) +
    geom_stratum(width = 1/8) +
    geom_text(stat="stratum", aes(label = after_stat(stratum)), size=4) +
    scale_x_discrete(limits=c("Cell", "Pathway", "Cell"), expand=c(.05,.05)) +
    labs(title = ttl, y = "Sum(prob)") +
    theme_minimal(base_size = 12)
  
  ggsave(paste0("CellChat_OSCC_Neurons/Alluvial_", prefix, "_", toupper(mode), "_OSCC_Neurons.pdf"),
         p, width=14, height=6, units="in")
  
  p
}

make_alluvial_oscc_neurons(dfF_all, "Female", mode="outgoing", top_pathways=20)
make_alluvial_oscc_neurons(dfF_all, "Female", mode="incoming", top_pathways=20)

make_alluvial_oscc_neurons(dfM_all, "Male", mode="outgoing", top_pathways=20)
make_alluvial_oscc_neurons(dfM_all, "Male", mode="incoming", top_pathways=20)

message("DONE. Outputs are in: CellChat_OSCC_Neurons/")
