###cellchat
# Packages
library(Seurat)
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(ggplot2)
library(data.table)

# Choose species DB
CellChatDB <- CellChatDB.mouse # or CellChatDB.human
cellchat@DB <- CellChatDB
# Name your list for clarity (edit to match your order)
names(cellchat_sex_condition_list) <- c("Normal_Male", "Tumor_Male", "Normal_Female", "Tumor_Female")

run_cellchat <- function(cellchat, species = c("mouse","human"),
                         categories = "all",
                         nboot = 100, min.cells = 10, thresh = 0.1,
                         use.mouse.PPI = TRUE) {
  species <- match.arg(species)
  cellchat@DB <- if (species == "mouse") CellChatDB.mouse else CellChatDB.human
  
  # 1) subset DB
  cellchat <- subsetDB(cellchat, search = categories)
  
  # 2) overexpression tests
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # 3) project to PPI
  if (species == "mouse") {
    cellchat <- projectData(cellchat, PPI.mouse)
  } else {
    cellchat <- projectData(cellchat, PPI.human)
  }
  
  # 4) communication probability + filtering
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, nboot = nboot, population.size = TRUE, type = "triMean", thresh = thresh)
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)
  
  # 5) pathway-level probability + aggregation + centrality
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  return(cellchat)
}


processed_list <- lapply(cellchat_sex_condition_list, function(obj) {
  run_cellchat(obj, species = "mouse", categories = "all", nboot = 100, min.cells = 10, thresh = 0.1)
})

# Assign names and save each
names(processed_list) <- names(cellchat_sex_condition_list)
for (nm in names(processed_list)) {
  saveRDS(processed_list[[nm]], paste0("cellchat_", nm, "_processed.rds"))
}

# Merge
cellchat.merge <- mergeCellChat(processed_list, add.names = names(processed_list))
saveRDS(cellchat.merge, "cellchat_merge_4conditions.rds")

# Global comparison of interaction strength
p_strength <- compareInteractions(cellchat.merge, measure = "strength")
print(p_strength)

# Rank pathways per condition
rankNet(cellchat.merge, mode = "comparison", measure = "strength")

library(dplyr)

make_top20_barplot <- function(obj, title_suffix) {
  path_scores <- data.frame(
    pathway = names(obj@netP$prob),
    total_prob = sapply(obj@netP$prob, function(x) sum(x, na.rm = TRUE))
  ) %>% arrange(desc(total_prob)) %>% head(20)
  
  ggplot(path_scores, aes(x = reorder(pathway, total_prob), y = total_prob)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste("Top 20 Pathways -", title_suffix),
         x = "Pathway", y = "Total Interaction Strength (prob)") +
    theme_minimal(base_size = 14)
}

plots_bar <- lapply(names(processed_list), function(nm) make_top20_barplot(processed_list[[nm]], nm))
# Example: print first
print(plots_bar[[1]])


netVisual_bubble(cellchat.merge, signaling = c("VEGF","TGFb","COLLAGEN","CXCL", "TNF", "LIGHT", "LTA", "LTB", "TRAF2", "NIK", "	NF-kappa B"),
                 remove.isolate = TRUE, title.name = "Selected pathways across conditions")

for (nm in names(processed_list)) {
  obj <- processed_list[[nm]]
  netVisual_circle(obj@net$count,
                   vertex.weight = as.numeric(table(obj@idents)),
                   weight.scale = TRUE, label.edge = FALSE, alpha.edge = 0.6,
                   title.name = paste("Global network -", nm))
}


netVisual_heatmap(cellchat.merge, measure = "pathway", color.heatmap = "RdBu")


# Export ligand–receptor communication tables per condition
for (nm in names(processed_list)) {
  df.net <- subsetCommunication(processed_list[[nm]])
  fwrite(df.net, paste0("CellChat_interactions_all_", nm, ".csv"))
  df.net.sig <- subsetCommunication(processed_list[[nm]], pval = 0.05)
  fwrite(df.net.sig, paste0("CellChat_interactions_sig_", nm, ".csv"))
}

# Export pathway matrices per condition
for (nm in names(processed_list)) {
  obj <- processed_list[[nm]]
  for (p in names(obj@netP$prob)) {
    mat <- obj@netP$prob[[p]]  # rows = sources, cols = targets
    fwrite(data.table(Source = rownames(mat), mat),
           paste0("CellChat_pathway_prob_", p, "_", nm, ".csv"))
  }
}


# Helper to compute Tumor - Normal per sex
diff_top20 <- function(obj_tumor, obj_normal, title) {
  paths <- intersect(names(obj_tumor@netP$prob), names(obj_normal@netP$prob))
  df <- data.frame(
    pathway = paths,
    diff = sapply(paths, function(p) sum(obj_tumor@netP$prob[[p]], na.rm = TRUE) -
                    sum(obj_normal@netP$prob[[p]], na.rm = TRUE))
  ) %>% arrange(desc(diff)) %>% head(20)
  
  ggplot(df, aes(x = reorder(pathway, diff), y = diff)) +
    geom_bar(stat = "identity", fill = ifelse(df$diff > 0, "firebrick", "steelblue")) +
    coord_flip() +
    labs(title = paste("Tumor – Normal (", title, ")", sep=""),
         x = "Pathway", y = "Δ Total Interaction Strength") +
    theme_minimal(base_size = 14)
}

p_male_diff   <- diff_top20(processed_list[["Tumor_Male"]],   processed_list[["Normal_Male"]],   "Male")
p_female_diff <- diff_top20(processed_list[["Tumor_Female"]], processed_list[["Normal_Female"]], "Female")
print(p_male_diff); print(p_female_diff)

# Suppose your processed CellChat object is called cellchat
saveRDS(cellchat, file = "cellchat_results.rds")