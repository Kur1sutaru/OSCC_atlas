# install.packages( "ggthemes" )
library( ggthemes )
source( "../SC_CommonFunctions.R")
source("../Shivani_Marmoset/marmoset_mouse_tongue_markers.r")

SC_GetSampleRDS <- function( fileName, sampleName ) {
  scData.norm.scaled  = readRDS( file = fileName )

  numClusters <- nlevels( scData.norm.scaled )
  new.cluster.ids <- paste0( sampleName, "_", c(0:(numClusters-1)) ) 
  names(new.cluster.ids) <- levels(scData.norm.scaled)
  scData.norm.scaled <- RenameIdents(scData.norm.scaled, new.cluster.ids)
}

reference.list <- list()
sampleNames = c( "FC", "FC2", "FC3", "MC", "MC2", "MC3" )
for (i in 1:length( sampleNames )) {
  reference.list[[i]] <- SC_GetSampleRDS( fileName = paste0( tolower(sampleNames[i]), '.norm.scaled.singlets.rds'),
                                          sampleName = sampleNames[i] )
}

# find anchors
anchors <- FindIntegrationAnchors( object.list = reference.list, dims = 1:50 )

# normal control merged.
OSCC.nc.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
DefaultAssay(OSCC.nc.integrated) <- "integrated"

#saveRDS(OSCC.nc.integrated, file = "OSCC.nc.integrated.rds")
#OSCC.nc.integrated = readRDS( "OSCC.nc.integrated.rds" )


# Run the standard workflow for visualization and clustering
OSCC.nc.integrated <- ScaleData(OSCC.nc.integrated, verbose = FALSE)
OSCC.nc.integrated <- RunPCA(OSCC.nc.integrated, npcs = 50, verbose = FALSE)
OSCC.nc.integrated <- RunUMAP(OSCC.nc.integrated, reduction = "pca", dims = 1:50, verbose = FALSE)

p1 <- DimPlot(OSCC.nc.integrated, reduction = "umap", group.by = "orig.ident")
ggsave("OSCC_NC_integrated_Umap_orig_ident.PDF", p1, width = 10, height = 10)
p2 <- DimPlot(OSCC.nc.integrated, reduction = "umap", label = TRUE, repel = TRUE)
ggsave("OSCC_NC_integrated_Umap.PDF", p2, width = 40, height = 20)
ggsave("OSCC_NC_integrated_Umap_combined.PDF", p1+p2, width = 40, height = 10)

saveRDS(OSCC.nc.integrated, file = "OSCC.nc.integrated.rds")
#OSCC.nc.integrated = readRDS( "OSCC.nc.integrated.rds" )

OSCC.nc.integrated <- RunTSNE(OSCC.nc.integrated, reduction = "pca", dims = 1:50, perplexity = 30)
p2 <- DimPlot(OSCC.nc.integrated, reduction = "tsne", label = TRUE, repel = TRUE)
ggsave("OSCC_NC_integrated_TSNE.PDF", p2, width = 40, height = 20)
OSCC.nc.integrated <- AddMetaData( object = OSCC.nc.integrated, metadata = as.character( OSCC.nc.integrated@active.ident ), col.name = "clusterIDs" )

AssignNewNames = function( cellType_toUpdate, sampleName, cellPrefix ) {
  sample.celltypes = read.xlsx( paste0( sampleName, "_custom_Mm-celltypes.xlsx" ) )
  for( i in 1:nrow(sample.celltypes) ) {
    idx = which( cellType_toUpdate == paste0( sampleName, "_", i-1 ) )
    cellType_toUpdate[idx] = paste0( cellPrefix, "_", sample.celltypes$customTypes[i] )
#    cellType_toUpdate[idx] = paste0( cellPrefix, "_", sample.celltypes$tumorCells[i] )
  }
  return( cellType_toUpdate )
}

#OSCC.nc.integrated = readRDS( "OSCC.nc.integrated.rds" )

#cellTypes = OSCC.nc.integrated@active.ident
cellTypes = OSCC.nc.integrated@meta.data$clusterIDs
newCellTypes = as.character( cellTypes )

newCellTypes = AssignNewNames( newCellTypes, "FC", "fc" )
newCellTypes = AssignNewNames( newCellTypes, "MC", "mc" )
newCellTypes = AssignNewNames( newCellTypes, "MC2", "mc2" )
newCellTypes = AssignNewNames( newCellTypes, "MC3", "mc3" )
newCellTypes = AssignNewNames( newCellTypes, "FC2", "fc2" )
newCellTypes = AssignNewNames( newCellTypes, "FC3", "fc3" )

OSCC.nc.integrated@active.ident = as.factor( newCellTypes )
p2 <- DimPlot(OSCC.nc.integrated, reduction = "umap", pt.size = 2, label = TRUE, repel = TRUE )
ggsave("OSCC_NC_integrated_CellTypes_UMAP.PDF", p2, width = 40, height = 20, limitsize = FALSE )
p2 <- DimPlot(OSCC.nc.integrated, reduction = "tsne", pt.size = 2, label = TRUE, repel = TRUE )
ggsave("OSCC_NC_integrated_CellTypes_TSNE.PDF", p2, width = 40, height = 20, limitsize = FALSE )

DefaultAssay( OSCC.nc.integrated ) = "RNA"
g2 = FeaturePlot( OSCC.nc.integrated, max.cutoff = 'q90', min.cutoff = 'q10', features = c("Pecam1", "Vwf", "Plvap", "Il6st", "Pdgfb", "Cldn5",  "Cdh5", "Rflnb"), ncol = 3, split.by = "orig.ident" )
g2 = FeaturePlot( OSCC.nc.integrated, max.cutoff = 'q90', min.cutoff = 'q10', features = c("Rflnb", "Fabp4", "Sptbn1", "Gngt2"), ncol = 3, split.by = "orig.ident", reduction = "umap" )
g2
ggsave( "OSCC_NC_integrated_FeatureMap_Endo_markerGenes.PDF", g2, width = 15, height = 18, limitsize = FALSE )
ggsave( "OSCC_NC_integrated_FeatureMap_Endo_markerGenes.png", g2, width = 15, height = 18, limitsize = FALSE )

# prepare for dotplot.
tmp = strsplit( unique( newCellTypes ), "_" )
all_ct = ""
for( i in 1:length(tmp)) { all_ct[i] = tmp[[i]][2] }
all_ct = sort( unique( all_ct ) )

uCellTypes = unique( newCellTypes )
tmp = strsplit( uCellTypes, "_" )
all_ct = ""
for( i in 1:length(tmp)) { all_ct[i] = substr(uCellTypes[i], nchar(tmp[[i]][1])+2, stop = 1000 ) }
all_ct = sort( unique( all_ct ) )

sample_ct = unique( newCellTypes )
for( j in 1:length( all_ct ) ) {
  ct_string = ""
  for( i in 1:length(sampleNames) ) {
    ct = paste0( tolower( sampleNames[i] ), "_", all_ct[j])
    if( ct %in% sample_ct ) {
      if( nchar( ct_string) > 0  ) {
        ct_string = paste0( ct_string, ",'", ct, "'" )
      } else {
        ct_string = paste0( "'", ct, "'" )
      }
    }
  }
  cat( ct_string )
  cat( "\n" )
}

OSCC.nc.integrated <- AddMetaData( object = OSCC.nc.integrated, metadata = newCellTypes, col.name = 'celltype' )
# plot dotplot with all 8 samples.
OSCC.nc.integrated$celltype <- factor(OSCC.nc.integrated$celltype, levels = c(
  'fc_Endothelial','fc2_Endothelial','fc3_Endothelial','mc_Endothelial','mc2_Endothelial','mc3_Endothelial',
  'fc_Endothelial2','fc2_Endothelial2','fc3_Endothelial2','mc_Endothelial2','mc2_Endothelial2','mc3_Endothelial2',
  'fc_Angiogenic ECs','fc2_Angiogenic ECs','fc3_Angiogenic ECs','mc_Angiogenic ECs','mc2_Angiogenic ECs','mc3_Angiogenic ECs',
  'fc_Arterial_Contractile ECs','fc2_Arterial_Contractile ECs','fc3_Arterial_Contractile ECs','mc_Arterial_Contractile ECs','mc2_Arterial_Contractile ECs','mc3_Arterial_Contractile ECs',
  'fc_Interferon_responsive ECs','fc3_Interferon_responsive ECs','mc_Interferon_responsive ECs','mc2_Interferon_responsive ECs','mc3_Interferon_responsive ECs',
  'fc_Lymphatic endothelial','fc2_Lymphatic endothelial','fc3_Lymphatic endothelial','mc_Lymphatic endothelial','mc2_Lymphatic endothelial','mc3_Lymphatic endothelial',
  'fc2_Pericytes','fc3_Pericytes','mc_Pericytes','mc2_Pericytes','mc3_Pericytes',
  'fc_Smooth muscle cell','fc2_Smooth muscle cell','fc3_Smooth muscle cell','mc_Smooth muscle cell','mc2_Smooth muscle cell','mc3_Smooth muscle cell',
  'fc_Skeletal muscle','fc2_Skeletal muscle','fc3_Skeletal muscle','mc2_Skeletal muscle','mc3_Skeletal muscle',
  'fc2_Muscle satellite cells','fc3_Muscle satellite cells','mc2_Muscle satellite cells','mc3_Muscle satellite cells',
  'fc_Fibroblasts','fc2_Fibroblasts','fc3_Fibroblasts','mc_Fibroblasts','mc2_Fibroblasts','mc3_Fibroblasts',
  'fc2_Basal keratinocyte','fc3_Basal keratinocyte','mc_Basal keratinocyte','mc2_Basal keratinocyte','mc3_Basal keratinocyte',
  'fc_Differentiating keratinocyte','fc2_Differentiating keratinocyte', 'fc3_Differentiating keratinocyte','mc_Differentiating keratinocyte','mc3_Differentiating keratinocyte',
  'fc2_T cell','mc3_T cell',
  'fc_Macrophages','fc2_Macrophages','fc3_Macrophages','mc_Macrophages','mc2_Macrophages','mc3_Macrophages',
  'fc_Myelinating Schwann cells','fc2_Myelinating Schwann cells','fc3_Myelinating Schwann cells','mc_Myelinating Schwann cells','mc2_Myelinating Schwann cells','mc3_Myelinating Schwann cells',
  'fc_Unmyelinating Schwann cells','fc3_Unmyelinating Schwann cells','mc_Unmyelinating Schwann cells',
  'fc2_Mast','mc_Mast',
  'fc2_Blood cell','mc2_Blood cell'))


DefaultAssay(OSCC.nc.integrated) <- "RNA"
Idents(OSCC.nc.integrated) <- "celltype"
marker_short_list_clean = SC_CleanMarkerList( OSCC.nc.integrated, marker_short_prog_list )
DotPlot( object = OSCC.nc.integrated, features = marker_short_list_clean, group.by = "celltype",
         cluster.idents=F, col.min = 0, scale.min = 0 ) + 
  theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) )
ggsave( "OSCC_NC_integrated_BeforeCleanUp_DotPlot.pdf", width = 30, height = 15 )

# move clustering code here. 
saveRDS(OSCC.nc.integrated, file = "OSCC.nc.integrated.rds")

# fix cluseters (according to integrated groups)
DefaultAssay(OSCC.nc.integrated ) <- "integrated"
OSCC.nc.integrated <- FindNeighbors(OSCC.nc.integrated, reduction = "pca", dims = 1:50)
OSCC.nc.integrated <- FindClusters(OSCC.nc.integrated, resolution = 2 )

OSCC.nc.integrated <- AddMetaData( object = OSCC.nc.integrated,  
                                   metadata = as.character( OSCC.nc.integrated@active.ident ), 
                                   col.name = "clusterIDs_integrated" )

p2 <- DimPlot( OSCC.nc.integrated, reduction = "umap", label = TRUE, repel = TRUE )
ggsave("OSCC_NC_integrated_ClusterAgain_Umap.PDF", p2, width = 10, height = 8)
p2 <- DimPlot( OSCC.nc.integrated, reduction = "tsne", label = TRUE, repel = TRUE )
ggsave("OSCC_NC_integrated_ClusterAgain_TSNE.PDF", p2, width = 10, height = 8)

Idents(OSCC.nc.integrated) <- "clusterIDs_integrated"
newIDs = read.delim( "Integrated_cellTypes-v2.txt", header = TRUE, stringsAsFactors = FALSE )
rownames( newIDs ) = newIDs$cluster
new.cluster.ids <- newIDs[as.character(levels( OSCC.nc.integrated )),2]
names(new.cluster.ids) <- levels( OSCC.nc.integrated )
OSCC.nc.integrated <- RenameIdents( OSCC.nc.integrated, new.cluster.ids )
p2 <- DimPlot(OSCC.nc.integrated, reduction = "umap", pt.size = 2, label = TRUE, repel = TRUE )
ggsave("OSCC_NC_integrated_merged_UMAP.PDF", p2, width = 11, height = 8, limitsize = FALSE )
p2 <- DimPlot(OSCC.nc.integrated, reduction = "tsne", pt.size = 2, label = TRUE, repel = TRUE )
ggsave("OSCC_NC_integrated_merged_TSNE.PDF", p2, width = 11, height = 8, limitsize = FALSE )

# Level 1 cell types
Idents(OSCC.nc.integrated) <- "clusterIDs_integrated"
newIDs = read.delim( "Integrated_cellTypes-level1.txt", header = TRUE, stringsAsFactors = FALSE )
rownames( newIDs ) = newIDs$cluster
new.cluster.ids <- newIDs[as.character(levels( OSCC.nc.integrated )),2]
names(new.cluster.ids) <- levels( OSCC.nc.integrated )
OSCC.nc.integrated <- RenameIdents( OSCC.nc.integrated, new.cluster.ids )
p2 <- DimPlot(OSCC.nc.integrated, reduction = "umap", pt.size = 2, label = TRUE, repel = TRUE )
ggsave("OSCC_NC_integrated_merged_level1_UMAP.PDF", p2, width = 11, height = 8, limitsize = FALSE )
p2 <- DimPlot(OSCC.nc.integrated, reduction = "tsne", pt.size = 2, label = TRUE, repel = TRUE )
ggsave("OSCC_NC_integrated_merged_level1_TSNE.PDF", p2, width = 11, height = 8, limitsize = FALSE )


DefaultAssay( OSCC.nc.integrated ) = "RNA"
g2 = FeaturePlot( OSCC.nc.integrated, features = c("Mpz"), max.cutoff = 'q90' )
ggsave( "OSCC_NC_integrated_FeatureMap_Mpz_schwannCells.PDF", g2, width = 7, height = 8 )

# remove blood cells.
bloodCell_ids = 49
integrated_clusterIDs = OSCC.nc.integrated@meta.data$clusterIDs_integrated
sample_clusterIDs = OSCC.nc.integrated@meta.data$clusterIDs
ids_to_remove = which( sample_clusterIDs == "FC2_26" | sample_clusterIDs == "MC2_25" | integrated_clusterIDs == 49 )

OSCC.nc.integrated <- subset(OSCC.nc.integrated, cells = ids_to_remove, invert = TRUE )


# fix fc2_Late differentiated non−keratinized or partially keratinized keratinocyte
#endo_ids = c(0, 1, 2, 4, 6, 30, 7, 10, 12, 14, 17, 23, 25, 31, 38 )
#angioEC_ids = 28
#bkcyte_ids = c( 18, 37, 40)
#integrated_clusterIDs = OSCC.nc.integrated@meta.data$clusterIDs_integrated
#sampleIDs = OSCC.nc.integrated@meta.data$orig.ident
#ids = which( sampleIDs == "FC2" & integrated_clusterIDs %in% endo_ids )
#newCellTypes[ids] = "fc2_Endothelial"
#ids = which( sampleIDs == "FC2" & integrated_clusterIDs %in% bkcyte_ids )
#newCellTypes[ids] = "fc2_Basal keratinocyte"
#ids = which( sampleIDs == "FC2" & integrated_clusterIDs %in% angioEC_ids )
#newCellTypes[ids] = "fc2_Angiogenic ECs"

DefaultAssay( OSCC.nc.integrated ) = "RNA"
OSCC.nc.integrated <- AddMetaData( object = OSCC.nc.integrated, metadata = newCellTypes, col.name = 'celltype' )
# plot dotplot with all 8 samples.
OSCC.nc.integrated$celltype <- factor(OSCC.nc.integrated$celltype, levels = c(
  'fc_Endothelial','fc2_Endothelial','fc3_Endothelial','mc_Endothelial','mc2_Endothelial','mc3_Endothelial',
  'fc_Endothelial2','fc2_Endothelial2','fc3_Endothelial2','mc_Endothelial2','mc2_Endothelial2','mc3_Endothelial2',
  'fc_Angiogenic ECs','fc2_Angiogenic ECs','fc3_Angiogenic ECs','mc_Angiogenic ECs','mc2_Angiogenic ECs','mc3_Angiogenic ECs',
  'fc_Arterial_Contractile ECs','fc2_Arterial_Contractile ECs','fc3_Arterial_Contractile ECs','mc_Arterial_Contractile ECs','mc2_Arterial_Contractile ECs','mc3_Arterial_Contractile ECs',
  'fc_Interferon_responsive ECs','fc3_Interferon_responsive ECs','mc_Interferon_responsive ECs','mc2_Interferon_responsive ECs','mc3_Interferon_responsive ECs',
  'fc_Lymphatic endothelial','fc2_Lymphatic endothelial','fc3_Lymphatic endothelial','mc_Lymphatic endothelial','mc2_Lymphatic endothelial','mc3_Lymphatic endothelial',
  'fc2_Pericytes','fc3_Pericytes','mc_Pericytes','mc2_Pericytes','mc3_Pericytes',
  'fc_Smooth muscle cell','fc2_Smooth muscle cell','fc3_Smooth muscle cell','mc_Smooth muscle cell','mc2_Smooth muscle cell','mc3_Smooth muscle cell',
  'fc_Skeletal muscle','fc2_Skeletal muscle','fc3_Skeletal muscle','mc2_Skeletal muscle','mc3_Skeletal muscle',
  'fc2_Muscle satellite cells','fc3_Muscle satellite cells','mc2_Muscle satellite cells','mc3_Muscle satellite cells',
  'fc_Fibroblasts','fc2_Fibroblasts','fc3_Fibroblasts','mc_Fibroblasts','mc2_Fibroblasts','mc3_Fibroblasts',
  'fc2_Basal keratinocyte','fc3_Basal keratinocyte','mc_Basal keratinocyte','mc2_Basal keratinocyte','mc3_Basal keratinocyte',
  'fc_Differentiating keratinocyte','fc2_Differentiating keratinocyte', 'fc3_Differentiating keratinocyte','mc_Differentiating keratinocyte','mc3_Differentiating keratinocyte',
  'fc2_T cell','mc3_T cell',
  'fc_Macrophages','fc2_Macrophages','fc3_Macrophages','mc_Macrophages','mc2_Macrophages','mc3_Macrophages',
  'fc_Myelinating Schwann cells','fc2_Myelinating Schwann cells','fc3_Myelinating Schwann cells','mc_Myelinating Schwann cells','mc2_Myelinating Schwann cells','mc3_Myelinating Schwann cells',
  'fc_Unmyelinating Schwann cells','fc3_Unmyelinating Schwann cells','mc_Unmyelinating Schwann cells',
  'fc2_Mast','mc_Mast'))

DefaultAssay(OSCC.nc.integrated) <- "RNA"
Idents(OSCC.nc.integrated) <- "celltype"
marker_short_list_clean = SC_CleanMarkerList( OSCC.nc.integrated, marker_short_prog_list )
DotPlot( object = OSCC.nc.integrated, features = marker_short_list_clean, group.by = "celltype",
         cluster.idents=F, col.min = 0, scale.min = 0 ) + 
  theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) )
ggsave( "OSCC_NC_integrated_AfterCleanUp_DotPlot.pdf", width = 30, height = 18 )

saveRDS(OSCC.nc.integrated, file = "OSCC.nc.integrated_clean.rds")

# collapse all tissues. 
cellTypeOnly = newCellTypes
tmp = strsplit( cellTypeOnly, "_" )
for( i in 1:length(tmp)) { cellTypeOnly[i] = substr(cellTypeOnly[i], nchar(tmp[[i]][1])+2, stop = 1000 ) }

OSCC.nc.integrated <- AddMetaData( object = OSCC.nc.integrated, metadata = cellTypeOnly, col.name = 'celltype_merged' )
# plot dotplot with all 8 samples.
OSCC.nc.integrated$celltype_merged <- factor(OSCC.nc.integrated$celltype_merged, levels = c(
  "Endothelial",
  "Endothelial2",
  "Angiogenic ECs",
  "Arterial_Contractile ECs",
  "Interferon_responsive ECs",
  "Lymphatic endothelial",
  "Pericytes",
  "Smooth muscle cell",
  "Skeletal muscle",
  "Muscle satellite cells",
  "Fibroblasts",
  "Basal keratinocyte",
  "Differentiating keratinocyte",
  "T cell",
  "Macrophages",
  "Myelinating Schwann cells",
  "Unmyelinating Schwann cells",
  "Mast" ))

DefaultAssay(OSCC.nc.integrated) <- "RNA"
Idents(OSCC.nc.integrated) <- "celltype_merged"
marker_short_list_clean = SC_CleanMarkerList( OSCC.nc.integrated, marker_short_prog_list )
DotPlot( object = OSCC.nc.integrated, features = marker_short_list_clean, group.by = "celltype_merged",
         cluster.idents=F, col.min = 0, scale.min = 0 ) + 
  theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) )
ggsave( "OSCC_NC_integrated_merged_prog_DotPlot.pdf", width = 26, height = 6 )


marker_endo_list <-
  list( "Endothelial"               = c( "Pecam1", "Prom1", "Cd34", "Fabp4", "Tie1", "Itga6", "Cdh5", "Cldn5", "Vwf", "Plvap" ),
        "Venus_inflamed ECs"        = c( "Ackr1", "Vcam1", "Sele", "Selp", "Icam1" ),
        "Angiogenic ECs"            = c( "Esm1", "Tnfaip2", "Penk", "Kdr", "Igfbp7", "Eng" ),
        "Arterial_Contractile ECs"  = c( "Rflnb", "Ctnnb1", "Ctnnbip1", "Tle5", "Epas1", "Gpihbp1", "Cox7a2l", "Timp4" ),
        "Interferon-Responsive ECs" = c( "Ifit1", "Isg15", "Irf7", "Rsad2", "Usp18", "Stat1" ),
        "Lymphatic endothelial"     = c( "Prox1", "Lyve1", "Flt4" )
  )
DefaultAssay(OSCC.nc.integrated) <- "RNA"
Idents(OSCC.nc.integrated) <- "celltype"
marker_endo_list_clean = SC_CleanMarkerList( OSCC.nc.integrated, marker_endo_list )
gg = DotPlot( object = OSCC.nc.integrated, features = marker_endo_list_clean, group.by = "celltype",
         cluster.idents=F, col.min = 0, scale.min = 0 ) + 
  theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) )
ggsave( "OSCC_NC_integrated_Endothelial_markers_DotPlot.pdf", gg, width = 20, height = 18 )


g2 = FeaturePlot( OSCC.nc.integrated, max.cutoff = 'q90', min.cutoff = 'q10', 
                  features = c( "Ntrk1", "Ntrk2", "Ngfr", "Ngf", "Nlrp3", "Trpm2", "Bdnf" ), label = TRUE )
ggsave( "OSCC_integrated_FeatureMap_selectedGenes.PDF", g2, width = 25, height = 15, limitsize = FALSE )
VlnPlot( OSCC.nc.integrated,
        features = c( "Ntrk1", "Ntrk2", "Ngfr", "Ngf", "Nlrp3", "Trpm2", "Bdnf" ) )
ggsave( "OSCC_integrated_ViolinPlot_selectedGenes.PDF", width = 25, height = 15, limitsize = FALSE )


saveRDS(OSCC.nc.integrated, file = "OSCC.nc.integrated_clean.rds")
#OSCC.nc.integrated = readRDS( "OSCC.nc.integrated_clean.rds" )


#
#
# plot barcharts, excel files, etc.
#
#OSCC.nc.integrated = readRDS( "OSCC.nc.integrated_clean.rds" )

scMatrixSave <- function( matrix, xlsName, htmlName ) {
  library( writexl )
  library( htmlwidgets )
  library( DT )
  
  write_xlsx( matrix, path = xlsName )
  
  dtable <- datatable(
    matrix, extensions = 'Buttons',
    options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
      pageLength = 20
    )
  )
  
  # Save as an HTML file
  saveWidget(dtable, htmlName )
}

library(dplyr)
df = data.frame( id = as.character( OSCC.nc.integrated$celltype_merged ), sample = OSCC.nc.integrated@meta.data$orig.ident)
cTable = dplyr::count( df, id, sample )

cTable_long = pivot_longer( cTable, 2)
colnames( cTable_long ) = c( "cell_types", "value", "name", "Samples")

library( RColorBrewer )
library(Polychrome)
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(cTable_long$cell_types)))
#mycolors = createPalette(19,  c("#ff0000", "#00ff00", "#0000ff"))

ggbarplot(cTable_long, "Samples", "value", title = "cell counts per sample",
          fill = "cell_types", color = "cell_types", palette = mycolors, 
          xlab = "Samples", ylab = "# of cells") + 
  theme(legend.position='right', legend.key.height = unit(0.3, "cm"))
ggsave( "Mouse_NC_6sample_cellcounts.pdf", width = 7, height = 5 )

cTable = dplyr::count( df, id, sample )
cMatrix = pivot_wider( cTable, names_from = sample, values_from = n)
cMatrix[is.na( cMatrix )] = 0
scMatrixSave( cMatrix, "Mouse_NC_6sample_cellCounts_byCellTypes.xlsx", "Mouse_NC_6sample_cellCounts.html" )

cells = apply( cMatrix[,2:7], 2, sum)
for( i in 1:length(cells) ) {
  cMatrix[, i+1] = cMatrix[, i+1]/cells[i]
}
cMatrix_out = cMatrix
cMatrix_out[, 2:7] = round( cMatrix_out[, 2:7], 4)
scMatrixSave( cMatrix_out, "Mouse_NC_6sample_cellProportion_byCellTypes.xlsx", "Mouse_NC_6sample_cellcount_proportion.html" )

cTable_long2 = pivot_longer( cMatrix, 2:7 )
colnames( cTable_long2 ) = c( "cell_types", "Samples", "value")

ggbarplot(cTable_long2, "Samples", "value", title = "percents of cell types per sample",
          fill = "cell_types", color = "cell_types", palette = mycolors, 
          xlab = "Samples", ylab = "Percent of cells") + 
  theme(legend.position='right', legend.key.height = unit(0.3, "cm"))
ggsave( "Mouse_NC_6sample_celltypeProportions.pdf", width = 9, height = 5 )

cMatrix$FvsM_pval = 1
for( i in 1:nrow(cMatrix) ) {
  s = t.test( cMatrix[i, 2:4], cMatrix[i, 5:7] )
  cMatrix$FvsM_pval[i] = s$p.value
}
cMatrix_out = cMatrix
cMatrix_out[, 2:8] = round( cMatrix_out[, 2:8], 4)
scMatrixSave( cMatrix_out, "Mouse_NC_6sample_cellProportion_byCellTypes.xlsx", "Mouse_NC_6sample_cellcount_proportion.html" )

cTable = dplyr::count( df, id, sample )
cMatrix = pivot_wider( cTable, names_from = sample, values_from = n)
cMatrix[is.na( cMatrix )] = 0

cMatrixFM = cMatrix[, 1:3]
cMatrixFM[, 2] = apply( cMatrix[,2:4], 1, sum)
cMatrixFM[, 3] = apply( cMatrix[,5:7], 1, sum)
colnames(cMatrixFM) = c( "id", "Female", "Male")

cTable_long3 = pivot_longer( cMatrixFM, 2:3)
colnames( cTable_long3 ) = c( "cell_types", "gender", "value")

mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(cTable_long3$cell_types)))

ggbarplot(cTable_long3, "gender", "value", title = "cell counts by gender",
          fill = "cell_types", color = "cell_types", palette = mycolors, 
          xlab = "Gender", ylab = "# of cells") + 
  theme(legend.position='right', legend.key.height = unit(0.3, "cm"))
ggsave( "Mouse_NC_cellcounts_byGender.pdf", width = 7, height = 5 )

cells = apply( cMatrixFM[,2:3], 2, sum)
for( i in 1:length(cells) ) {
  cMatrixFM[, i+1] = cMatrixFM[, i+1]/cells[i]
}
cTable_long3 = pivot_longer( cMatrixFM, 2:3 )
colnames( cTable_long3 ) = c( "cell_types", "gender", "value")
ggbarplot(cTable_long3, "gender", "value", title = "cell count proportions by gender",
          fill = "cell_types", color = "cell_types", palette = mycolors, 
          xlab = "Gender", ylab = "Proportion of cells") + 
  theme(legend.position='right', legend.key.height = unit(0.3, "cm"))
ggsave( "Mouse_NC_celltypeProportions_byGender.pdf", width = 7, height = 5 )

logM = data.frame( name = cMatrixFM[,1], log2( cMatrixFM[, 2:3] ))
ggscatter( as.data.frame( logM ), x = "Female", y = "Male", 
           add = "reg.line", conf.int = TRUE, conf.int.level = 0.95,
           label = "id", repel = TRUE, 
           xlab = "Female - cell count proportion, log2", 
           ylab = "Male - cell count proportion, log2" )
ggsave("Mouse_NC_female_male_cellcount_scatterPlot.pdf", width = 8, height = 6 )

# save stats only.
cTable = dplyr::count( df, id, sample )
cMatrix = pivot_wider( cTable, names_from = sample, values_from = n)
cMatrix[is.na( cMatrix )] = 0

cMatrixFM = cMatrix[, 1:5]
cMatrixFM[, 2] = apply( cMatrix[,2:4], 1, sum)
cMatrixFM[, 3] = apply( cMatrix[,5:7], 1, sum)

cells = apply( cMatrixFM[,2:3], 2, sum)
for( i in 1:length(cells) ) {
  cMatrixFM[, i+3] = cMatrixFM[, i+3]/cells[i]
}
colnames(cMatrixFM) = c( "id", "Female_count", "Male_count", "Female_proportion", "Male_proportion")
cMatrixFM[, 4:5] = round( cMatrixFM[, 4:5], 3)
scMatrixSave( cMatrixFM, "Mouse_NC_cellCount_proportion_byGender.xlsx", "Mouse_NC_cellCount_proportion_byGender.html" )


#
#
# differential expression analysis

# Get DEGs into a file (only male vs female).
#BiocManager::install( "EnhancedVolcano" )
library( EnhancedVolcano )
CS_GetDEGs <- function( scData, cellType, assay = "integrated" ) {
  # create matrix to store the results. (each DEG list contains 5 columns)
  # if assay is RNA, we have to handle differently.
  
  mData = scData@meta.data
  if( length( which( mData$gender_celltype == paste0( "male_", cellType ) )) > 0 &
      length( which( mData$gender_celltype == paste0( "female_", cellType ) )) > 0 ) {
    
    if( assay == "RNA" ) {
      degs = FindMarkers( scData, ident.1 = paste0( "male_", cellType ), 
                                  ident.2 = paste0( "female_", cellType ), 
                          group.by = "gender_celltype",  min.prt = 0.1, logfc.threshold = 0.1, 
                          assay = "RNA" )
    } else { 
      degs = FindMarkers( scData, ident.1 = paste0( "male_", cellType ), 
                                  ident.2 = paste0( "female_", cellType ), 
                          group.by = "gender_celltype",  min.prt = 0.25, logfc.threshold = 0.25, 
                          assay = "integrated")
    }
  }
  
  degs = cbind( symbols = rownames(degs), degs)
  if( assay == "RNA" ) {
    write_xlsx( degs, paste0( "maleVSfemale_DEGs_full_prt0.1_logfc0.1_", cellType, ".xlsx" ) )
  } else {
    write_xlsx( degs, paste0( "maleVSfemale_DEGs_integrated_prt0.25_logfc0.25_", cellType, ".xlsx" ) )
  }
  
  xlim <- ceiling( max(abs(degs$avg_log2FC[is.finite(degs$avg_log2FC)])) )
  maxPvalue <- max( degs$p_val[degs$p_val_adj < 0.05] )
  keyvals <- ifelse( (degs$avg_log2FC < -1), 'royalblue',
                     ifelse((degs$avg_log2FC > 1), 'red', 'gray'))
  keyvals[is.na(keyvals)] <- 'gray'
  names(keyvals)[keyvals == 'red'] <- 'upReg'
  names(keyvals)[keyvals == 'gray'] <- 'noChange'
  names(keyvals)[keyvals == 'royalblue'] <- 'downReg'

  degs$log10p = -log10( degs$p_val )
  minPvalue <- ceiling( max( degs$log10p[ is.finite(degs$log10p) ] ) )
  degs$log10p = ifelse( degs$log10p > minPvalue, minPvalue, degs$log10p )
  degs$sigFlag = ifelse( abs(degs$avg_log2FC)>1 & degs$p_val_adj < 0.05, 1, 0 )
  
  EnhancedVolcano( degs,  lab = degs$symbols,
                   x = 'avg_log2FC', y = 'p_val',
                   xlab = 'Log2-fold change',
                   subtitle = bquote( italic(.(paste0( cellType, ": Male vs. Female" ))) ),
                   FCcutoff = 1, pCutoff = maxPvalue,
                   labSize = 3.0, colAlpha = 0.7,
                   colCustom = keyvals,
                   xlim = c( -xlim, xlim ),
                   pointSize = ifelse( degs$sigFlag == 1, 3, 1),
                   legendPosition = 'right', legendLabSize = 12, legendIconSize = 4.0,
                   drawConnectors = TRUE,  min.segment.length = 5, widthConnectors = 0.5, 
                   typeConnectors = "open", colConnectors = "black")
  
  if( assay == "RNA" ) {
    ggsave( paste0( "VolcanoPlot_maleVSfemale_DEGs_full_prt0.1_logfc0.1_", cellType, ".pdf" ), 
            width = 10, height = 10 )
  } else {
    ggsave( paste0( "VolcanoPlot_maleVSfemale_DEGs_integrated_prt0.25_logfc0.25_", cellType, ".pdf" ), 
            width = 10, height = 10 )
  }
  
  return( degs )
}


orig.ident = OSCC.nc.integrated@meta.data$orig.ident
genderTxt =   c( "female", "female", "female", "male", "male", "male" )
gender = rep( "", length.out = length(orig.ident))
for( i in 1:length(sampleNames) ) {
  idx = which( orig.ident == sampleNames[i] )
  gender[idx] = genderTxt[i]
}
celltype = OSCC.nc.integrated@meta.data$celltype_merged
gender_celltype = paste( gender, celltype, sep = "_" )
OSCC.nc.integrated <- AddMetaData( object = OSCC.nc.integrated, metadata = gender, col.name = 'gender')
OSCC.nc.integrated <- AddMetaData( object = OSCC.nc.integrated, metadata = gender_celltype, col.name = 'gender_celltype')

celltype = unique( as.character( celltype ) )
for( c in celltype ) {
  deTable = CS_GetDEGs( OSCC.nc.integrated, c, assay = "RNA"  )
}








# add sample annotation, sample names as reference
#sampleNames = c( "FC", "FC2", "FC3", "FT", "FT2", "FT3", 
#                 "MC", "MC2", "MC3", "MTT", "MT2", "MT3" )

orig.ident = OSCC.integrated@meta.data$orig.ident
sampleNames = c( "FC", "FC2", "FC3", "FT", "FT2", "FT3", 
                 "MC", "MC2", "MC3", "MTT", "MT2", "MT3" )
genderTxt =   c( "female", "female", "female", "female", "female", "female",
                 "male",   "male",   "male",   "male",   "male",   "male" )
diagTxt =     c( "normal", "normal", "normal", "tumor", "tumor", "tumor",
                 "normal", "normal", "normal", "tumor", "tumor", "tumor" )
groupTxt =     c( "normal-female", "normal-female", "normal-female", "tumor-female", "tumor-female", "tumor-female",
                  "normal-male",   "normal-male",   "normal-male",   "tumor-male",   "tumor-male",   "tumor-male" )

gender = rep( "", length.out = length(orig.ident))
diagnosis = rep( "", length.out = length(orig.ident))
sampleGroup = rep( "", length.out = length(orig.ident))
for( i in 1:length(sampleNames) ) {
  idx = which( orig.ident == sampleNames[i] )
  gender[idx]      = genderTxt[i]
  diagnosis[idx]   = diagTxt[i]
  sampleGroup[idx] = groupTxt[i]
}

OSCC.integrated <- AddMetaData( object = OSCC.integrated, metadata = diagnosis, col.name = 'diagnosis')
OSCC.integrated <- AddMetaData( object = OSCC.integrated, metadata = gender, col.name = 'gender')
OSCC.integrated <- AddMetaData( object = OSCC.integrated, metadata = newCellTypes, col.name = 'celltype' )
OSCC.integrated <- AddMetaData( object = OSCC.integrated, metadata = sampleGroup, col.name = 'sampleGroup' )

g2 <- DimPlot(OSCC.integrated, reduction = "umap", label = TRUE, repel = TRUE, ncol = 3, 
              pt.size = 2, label.size = 4, split.by = "orig.ident") + NoLegend()
ggsave( "OSCC_integrated_CellTypes_splitbySample_UMAP.PDF", g2, width = 15, height = 20, limitsize = FALSE )

g2 <- DimPlot(OSCC.integrated, reduction = "umap", label = TRUE, repel = TRUE, ncol = 2, 
              pt.size = 2, label.size = 4, split.by = "orig.ident", group.by = "clusterIDs") + NoLegend()
ggsave( "OSCC_integrated_CellTypes_splitbySample_clusterIDs_UMAP.PDF", g2, width = 15, height = 20, limitsize = FALSE )

g2 <- DimPlot(OSCC.integrated, reduction = "tsne", label = TRUE, repel = TRUE, ncol = 2, 
              pt.size = 2, label.size = 4, split.by = "orig.ident", group.by = "seurat_clusters") + NoLegend()
ggsave( "OSCC_integrated_CellTypes_splitbySample_clusterIDs_TSNE.PDF", g2, width = 20, height = 20, limitsize = FALSE )

g2 <- DimPlot(OSCC.integrated, reduction = "umap", label = TRUE, repel = TRUE, 
              pt.size = 2, label.size = 4, split.by = "diagnosis") + NoLegend()
ggsave( "OSCC_integrated_CellTypes_splitbyDiagnosis_UMAP.PDF", g2, width = 20, height = 10, limitsize = FALSE )
g2 <- DimPlot(OSCC.integrated, reduction = "umap", label = TRUE, repel = TRUE, 
              pt.size = 2, label.size = 4, split.by = "gender" ) + NoLegend()
ggsave( "OSCC_integrated_CellTypes_splitbyGender_UMAP.PDF", g2, width = 20, height = 10, limitsize = FALSE )

g2 <- DimPlot(OSCC.integrated, reduction = "umap", label = TRUE, repel = TRUE, ncol = 2,
              pt.size = 2, label.size = 4, split.by = "sampleGroup" ) + NoLegend()
ggsave( "OSCC_integrated_CellTypes_splitbySampleGroup_UMAP.PDF", g2, width = 20, height = 20, limitsize = FALSE )

p2 <- DimPlot(OSCC.integrated, reduction = "umap", pt.size = 2, label = TRUE, repel = TRUE )
ggsave("OSCC_integrated_TumorCells_UMAP.PDF", p2, width = 45, height = 20, limitsize = FALSE )
g2 <- DimPlot(OSCC.integrated, reduction = "umap", label = TRUE, repel = TRUE, ncol = 2,
              pt.size = 2, label.size = 4, split.by = "sampleGroup" ) + NoLegend()
ggsave( "OSCC_integrated_CellTypes_splitbySampleGroup_tumorCells_UMAP.PDF", g2, width = 20, height = 20, limitsize = FALSE )


g2 <- DimPlot(OSCC.integrated, reduction = "umap", label = TRUE, repel = TRUE, ncol = 2,
              pt.size = 2, label.size = 4, split.by = "sampleGroup" ) 
#+ NoLegend()
ggsave( "OSCC_integrated_CellTypes_splitbySampleGroup_tumorCells_UMAP_test.PDF", g2, width = 40, height = 20, limitsize = FALSE )

saveRDS(OSCC.integrated, file = "OSCC.integrated.rds")





#
# stoped here.

# need to set default assay to RNA. 
DefaultAssay( OSCC.integrated ) = "RNA"
g2 = FeaturePlot( OSCC.integrated, max.cutoff = 'q80', min.cutoff = 'q20', features = c("Cd33", "Cnn1", "Mpz", "Cxcl2",
                                                                     "Cd34", "Tnc"), split.by = "orig.ident" )
ggsave( "OSCC_integrated_FeatureMap_6genes.PDF", g2, width = 13, height = 15, limitsize = FALSE )

FeaturePlot( OSCC.integrated, features = "Clec3b", split.by = "orig.ident"  )

ggsave("OSCC_Ptprc_UMAP.PDF", width = 30, height = 15)

g2 <- DimPlot(OSCC.integrated, reduction = "umap", label = TRUE, repel = TRUE, ncol = 2,
              pt.size = 3, label.size = 8, split.by = "orig.ident") + NoLegend()
ggsave( "OSCC_integrated_CellTypes_splitbySamples_UMAP.PDF", g2, width = 20, height = 20, limitsize = FALSE )


a = OSCC.integrated@active.ident
b = grep('^ft', a, value = TRUE)
sub_obj <- subset(OSCC.integrated, idents = unique( b))
g2 <- DimPlot(sub_obj, reduction = "umap", label = TRUE, repel = TRUE, ncol = 2,
              pt.size = 3, label.size = 8, split.by = "orig.ident")




mData = OSCC.integrated@meta.data
mData$lungCellTypes = OSCC.integrated@active.ident
OSCC.integrated@meta.data = mData

# Get DEGs into a file.
CS_GetDEGs <- function( scData, cellType, sName1, sName2, assay = "integrated" ) {
  # create matrix to store the results. (each DEG list contains 5 columns, 4 samples, +1 counts DEGs)
  # if assay is RNA, we have to handle differently.
  
  # we do total of 1 comparison (male vs female).
  if( assay == "integrated" ) {
    numGenes = nrow( scData )
    DE.table = data.frame( matrix( rep(0, numGenes * 6), nrow = numGenes ) )
    rownames( DE.table ) <- rownames( scData )
  } else {
    numGenes = nrow( scData@assays$RNA$counts)
    DE.table = data.frame( matrix( rep(0, numGenes * 6), nrow = numGenes ) )
    rownames( DE.table ) <- rownames( scData@assays$RNA$counts )
  }
  colnames( DE.table )[1] <- "Counts"
  
  mData = scData@meta.data
  if( length( which( mData$celltype == paste0( tolower(sName1), "_", cellType ) )) > 0 &
      length( which( mData$celltype == paste0( tolower(sName2), "_", cellType ) )) > 0 ) {

    if( assay == "RNA" ) {
      degs = FindMarkers( scData, ident.1 = paste0( tolower(sName1), "_", cellType ), 
                          ident.2 = paste0( tolower(sName2), "_", cellType ), 
                          group.by = "celltype",  min.prt = 0.1, logfc.threshold = 0.1, assay = "RNA" )
    } else { 
      degs = FindMarkers( scData, ident.1 = paste0( tolower(sName1), "_", cellType ), 
                          ident.2 = paste0( tolower(sName2), "_", cellType ), 
                          group.by = "celltype",  min.prt = 0.25, logfc.threshold = 0.25)
    }
    colnames( DE.table )[2:6] <- paste0( sName1, "vs", sName2, "_", colnames(degs) )
    DE.table[rownames(degs), 1] = DE.table[rownames(degs), 1] + 1
    DE.table[rownames(degs), 2:6] = degs
  }
  
  # dx = which( DE.table$Counts == 0 )
  # DE.table = DE.table[-dx, ]
  # for( i in 2:21 ) {
  #   DE.table[ which( DE.table[,i] == 0 ), i] = ""
  # }
  DE.table = cbind( symbols = rownames(DE.table), DE.table)
  DE.table = arrange( DE.table, -Counts )
  if( assay == "RNA" ) {
    write_xlsx( DE.table, paste0( sName1, "vs", sName2, "_DEGs_full_prt0.1_logfc0.1_", cellType, ".xls" ) )
  } else {
    write_xlsx( DE.table, paste0( sName1, "vs", sName2, "_DEGs_integrated_prt0.25_logfc0.25_", cellType, ".xls" ) )
  }
  return( DE.table )
}

# Get DEGs into a file.
CS_GetDEGs_byCluster <- function( scData, cid1, cid2, assay = "integrated" ) {
  # create matrix to store the results. (each DEG list contains 5 columns, 4 samples, +1 counts DEGs)
  # if assay is RNA, we have to handle differently.
  
  # we do total of 4 comparisons.
  if( assay == "integrated" ) {
    numGenes = nrow( scData )
    DE.table = data.frame( matrix( rep(0, numGenes * 6), nrow = numGenes ) )
    rownames( DE.table ) <- rownames( scData )
  } else {
    numGenes = nrow( scData@assays$RNA$counts)
    DE.table = data.frame( matrix( rep(0, numGenes * 6), nrow = numGenes ) )
    rownames( DE.table ) <- rownames( scData@assays$RNA$counts )
  }
  colnames( DE.table )[1] <- "Counts"
  
  mData = scData@meta.data
  if( length( which( mData$clusterIDs == cid1 )) > 0 &
      length( which( mData$clusterIDs == cid2 )) > 0 ) {
    
    if( assay == "RNA" ) {
      degs = FindMarkers( scData, ident.1 = cid1, ident.2 = cid2, 
                          group.by = "clusterIDs",  min.prt = 0.1, logfc.threshold = 0.1, assay = "RNA" )
    } else { 
      degs = FindMarkers( scData, ident.1 = cid1, ident.2 = cid2, 
                          group.by = "clusterIDs",  min.prt = 0.25, logfc.threshold = 0.25)
    }
    colnames( DE.table )[2:6] <- paste0( cid1, "vs", cid2, "_", colnames(degs) )
    DE.table[rownames(degs), 1] = DE.table[rownames(degs), 1] + 1
    DE.table[rownames(degs), 2:6] = degs
  }
  
  # dx = which( DE.table$Counts == 0 )
  # DE.table = DE.table[-dx, ]
  # for( i in 2:21 ) {
  #   DE.table[ which( DE.table[,i] == 0 ), i] = ""
  # }
  DE.table = cbind( symbols = rownames(DE.table), DE.table)
  DE.table = arrange( DE.table, -Counts )
  if( assay == "RNA" ) {
    write_xlsx( DE.table, paste0( cid1, "vs", cid2, "_DEGs_full_prt0.1_logfc0.1.xls" ) )
  } else {
    write_xlsx( DE.table, paste0( cid1, "vs", cid2, "_DEGs_integrated_prt0.25_logfc0.25.xls" ) )
  }
  return( DE.table )
}

deTable = CS_GetDEGs( OSCC.integrated, "Fibroblasts", "FC", "FT"  )


deTable = CS_GetDEGs( OSCC.integrated, "Fibroblasts", "FC", "FT", assay = "RNA"  )

# neutrophil vs macrophage
deTable = CS_GetDEGs_byCluster( OSCC.integrated, "FT_12", "FT_2", assay = "RNA"  )
deTable = CS_GetDEGs_byCluster( OSCC.integrated, "MTT_0", "FT_6", assay = "RNA"  )

# macrophage (normal) vs. macrophage (tumor)
deTable = CS_GetDEGs_byCluster( OSCC.integrated, "FC_18", "FT_12", assay = "RNA"  )
deTable = CS_GetDEGs_byCluster( OSCC.integrated, "MC_20", "MTT_0", assay = "RNA"  )

# fibroblast vs myofibroblast
deTable1 = CS_GetDEGs_byCluster( OSCC.integrated, "FC_1", "FT_10", assay = "RNA"  )
deTable2 = CS_GetDEGs_byCluster( OSCC.integrated, "MC_4", "MTT_3", assay = "RNA"  )

# endothelial (normal) vs. endothelial (tumor)
deTable = CS_GetDEGs_byCluster( OSCC.integrated, "FC_8", "FT_28", assay = "RNA"  )
deTable = CS_GetDEGs_byCluster( OSCC.integrated, "MC_3", "MTT_18", assay = "RNA"  )

# lymphatic endo (normal) vs. lym endo (tumor)
deTable = CS_GetDEGs_byCluster( OSCC.integrated, "FC_0", "FT_6", assay = "RNA"  )
deTable = CS_GetDEGs_byCluster( OSCC.integrated, "MC_0", "MTT_10", assay = "RNA"  )

marker_short_list_clean = SC_CleanMarkerList_Dedup( OSCC.integrated, marker_short_list )
DotPlot( object = OSCC.integrated, features = marker_short_list_clean, 
         cluster.idents=T, col.min = 0, scale.min = 0 ) + 
  theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6) )
ggsave( "OSCC_integrated_DotPlot_short_grouped.pdf", width = 25, height = 20)

cIDs = c( "FC_1", "FC_3", "FC_4", "FC_13", "FC_15", "FT_10", 
          "MC_4", "MC_5", "MC_7", "MC_9", "MC_12", "MC_17", "MTT_3")
a = subset( x = OSCC.integrated, idents = cIDs )
a <- ScaleData(a)
genes = c( "Mmp13", "Spp1", "Mmp3", "Col12a1", "Crabp1", "Mmp10", "Saa3", "Serpine2", "Tnc", "Thbs2", 
           "Col7a1", "Igfbp3", "Col1a1", "Cdh11", "Bgn", "Spon1", "Pla2g7", "Cacna1c", "Stc1", "Mgp", 
           "Ablim1", "Rcan2", "Ifi27l2a", "Nova1", "Aspn", "Tppp3", "C3", "Efemp1", "Ebf1", "Ogn", 
           "Ppp2r2b", "Cxcl13", "Ly6a", "Dpt", "Col14a1", "Opcml", "Ltbp4", "Gstm1", "Igfbp6", "Creb5" )
DoHeatmap(object = a, features = genes)
ggsave( "Fibroblast_vs_Myofibroblasts.pdf", width = 20, height = 10 )


rnk1 = deTable1$FC_1vsFT_10_avg_log2FC
rnk2 = deTable2$MC_4vsMTT_3_avg_log2FC
names(rnk1) = rownames( deTable1 )
names(rnk2) = rownames( deTable2 )

gsea1 = SC_GSEA( rnk1, speciesTax = "Mus musculus", mSigDBCat = "H", mSigDBsubCat = NULL ) 
nes = data.frame( pathways = gsea1$pathway, 
                  value = gsea1$NES, 
                  adjp = -log10(gsea1$padj), 
                  name = "female" )

gsea2 = SC_GSEA( rnk2, speciesTax = "Mus musculus", mSigDBCat = "H", mSigDBsubCat = NULL ) 
nes = rbind( nes, data.frame( pathways = gsea2$pathway, 
                              value = gsea2$NES, 
                              adjp = -log10(gsea2$padj), 
                              name = "male" ) )

gsea1 = gsea1[order(gsea1$NES, decreasing = FALSE ), ]
nes$pathways =factor(nes$pathways, levels=gsea1$pathway )

nes2 = data.frame( pathways = gsea1$pathway, matrix( 0, nrow = nrow(gsea1), ncol =  2) )
colnames(nes2)[2:3] <- c("female", "male")
rownames(nes2) = gsea1$pathway
nes2[gsea1$pathway, 2] = gsea1$NES
nes2[gsea2$pathway, 3] = gsea2$NES
nes2$pathways =factor(nes2$pathways, levels=gsea1$pathway )

gg <- ggplot( nes )+
  geom_segment( data = nes2,
                aes( x = male, y = pathways,
                     yend = pathways, xend = female ),
                color = "#aeb6bf",
                size = 2.5, alpha = .5) +
  xlab( "NES") + 
  geom_point( aes( x = value, y = pathways, color = name, size = adjp), 
              show.legend = TRUE ) +
  ggtitle( "Difference between Male and Female, expression program changes from Fibro to myofibro" ) + 
  theme_clean()
ggsave( "Fibroblast_to_Myofibroblast_difference_in_male_vs_female.pdf", width = 10, height = 7)


# combined results.
for( i in 1:length(sampleNames) ) {
  if( i == 1 ) {
    cell.types = read.xlsx( paste0( sampleNames[i], "_custom_Mm-celltypes.xlsx") )
  } else {
    cell.types = rbind( cell.types, read.xlsx( paste0( sampleNames[i], "_custom_Mm-celltypes.xlsx") ) )
  }
}
cellTypes = unique( cell.types$customTypes )
prct = matrix( rep(0, length( sampleNames) * length( cellTypes ) ), ncol = length( sampleNames) )
colnames( prct ) <- sampleNames
cell.percent = data.frame( cell_types = cellTypes, prct )
rownames( cell.percent ) = cellTypes

cell_counts = matrix( rep(0, length( sampleNames) * length( cellTypes ) ), ncol = length( sampleNames) )
colnames( cell_counts ) <- sampleNames
cell_counts = data.frame( cell_types = cellTypes, cell_counts )
rownames( cell_counts ) = cellTypes

for( i in 1:length(sampleNames) ) {
  tmp = read.xlsx( paste0( sampleNames[i], "_custom_Mm-celltypes.xlsx") )
  tmp1 = aggregate( tmp$percent, by = list(sctypeType = tmp$customTypes), FUN = sum )
  cell.percent[tmp1$sctypeType, sampleNames[i]] = tmp1$x
  
  tmp2 = aggregate( tmp$counts, by = list(sctypeType = tmp$customTypes), FUN = sum )
  cell_counts[tmp2$sctypeType, sampleNames[i]] = tmp2$x
}
write_xlsx( cell.percent, path = "OSCC_4samples_cell_percent.xlsx")
write_xlsx( cell_counts, path = "OSCC_4samples_cell_counts.xlsx")

cell.percent$cell_types = factor( cell.percent$cell_types, levels = cell.percent$cell_types)
cell.percent.long = pivot_longer( cell.percent, 2:5)

library( gplots )
library( ggpubr )
library( RColorBrewer )
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(length(cellTypes))
ggbarplot(cell.percent.long, "name", "value", title = "Proportion of cell types",
          fill = "cell_types", color = "cell_types", palette = mycolors, 
          xlab = "Samples", ylab = "Percent of cells")
ggsave( "OSCC_4sample_cell_percent.pdf", width = 8, height = 8 )




# Session Info.
sessionInfo()

