#' ---
#' title: "Single Cell RNA-seq Using Seurat, Tongue cancer, Dr. Shivani Ruparel"
#' author: "Yidong Chen"
#' date:  "Oct 2023"
#' file name: scRNA-Shinvani_TongueCancer_Clean_Analysis.R
#' ---
#' 
#install.packages( 'Seurat' )
#install.packages( 'patchwork' )

# for cell annotation.
#library(devtools)
#install_github("guokai8/rcellmarker")
#install.packages("openxlsx")
#install.packages("HGNChelper")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install( "scuttle")
#BiocManager::install("scRNAseq")
#BiocManager::install("singleCellTK")
##remotes::install_github("mojaveazure/seurat-disk")
#
#library( devtools )
#devtools::install_github("nicolerg/ssGSEA2")
#install_github("cmap/cmapR" )
#BiocManager::install("flowCore")
#install_github( "eddelbuettel/bh" )

library( dplyr)
library(tidyverse)
library( Seurat )
library( patchwork )
library( ggplot2 )

# for cell annotation
library( rcellmarker )
library( scSorter )

library( openxlsx)
library( SingleR )
library( celldex )
library( scuttle )
library( writexl )
library( ggpubr )

source( "../SC_CommonFunctions.R")

#remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = F)
suppressMessages(require(DoubletFinder))

setwd("/r_workspace/cheny8/Shivani_OSCC")

#' ***
#' ## First test: load public-domain data
#' If you prefer to load 10X Genomic output result, using the function Read10X() (following Seurat tutorial)
#' Here we read data matrix, commonly available through GEO website, or other single-cell databases
#' 

#' ***
#' ## Load data set provided (10X CellRanger output)
#' 
#' This the main program reading scRNA-seq data
markers_list <- list( "T cell"        = c("Cd3d", "Cd3g", "Cd3e", "Cd2", "Gzmb", "Nkg7", "Ccl5" ),
                      "monocyte"      = c( "Cd83", "Cd74", "Ccr2", "Ctss", "Fcgr2b" ),  
                      "DC"            = c( "Cytip", "Ccr7", "Lsp1", "Ly75", "Ccl17", "Ccl22" ),
                      "Immune cell"   = c( "Cy1bb", "Hbegf", "Pak1", "Ifi30", "Fn1", "Acp5" ),
                      "Epithelial"    = c( "Krt8", "Krt7", "Krt19", "Krt18", "Krt14", "Cldn4",  "Cldn3", "Cldn7", "Epcam" ), 
                      "neutrophils"   = c( "Csf3r", "Cxcl2", "Cxcl3", "Il1b", "Cd14", "S100a9", "Hdc" ),
                      "Schwann"       = c( "Sorbs1", "Lgi4", "Apod", "Cryab", "Fxyd1", "Ptn", "Plp1", "Adgrl3" ),
                      "smooth muscle" = c( "Acta2", "Myl9", "Myh11", "Mylk", "Tagln", "Abcc9", "Rgs5", "Notch3" ),
                      "myofibroblast" = c( "Fbn1", "Palld", "Sulf1", "Des" ),
                      "fibroblast"    = c( "Postn", "Col1a1", "Col1a2", "Col3a1", "Col6a1", "Col12a1", "Mgp", "Dcn" ),
                      "Endothelial"   = c( "Timp3", "Flt4", "Cldn5", "Egfl7", "Lyve1", "Mmrn1", "Rrox1", "Sema3d", "Sema6a",
                                           "Fabp4", "Ptprb", "Flt1", "Pecam1", "Cd34", "Ly6a" ),
                      "keratinocyte"  = c( "Perp", "Fabp5", "Krt5" ) )

markers_HNC_list <- list( "Epithelial"   = c("Epcam", "Trp63","Krt6a", "Krt6b", "Krt6c", "Krt5", "Sfn", "Krt14", "Krt15", "Krt16", "Krt17",  "Krt18", "Krtcap3" ),
                      "Endothelial"      = c( "Vwf", "Cdh5", "Sele", "Cd34", "Eng", "Pecam1" ),  
                      "T cell"            = c( "Cd3d", "Cd3e", "Cd3g", "Cd6", "Cd2"  ),
                      "B cell"            = c( "Slamf7", "Cd79a", "Cd79b", "Blnk", "Fcrl5", "Cd19"),
                      "NK cell"            = c( "Nkg7", "Ncam1", "Klrb1", "Klrc1", "Fcgr3"),
                      "Macrophages"       = c( "Cd14", "Cd163", "Cd68", "Fcgr2a", "Csf1r" ), 
                      "DC"                = c( "Cd40", "Cd80", "Cd83", "Ccr7", "Cd1c", "Clec4d", "IL3RA" ),
                      "Fibroblast"       = c( "Fap", "Thy1", "Pdpn", "Mmp2", "Mmp11", "Pdgfa", "Pdgfrl", "Tgfbr3", "Acta2", "Mcam", "Mylk", "MYL9", "Il6", "Pdgfra"),
                      "Mast"             = c( "Ms4a2", "Tpsb2", "Tpsab1", "Cma1" ) )

marker_custom_list_old<-  
  list( "Endothelial"               = c( "Pecam1", "Cd34", "Prom1", "Kdr", "Vwf", "Fabp4", "Itag6" ),
        "Pericytes"                 = c( "Pdgfrb", "Acta2", "Ng2", "Mcam", "Rgs5", "Tagln", "Ccdc102b", "Ccl11", "Cspg4"),
        "Smooth Muscle"             = c( "Myh11", "Acta2", "Tgln", "Cnn1" ),
        "Early Smooth Muscle"       = c( "Pdgfrb", "Acta2", "Tgln", "Mcam", "Rgs5", "Myl9", "Myl6", "Myo1b" ),
        "Skeletal Muscle"           = c( "Des", "Acta1", "Myog", "Myod1", "Myh1", "Mb", "Ckm"),
        "skeletal muscle satellite cells" = c( "Pax7", "Myf5",  "Des", "Tff3" ),
        "Embryonic myogenic progenitor"   = c( "Mef2a", "Pax3", "Ncam", "Lama2" ),
        "Fibroblasts"               = c( "Ncam1", "Col1a1", "Col1a2", "Col3a1", "Dcn", "Fn1", "S100a4", "Thy1", "Scx", "Ctsk", "Thbs1", "Timp1", "ccn2", "Tnc", "Cd44"),
        "myofibroblasts"            = c( "Dcn", "Acta2", "Cnn1", "Fn1" ),
        "keratinocyte"              = c("Krt5", "Krt14", "Cdh1", "Tp63", "Itgb4" ),
        "differentiating keratinocyte" = c(	"Krt1", "Krt10", "Krt18", "Cdh1" ),  
        "Mucosal Epithelial Cells"  = c( "Muc1", "Krt7", "Cldn1", "Cdh1", "Epcam"), 
        "Secretary Epithelial cell" = c("Epcam", "Cd24", "Muc5b", "Krt19", "Scgb3a", "Pigr"), 
        "Goblet Cells"              = c( "Muc2", "Muc5ac", "Krt20", "Cd24", "Tff3", "Tff1", "Muc5b" ),      
        "T cell"                    = c(  "Ptprc", "Cd3d", "Cd3g", "Cd3e", "Cd4","Cd8" ),
        "Macrophages"               = c( "Ptprc", "Cd163", "Cd68", "Clec4a", "Ccr2", "Fcgr3a", "Cd14"),
        "Schwann Cells"             = c( "Lgi4", "Ptn", "Plp1", "Mpz" )
  )

marker_custom_list <- 
  list( "Endothelial"               = c( "Pecam1", "Cd34", "Prom1", "Kdr", "Vwf", "Fabp4", "Itga6", "Egfl7", "Tie1" ),
        "Lymphatic endothelial"     = c( "Prox1", "Lyve1", "Flt4" ),
        "Pericytes"                 = c( "Pdgfrb", "Acta2", "Cspg4", "Mcam", "Rgs5", "Tagln", "Ccdc102a", "Ccl11" ),
        "Smooth Muscle"             = c( "Myh11", "Acta2", "Tagln", "Cnn1" ),
        "Early Smooth Muscle"       = c( "Pdgfrb", "Acta2", "Tagln", "Mcam", "Rgs5", "Myl9", "Myl6", "Myo1b" ),
        "Skeletal Muscle"           = c( "Des", "Acta1", "Mb", "Ckm"),
        "skeletal muscle satellite cells" = c( "Pax7", "Myf5",  "Des" ),
#        "Embryonic myogenic progenitor"   = c( "Mef2a", "Pax3", "Ncam1", "Lama2" ),
        "neural progenitor cell"    = c("Sox5", "Fgfr1", "Robo1", "Ncam1"),
        "Fibroblasts"               = c( "Ncam1", "Col1a1", "Col1a2", "Col3a1", "Dcn", "Fn1", "S100a4", "Thy1", "Scx", "Ctsk", "Thbs1", "Timp1", "Tnc", "Cd34"),
        "myofibroblasts"            = c( "Dcn", "Acta2", "Cnn1", "Fn1" ),
        "keratinocyte"              = c("Krt5", "Krt14", "Cdh1", "Trp63", "Itgb4" ),
        "differentiating keratinocyte" = c( "Krt10", "Krt18", "Cdh1" ), 
        "differentiating suprabasal keratinocyte" = c("Krt19", "Cdh11", "Klf5" ), 
        "Mucosal Epithelial Cells"  = c( "Muc1", "Krt7", "Cldn1", "Cdh1", "Epcam", "Tff3", "Muc5b"), 
        "Secretary Epithelial cell" = c("Epcam", "Cd24a", "Muc5b", "Krt19", "Scgb3a1", "Pigr"), 
        "Goblet Cells"              = c( "Muc5ac", "Cd24a", "Tff3", "Tff1", "Muc5b" ),      
        "T cell"                    = c(  "Ptprc", "Cd3d", "Cd3g", "Cd3e", "Cd4","Cd8a" ),
        "Macrophages"               = c( "Ptprc", "Cd163", "Cd68", "Clec4a2", "Ccr2", "Fcgr4", "Cd14"),
        "neutrophils"               = c( "Csf3r", "S100a8" ),
        "DC"                        = c( "Itgam", "Itgax", "Zbtb46", "Ly75", "Cd80", "Cd83"),
        "Schwann Cells"             = c( "Lgi4", "Ptn", "Plp1", "Mpz", "Scn7a"),
        "Mast"                      = c( "Ms4a2", "Tpsb2", "Mcpt4", "Mrgprb2" ),
        "EMT"                       = c( "Cdh1", "Cdh2", "Vim", "Snai1", "Twist1" ),
        "immune-like"               = c( "Cd274", "Ido1", "Gata3", "Cd44" )
  )

marker_short_list <- 
  list( "Endothelial"               = c( "Pecam1", "Prom1", "Kdr", "Vwf", "Tie1" ),
        "Lymphatic endothelial"     = c( "Prox1", "Lyve1", "Flt4" ),
        "Pericytes"                 = c( "Pdgfrb", "Rgs5", "Tagln",  "Cnn1" ),
        "Smooth Muscle"             = c( "Acta2", "Myh11", "Myl9", "Myl6", "Myo1b" ),
        "Skeletal Muscle"           = c( "Des", "Acta1", "Mb", "Ckm"),
        "skeletal muscle satellite cells" = c( "Pax7", "Myf5",  "Des" ),
        "Fibroblasts"               = c( "Col1a1", "Col1a2", "Col3a1", "Dcn", "Fn1" ),
        "Basal keratinocyte"        = c("Krt5", "Krt14", "Cdh1", "Itgb4" ),
        "differentiating keratinocyte" = c(	"Krt7", "Krt10", "Krt18" ),
        "differentiating suprabasal keratinocyte" = c("Krt19", "Cdh11", "Klf5" ), 
        "T cell"                    = c(  "Ptprc", "Cd3d", "Cd3g", "Cd3e" ),
        "Macrophages"               = c( "Cd163", "Cd68", "Ccr2", "Cd14"),
        "neutrophils"               = c( "S100a8", "Csf3r", "Cxcl2", "Cxlc3"),
        "DC"                        = c( "Zbtb46", "Ly75", "Cd80", "Cd83"),
        "Schwann Cells"             = c( "Lgi4", "Ptn", "Plp1", "Mpz"),
        "Mast"                      = c( "Ms4a2", "Tpsb2", "Mcpt4", "Mrgprb2" ), 
        "OSCC"                      = c( "Egfr", "Scca", "Kras", "Ccnd1")
  )

# making final figures. 
marker_OSCC_list <- 
  list( "OSCC"          = c("Krt7", "Krt19", "Krt14", "Cldn4",  "Cldn3", "Cldn7", "Epcam"),
        "Endothelial"   = c("Flt4", "Cldn5", "Egfl7", "Lyve1", "Fabp4", "Flt1", "Pecam1", "Cd34"),
        "Monocyte"      = c("Cd83", "Cd74", "Ccr2", "Ctss", "Fcgr2b", "Cd68", "Cybb", "Adgre1"),
        "T cell"        = c("Cd3d", "Cd3g", "Cd3e"),
        "DC"            = c("Ccr7", "Ccl17", "Ccl22"),
        "Neutrophils"   = c("Csf3r", "Cxcl2", "Cxcl3", "S100a9", "Hdc"),
        "Fibroblast"    = c("Col1a1", "Col1a2", "Col3a1", "Mgp", "Dcn"),
        "Myofibroblast" = c("Fbn1", "Palld", "Sulf1", "Des"),
        "Keratinocyte"  = c("Perp", "Fabp5", "Krt5"),
        "Smooth muscle" = c("Myh11", "Rgs5", "Notch3"),
        "Schwann cell"  = c("Lgi4", "Apod","Ptn", "Plp1") )


path2files = c( "/r_workspace/data/Shivani_OSCC/FC_Count/FC_Count_outs/filtered_feature_bc_matrix", 
                "/r_workspace/data/Shivani_OSCC/FT_Count/FT_Count_outs/filtered_feature_bc_matrix",
                "/r_workspace/data/Shivani_OSCC/MC_Count/MC_Count_outs/filtered_feature_bc_matrix",
                "/r_workspace/data/Shivani_OSCC/MTT_Count/MTT_Count_outs/filtered_feature_bc_matrix",
                "/r_workspace/data/Shivani_OSCC/FN1_count/FN1_count_outs/filtered_feature_bc_matrix", # new data.
                "/r_workspace/data/Shivani_OSCC/FN2_count/FN2_count_outs/filtered_feature_bc_matrix", # new data.
                "/r_workspace/data/Shivani_OSCC/F1_count/F1_count_outs/filtered_feature_bc_matrix", # new data.
                "/r_workspace/data/Shivani_OSCC/F2_count/F2_count_outs/filtered_feature_bc_matrix", # new data.
                "/r_workspace/data/Shivani_OSCC/MN1_count/MN1_count_outs/filtered_feature_bc_matrix", # new data.
                "/r_workspace/data/Shivani_OSCC/MN2_count/MN2_count_outs/filtered_feature_bc_matrix", # new data.
                "/r_workspace/data/Shivani_OSCC/M1_count/M1_count_outs/filtered_feature_bc_matrix", # new data.
                "/r_workspace/data/Shivani_OSCC/M2_count/M2_count_outs/filtered_feature_bc_matrix" # new data.
)
uppThr = c( 5000, 7500, 5000, 7500, 
            5000, 5000, 6000, 5000, 5000, 5000, 7000, 5000 )
mtThr =  c( 10, 10, 10, 10, 
            15, 15, 10, 10, 15, 15, 10, 10 )
sampleNames = c( "FC", "FT", "MC", "MTT", 
                 "FC2", "FC3", "FT2", "FT3", "MC2", "MC3", "MT2", "MT3" )

procParams = data.frame( matrix( nrow = 0, ncol = 13 ) )
colnames( procParams ) = c( "parameters", sampleNames )
procParams[1, 1] = "min.cells"
procParams[2, 1] = "min.features"
procParams[3, 1] = "num.cells"
procParams[4, 1] = "Feature.upper.threshold"
procParams[5, 1] = "min.mitochondrial"
procParams[6, 1] = "resolution"
procParams[7, 1] = "num.cells.filtered"
procParams[8, 1] = "num.cells.singlelets"
procParams[9, 1] = "num.clusters"
procParams[10, 1] = "markers.min.pct"
procParams[11, 1] = "markers.logfc.threshold"
procParams[12, 1] = "markers.num.reported"

# just update the new ones. 
updateSamples = c( 8, 10, 12)
for( i in updateSamples ) {
#for( i in 5:length( sampleNames ) ) {
  scData.raw <- Read10X(data.dir = path2files[i])  
  scData <- CreateSeuratObject( counts = scData.raw,  project = sampleNames[i], min.cells = 3, min.features = 200 )
  procParams[1, sampleNames[i]] = 3
  procParams[2, sampleNames[i]] = 200
  procParams[3, sampleNames[i]] = ncol( scData )
  
  dir.create( paste0( sampleNames[i], "_results" ) )
  samplePrefix = paste0( sampleNames[i], "_results/", sampleNames[i] )
  
  scData.norm.scaled <- SC_Mm_PreprocessData( scData, samplePrefix, nFeatureLowerThr = 200, nFeatureUpperThr = uppThr[i], 
                                              mtThr = mtThr[i], resThr = 1.5 )
  saveRDS(scData.norm.scaled, file = paste0( sampleNames[i], "_results/", tolower(sampleNames[i]), ".norm.scaled.rds" ) )
  #saveRDS(scData.norm.scaled, file = paste0( sampleNames[i], "_results/", tolower(sampleNames[i]), ".norm.scaled.rds" ) )
  
  procParams[4, sampleNames[i]] = uppThr[i]
  procParams[5, sampleNames[i]] = mtThr[i]
  procParams[6, sampleNames[i]] = 1.5
  procParams[7, sampleNames[i]] = ncol( scData.norm.scaled )
  
  scData.singlets = SC_RemoveDoublets( scData.norm.scaled, samplePrefix )
  # repeat QC process.
  scData.norm.scaled <- SC_Mm_PreprocessData( scData.singlets, samplePrefix, nFeatureLowerThr = 200, nFeatureUpperThr = uppThr[i], 
                                              mtThr = mtThr[i], resThr = 1.5, qcFigName = "Preprocessing_QCPlots_2nd" )
  saveRDS(scData.norm.scaled, file = paste0( sampleNames[i], "_results/", tolower(sampleNames[i]), ".norm.scaled.singlets.rds" ) )
  #scData.norm.scaled = readRDS( paste0( sampleNames[i], "_results/", tolower(sampleNames[i]), ".norm.scaled.singlets.rds" ) )
  
  SC_Preprocessing_ExtraOutput( scData.norm.scaled, 'mouse', samplePrefix )
  
  procParams[8, sampleNames[i]] = ncol( scData.norm.scaled )
  procParams[9, sampleNames[i]] = nlevels( scData.norm.scaled )
  procParams[10, sampleNames[i]] = 0.25
  procParams[11, sampleNames[i]] = 0.25
  procParams[12, sampleNames[i]] = 50
  
  numClusters <- nlevels( scData.norm.scaled )
  scData.markers <- FindAllMarkers( scData.norm.scaled, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25 )
  saveRDS( scData.markers, file = paste0( sampleNames[i], "_results/", tolower(sampleNames[i]), ".markers.rds") )
  #scData.markers = readRDS( paste0( sampleNames[i], "_results/", tolower(sampleNames[i]), ".markers.rds") )
  
  SC_MarkerGenes_Processing( scData.norm.scaled, scData.markers, samplePrefix, numDEGs = 50, species = "mouse" )
  
  markers_list_clean = SC_CleanMarkerList( scData.norm.scaled, markers_list )
  cell.types = SC_Mm_CellTypeAnnotation( scData.norm.scaled, markers_list_clean, scData.markers, paste0( samplePrefix, "_enrichr_mouse" )  )
  lungMap_cluster_ids <- cell.types$sctypeType
  SC_MarkerGenes_Plots( scData.norm.scaled, markers_list_clean, lungMap_cluster_ids, paste0( samplePrefix, "_enrichr_mouse" ),
                         heightAt = 8, widthAt = 40 )
  
  markers_HNC_list_clean = SC_CleanMarkerList( scData.norm.scaled, markers_HNC_list )
  cell.types = SC_Mm_CellTypeAnnotation( scData.norm.scaled, markers_HNC_list_clean, scData.markers, paste0( samplePrefix, "_HNC_mouse" )  )
  HNC_cluster_ids <- cell.types$sctypeType
  SC_MarkerGenes_Plots( scData.norm.scaled, markers_HNC_list_clean, HNC_cluster_ids, paste0( samplePrefix, "_HNC_mouse" ),
                        heightAt = 8, widthAt = 40 )

  marker_custom_list_clean = SC_CleanMarkerList( scData.norm.scaled, marker_custom_list )
  cell.types = SC_Mm_CellTypeAnnotation( scData.norm.scaled, marker_custom_list_clean, scData.markers, paste0( samplePrefix, "_custom" )  )
  custom_cluster_ids <- cell.types$sctypeType
  SC_MarkerGenes_Plots( scData.norm.scaled, marker_custom_list_clean, custom_cluster_ids,  paste0( samplePrefix, "_custom" ),  
                        heightAt = 8, widthAt = 40 )

  marker_OSCC_list_clean = SC_CleanMarkerList( scData.norm.scaled, marker_OSCC_list )
  cell.types = SC_Mm_CellTypeAnnotation( scData.norm.scaled, marker_OSCC_list_clean, scData.markers, paste0( samplePrefix, "_OSCC" )  )
  custom_cluster_ids <- cell.types$sctypeType
  SC_MarkerGenes_Plots( scData.norm.scaled, marker_OSCC_list_clean, custom_cluster_ids,  paste0( samplePrefix, "_OSCC" ),  
                        heightAt = 8, widthAt = 40 )
}
write_xlsx( procParams, path = "OSCC_procParams_updated.xlsx"  )

# split cluster 16 (of MC sample) into two groups. 
scData.norm.scaled = readRDS( paste0( "mc.norm.scaled.singlets.rds" ) )
scData.norm.scaled.c16 <- subset( scData.norm.scaled, idents = 16  )
scData.norm.scaled.c16 <- FindNeighbors(scData.norm.scaled.c16, dims = 1:50, k.param = 3)
scData.norm.scaled.c16 <- FindClusters(scData.norm.scaled.c16, resolution = 0.1, verbose = FALSE)
DimPlot( scData.norm.scaled.c16, reduction = "umap")
dx = which( scData.norm.scaled.c16@active.ident == 1 )
cNames = names( scData.norm.scaled.c16@active.ident[dx] )

new_levels <- c(levels(scData.norm.scaled@active.ident), "26")
scData.norm.scaled@active.ident <- factor(scData.norm.scaled@active.ident, levels = new_levels)
scData.norm.scaled@active.ident[cNames] = 26
scData.norm.scaled@meta.data$seurat_clusters = scData.norm.scaled@active.ident
saveRDS( scData.norm.scaled, file = "mc.norm.scaled.singlets.rds" )
DimPlot( scData.norm.scaled, reduction = "umap")

#
#
# fast update (shorter code comparing to above loop)
for( i in 10:length( sampleNames ) ) {
  samplePrefix = paste0( sampleNames[i], "_results/", sampleNames[i] )
  
#  scData.norm.scaled = readRDS( paste0( sampleNames[i], "_results/", tolower(sampleNames[i]), ".norm.scaled.singlets.rds" ) )
  scData.norm.scaled = readRDS( paste0( tolower(sampleNames[i]), ".norm.scaled.singlets.rds" ) )
  
  #SC_Preprocessing_ExtraOutput( scData.norm.scaled, 'mouse', samplePrefix )
  
  numClusters <- nlevels( scData.norm.scaled )
#  scData.markers <- FindAllMarkers( scData.norm.scaled, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25 )
#  saveRDS( scData.markers, file = paste0( sampleNames[i], "_results/", tolower(sampleNames[i]), ".markers.rds") )
#  scData.markers = readRDS( paste0( sampleNames[i], "_results/", tolower(sampleNames[i]), ".markers.rds") )
  scData.markers = readRDS( paste0( tolower(sampleNames[i]), ".markers.rds") )
  #SC_MarkerGenes_Processing( scData.norm.scaled, scData.markers, samplePrefix, numDEGs = 50, species = "mouse" )
  
  marker_custom_list_clean = SC_CleanMarkerList( scData.norm.scaled, marker_custom_list )
  cell.types = SC_Mm_CellTypeAnnotation( scData.norm.scaled, marker_custom_list_clean, scData.markers, paste0( samplePrefix, "_custom" )  )
  custom_cluster_ids <- cell.types$sctypeType
  SC_MarkerGenes_Plots( scData.norm.scaled, marker_custom_list_clean, custom_cluster_ids,  paste0( samplePrefix, "_custom" ),  
                        heightAt = 8, widthAt = 40 )
  
  marker_short_list_clean = SC_CleanMarkerList( scData.norm.scaled, marker_short_list )
  cell.types = SC_Mm_CellTypeAnnotation( scData.norm.scaled, marker_short_list_clean, scData.markers, paste0( samplePrefix, "_short" )  )
  custom_cluster_ids <- cell.types$sctypeType
  SC_MarkerGenes_Plots( scData.norm.scaled, marker_short_list_clean, custom_cluster_ids,  paste0( samplePrefix, "_short" ),  
                        heightAt = 8, widthAt = 25 )
}

marker.genes = c( "Krt7", "Egfr", "Kras", "Ccnd1", "Cdh2", "Vim", "Twist1", "Snai1", "Ido1", "Cd55" )
for( i in 1:length( sampleNames ) ) {
  samplePrefix = paste0( sampleNames[i], "_results/", sampleNames[i] )
  
  scData.norm.scaled = readRDS( paste0( tolower(sampleNames[i]), ".norm.scaled.singlets.rds" ) )

  DotPlot( scData.norm.scaled, features = marker.genes, col.min = 0.1, dot.min = 0 ) +  RotatedAxis()
  ggsave(  paste0( samplePrefix, "_DotPlot_OSCC_shortList.pdf" ), width = 6, height = 8 )
}



# using inferCNV, other than install infercnv, you may need to install package JAGS (Just Another Gibbs Sampler Files)
# from https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/
#BiocManager::install( "infercnv")
library(infercnv)

# use endothelial cells as reference.
ref_groups = c( "0", "6", "0", "10",
                "2", "1", "15", "13",
                "0", "0", "10", "13")
for( i in 5:length( sampleNames ) ) {
  scData.norm.scaled = readRDS( paste0( tolower(sampleNames[i]), ".norm.scaled.singlets.rds" ) )
  meta.data <- scData.norm.scaled@meta.data[, 9, drop = F]
  
  annotFileName = paste0( "annotation_", tolower(sampleNames[i]),  ".csv")
  write.table( meta.data, annotFileName, sep=",", quote=F, col.names=F )
  
  # infercnv object.
  infercnv_obj <- CreateInfercnvObject( raw_counts_matrix= scData.norm.scaled@assays$RNA@counts,
                                        annotations_file= annotFileName,
                                        delim=",",
                                        gene_order_file="../mm10_gencode_vM25_pos.txt",
                                        ref_group_names= ref_groups[i])
  
  infercnv_obj = infercnv::run(infercnv_obj = infercnv_obj, 
                               cutoff=0.1, out_dir=paste0( sampleNames[i], "_results/" ), 
                               cluster_by_groups=T,  denoise=T,
                               num_threads=8, HMM=F )
}


# feature plot.
featureMarkers = c("Cd3d", "Fcgr2b", "Ccr7", "Mki67", "Fn1", "Spp1", "Cxcl2", "Ptn", "Rgs5", "Dcn", "Fabp4", "Krt5")
FeaturePlot( tongue.mtt.norm.scaled, features = featureMarkers )
ggsave( "tongue.mtt_featurePlot_selectedMarkers.pdf", width = 15, height = 10)

# making final figures. 

new.cluster.ids <- c("OSCC",          "Monocyte/Mac",    "Fibroblast",   "Monocyte/Mac",  "OSCC", 
                     "OSCC",          "Endothelial",     "Neutrophils",  "Neutrophils",   "OSCC",
                     "Smooth muscle", "OSCC",            "Monocyte/Mac", "T cell",        "Lymph vessel endothelial", 
                     "Monocyte/Mac",   "Schwann cells",   "DC",           "Myofibroblast", "Keratinocytes" )
names(new.cluster.ids) <- levels(tongue.mtt.norm.scaled)
tongue_tmp <- RenameIdents(tongue.mtt.norm.scaled, new.cluster.ids)
p1 = DimPlot(tongue.mtt.norm.scaled, reduction = "umap", label = TRUE, pt.size = 0.5)
p2 = DimPlot(tongue_tmp, reduction = "umap", label = TRUE, pt.size = 0.5)
p1 + p2
ggsave( "tongue.mtt_UMap_labeled.pdf", width = 15, height = 7)

# generate average expression
avgExpr = AverageExpression( tongue.mtt.norm.scaled )
out_data = data.frame( gene_name = rownames(avgExpr$RNA), avgExpr$RNA )
colnames(out_data)[2:21] = colnames(avgExpr$RNA)
write.table( out_data, file = "tongue_mtt_avgExpr_byCluster.txt", 
             row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE )
write_xlsx( out_data, path = "tongue_mtt_avgExpr_byCluster.xlsx" )



# Cell type annotation. We will try multiple packages.
# 1) CellMarker
res <- cellMarker(hcc_Cntl.markers, type = "seurat", species = "human", 
                  keytype = "SYMBOL", weight = 1, padj = 0.05,  
                  tissue = "liver", minSize = 3)
cell.types = data.frame( cluster = unique(hcc_Cntl.markers$cluster))
cellmarkerType = rep( "unassigned", nrow( cell.types))
cellmarkerGenes = rep("", nrow( cell.types))
cellmarkerPadj = rep(1, nrow( cell.types))
for( i in 1:nrow(cell.types) ) {
  idx = which( res$cluster == i )
  if( length( idx ) > 0 ) {
    cellmarkerType[i] = res$cellType[idx[1]]
    cellmarkerGenes[i] = res$GeneID[idx[1]]
    cellmarkerPadj[i] = res$Padj[idx[1]]
  }
}
cell.types$cellmarkerType = cellmarkerType
cell.types$cellmarkerGenes = cellmarkerGenes
cell.types$cellmarkerPadj = cellmarkerPadj

# 2) scType
library( HGNChelper )
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); 
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# get cell-type-specific gene sets from our in-built database (DB)
# e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
gs_list = gene_sets_prepare( db_,  "Liver" ) 

# assign cell types
es.max = sctype_score(scRNAseqData = hcc.Cntl.norm.scaled[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative )

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(hcc.Cntl.norm.scaled@meta.data$seurat_clusters), function(cl){
  metaData <- hcc.Cntl.norm.scaled@meta.data
  es.max.cl = sort( rowSums( es.max[ ,rownames(metaData[metaData$seurat_clusters==cl, ])]), decreasing = !0 )
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, 
                  ncells = sum(metaData$seurat_clusters==cl)), 10 )
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
sctype_scores = sctype_scores[order(sctype_scores$cluster), ]
colnames(sctype_scores) = c("cluster", "sctypeType", "sctypeScore", "sctypeNcells")
cell.types = cbind( cell.types, sctype_scores[,2:4])
write.table( cell.types, file = "HCC-Cntl-celltypes.txt", row.names = FALSE, sep = "\t" )


hcc.1N.norm.scaled <- SetIdent(hcc.1N.norm.scaled, value = "RNA_snn_res.1")
new.cluster.ids <- cell.types$sctypeType
names(new.cluster.ids) <- levels(hcc.1N.norm.scaled)
hcc.1N.norm.scaled <- RenameIdents(hcc.1N.norm.scaled, new.cluster.ids)
hcc.1N.norm.scaled[["sctype_liver"]] <- Idents(object = hcc.1N.norm.scaled)
DimPlot(hcc.1N.norm.scaled, label = T , repel = T, label.size = 3) + NoLegend()
ggsave( "hcc.1N_UMap_sctype.pdf", width = 7, height = 5 )


# using SingleR
# 1) using human primary cell atlas
hpca.se <- HumanPrimaryCellAtlasData()
hcc.Cntl.sce = as.SingleCellExperiment(hcc.Cntl.norm.scaled )
hcc.Cntl <- logNormCounts(hcc.Cntl.sce)
pred.Cntl <- SingleR(test = hcc.Cntl, ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.main )
table(pred.Cntl$labels)
hcc.13N.norm.scaled@meta.data$hpca.pruned <- pred.13N$pruned.labels
hcc.13N.norm.scaled <- SetIdent(hcc.13N.norm.scaled, value = "hpca.pruned")
DimPlot(hcc.13N.norm.scaled, label = T , repel = T, label.size = 3) + NoLegend()
ggsave( "hcc.13N_UMap_HPCA.pdf", width = 7, height = 5 )

# 2) using HCC data (the code works, but the results are not correct) 
library(scRNAseq)
#library( singleCellTK )
sceLiver = ZhaoImmuneLiverData()
sceLiver <- sceLiver[,!is.na(sceLiver@colData$broad)]
sceLiver <- logNormCounts(sceLiver)

# sceLiver data set with gene id in ensembl, need to change.
source( "convertGeneIDs.R" )
library("org.Hs.eg.db")
sceLiver <- convertGeneIDs(inSCE = sceLiver,
                           inSymbol = "ENSEMBL", outSymbol = "SYMBOL",
                           database = "org.Hs.eg.db")

pred.1T.liver <- SingleR(test=hcc.1T, ref=sceLiver, labels=sceLiver@colData$broad, de.method="wilcox")
table(pred.1T.liver$labels)

hcc.1T.norm.scaled@meta.data$hcc_zhao.pruned <- pred.1T.liver$pruned.labels
hcc.1T.norm.scaled <- SetIdent(hcc.1T.norm.scaled, value = "hcc_zhao.pruned")
DimPlot(hcc.1T.norm.scaled, label = T , repel = T, label.size = 3)

# We try to load another HCC dataset
library(SeuratDisk)
Convert("HCC.h5ad", "h5seurat")
HCC <- LoadH5Seurat("HCC.h5Seurat")
DimPlot(HCC, label = T , repel = T, label.size = 3)






# 3) scSorter using cell markers
library( scSorter )

#annot = openxlsx::read.xlsx( "/Users/cheny8/Library/CloudStorage/Box-Box/ChenLab_Bioinformatics/BioinfoServices/SingleCellProjects/CellTypeMarkerSets/Cell_marker_Human_v2.0.xlsx")
annot = openxlsx::read.xlsx( "Cell_marker_Human_v2.0.xlsx")
idx = which( annot$Tissue.type == "Brain")
brain_annot = annot[idx, ]
brain_marker = brain_annot[, c(5,6,7,10,11,13)]
cn <- colnames( brain_marker)
cn[3] <- "Type"
cn[4] <- "Marker"
colnames(brain_marker ) <- cn

# select genes for the matrix
topgenes <- head(VariableFeatures(bt40.hg.norm), 2000)
expr = GetAssayData(bt40.hg.norm)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]

picked_genes = unique( c(brain_marker$Marker, topgenes) )
expr = expr[rownames(expr) %in% picked_genes, ]

# remove genes that are not matched with expr matrix
dx <- brain_marker$Marker %in% rownames( counts.hg )
idx <- which( dx == TRUE )
brain_marker <- brain_marker[idx, ]

# Running scSorter (this command takes about 6 hours,  be aware)
rts <- scSorter(expr, brain_marker)
print(table(rts$Pred_Type))
table(bt40.hg.norm$seurat_clusters, rts$Pred_Type)

# using inferCNV, other than install infercnv, you may need to install package JAGS (Just Another Gibbs Sampler Files)
# from https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/
#BiocManager::install( "infercnv")
library(infercnv)
meta.data <- hcc.13N.norm.scaled@meta.data[, 8, drop = F]
write.table( meta.data, "annotation.csv", sep=",", quote=F, col.names=F )

# infercnv object.
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix= hcc.13N.norm.scaled@assays$RNA@counts,
                                     annotations_file= "annotation.csv",
                                     delim=",",
                                     gene_order_file="../gencode_v19_gene_pos.txt",
                                     ref_group_names=c('6'))

infercnv_obj = infercnv::run(infercnv_obj = infercnv_obj, 
                             cutoff=0.1, out_dir="./", cluster_by_groups=T,  denoise=T,
                             num_threads=8, HMM=F )



# Session Info.
sessionInfo()
