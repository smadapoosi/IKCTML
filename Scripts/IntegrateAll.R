library(Seurat)
args <- commandArgs(TRUE)

sobjm <- readRDS(args[1])

rename <- function(sobj){
  sobj[['Leveled_Names']] <- rep(NA, length(sobj$nCount_RNA))
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('Glomerular_endothelium_Young','EC-AEA_Menon','Endothelial Cells - AEA & DVR _Lake','GC-EC_Menon','Descending_vasa_recta_Young', 'Endothelial Cells - glomerular capillaries_Lake','EC-PT_Menon','Ascending_vasa_recta_Young','Endothelial Cells - AVR_Lake','EC_Wu','Endothelial Cells (unassigned)_Lake')] <- 'Vascular Endothelium'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('cSMC/MC_Menon','Mesangial_Young','Vascular Smooth Muscle Cells and pericytes_Lake','Pericyte_Wu','Mesangial Cells_Lake')] <- 'Perivascular Cells and Mesangium'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('Fibroblast_ureter_Young','FIB_Menon','Myofibroblast_Wu','Interstitium_Lake','Fibroblast_Wu')] <- 'Fibroblasts'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('PODO_Menon','Glomerular_epithelium_and_podocytes_Young','Podocytes_Lake')] <- 'Glomeruar Epithelium and Podocytes'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('Proximal Tubule Epithelial Cells (S2)_Lake','DTL_Menon','Unknown - Novel PT CFH+ Subpopulation (S2)_Lake','LOH (DL)_Wu','Decending thin Limb_Lake','PEC_Menon','Proximal Tubule Epithelial Cells (S3)_Lake', 'Glomerular parietal eptihelial_Liao')] <- 'Descending Thin Limb, Parietal Epithelium, and Late Proximal Tubule'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('Mono1_Wu','MYL-2_Menon','Immune Cells - Macrophages_Lake','Monocytes_Liao','MNP_Young','MYL-1_Menon','Mono2_Wu','B Cells_Wu')] <- 'Myeloid'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('NK_Young','NK_Menon','NK_proliferating_Young','NKT_Young','NK-T_Liao','T cells_Wu','CD8T_Young','Th_Young','T_Menon')] <- 'Natural Killer and T'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('B_Menon','B_Liao','B_Young','Plasmacytoid_Young','Plasma2_Wu','Plasma1_Wu')] <- 'B, Plasma, and Plasmacytoid'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('Intercalated_Liao','Collecting_duct_type_B_Young','IC_Menon','Collecting Duct - Intercalated Cells Type A (medulla)_Lake','Collecting Duct - Intercalated Cells Type A (cortex)_Lake','IC-B_Menon','Collecting_duct_type_A_Young','Collecting Duct - Intercalated Cells Type B_Lake')] <- 'Intercalated'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('PC_Menon','CD_Wu','Distal_tubule_and_collecting_duct_principle_cells_Young','Principal_Liao','Collecting Duct - Principal Cells (medulla)_Lake','Collecting system - Principal Cells (cortex)_Lake')] <- 'Principal'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('Urothelium_superficial_Young','Urothelium_Young','Urothelium_intermediate_and_basal_Young')] <- 'Urothelium'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('B_Menon','B_Liao','B_Young','Plasmacytoid_Young','Plasma2_Wu','Plasma1_Wu')] <- 'B, Plasma, and Plasmacytoid'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('Thick Ascending Limb_Lake','TAL_Menon','Distal tubule_Liao','Loop_of_Henle_Young','LOH (AL)_Wu','ATL_Menon','Thin ascending limb_Lake')] <- 'Ascending Loop of Henle'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('DCT_Menon','Distal Convoluted Tubule_Lake','Connecting Tubule_Lake','Connecting Tubule_Lake','CNT_Menon')] <- 'Distal Convoluted and Connecting Tubule'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('Mast_Young','Mast cells_Wu')] <- 'Mast'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('Neutrophil_Young')] <- 'Neutrophil'
  
  sobj$Leveled_Names[sobj$Original_Cell_Types %in% c('Proximal_tubule_convoluted_and_straight_Young','Proximal Tubule Epithelial Cells (S1/S2)_Lake','PT_Wu','Proximal_tubule_convoluted_tubular_Young','PT_Menon','Proximal convoluted tubule_Liao','Proximal tubule_Liao','Proximal straight tubule_Liao','Proximal Tubule Epithelial Cells - Fibrinogen+ (S3)_Lake')] <- 'Proximal Tubule'
  
  #print(unique(sobj$Original_Cell_Types[is.na(sobj$Leveled_Names)]))
  return(sobj)
}

sobjm <- rename(sobjm)

sobjmlist <- SplitObject(sobjm, split.by = 'Study')
rm(sobjm)

sobjmlist <- lapply(X = sobjmlist, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = sobjmlist)

sobjmlist <- lapply(X = sobjmlist, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})

anchors <- FindIntegrationAnchors(object.list = sobjmlist, anchor.features = features, reduction = "rpca")
rm(sobjmlist)

sobji <- IntegrateData(anchorset = anchors)
rm(anchors)

sobji <- ScaleData(sobji)

sobji <- RunPCA(sobji)

reticulate::use_python('/opt/venv/bin/python3', required = TRUE)

sceasy::convertFormat(sobji, from = 'seurat', to = 'anndata', outFile = 'data/integratedobject.h5ad')






