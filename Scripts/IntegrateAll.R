library(Seurat)
args <- commandArgs(TRUE)

sobjm <- readRDS(args[1])

sobjm[['Leveled_Names']] <- rep('Unknown', length(sobjm$nCount_RNA))

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

#reticulate::use_python('/opt/venv/bin/python3', required = TRUE)
# ^ This line was pointing to a nonexistant python... I think -- Adam 6/25/22 17:28

sceasy::convertFormat(sobji, from = 'seurat', to = 'anndata', outFile = 'data/integratedobject.h5ad')






