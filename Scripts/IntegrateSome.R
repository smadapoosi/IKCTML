library(Seurat)
args <- commandArgs(TRUE)
sobjm <- readRDS(args[1])

#This CSV is produced by the python notebook QualityControl.ipynb

cells <- read.csv(args[2])
cellsuse <- cells[[2]]
names(cellsuse) <- cells[[1]]

sobjm[['SVMQC']] <- cellsuse
sobjm <- subset(sobjm, subset = SVMQC == unique(sobjm$SVMQC)[1])

iterator <- c('Without Lake', 'Without Liao', 'Without Menon', 'Without Wu', 'Without Young')

sobjmlist <- SplitObject(sobjm, split.by = 'Study')
rm(sobjm)

sobjmlist <- lapply(X = sobjmlist, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

reticulate::use_python('/opt/venv/bin/python3', required = TRUE)

for(i in 1:5){
  sublist <- sobjmlist[-i]
  
  single <- sobjmlist[i]
  
  
  
  features <- SelectIntegrationFeatures(object.list = sublist)
  
  sublist <- lapply(X = sublist, FUN = function(x) {
    x <- ScaleData(x, features = features)
    x <- RunPCA(x, features = features)
  })
  
  print('The Testing Dataset is:')
  print(sobjmlist[i])
  
  print('Integrating Reference')
  anchors <- FindIntegrationAnchors(object.list = sublist, anchor.features = features, reduction = "rpca")
  rm(sublist)
  
  sobji <- IntegrateData(anchorset = anchors)
  rm(anchors)
  
  sobji <- ScaleData(sobji)
  
  sobji <- RunPCA(sobji)
  
  single <- single[[1]]
  single[['Study']] <- rep('Single', length(single$nCount_RNA))
  sublist <- list(single, sobji)
  
  features <- SelectIntegrationFeatures(object.list = sublist)
  
  sublist <- lapply(X = sublist, FUN = function(x) {
    x <- ScaleData(x, features = features)
    x <- RunPCA(x, features = features)
  })
  
  anchors <- FindIntegrationAnchors(object.list = sublist, anchor.features = features, reduction = "rpca")
  rm(sublist)
  
  sobji <- IntegrateData(anchorset = anchors)
  rm(anchors)
  
  sobji <- ScaleData(sobji)
  
  sobji <- RunPCA(sobji)
  
  sobji <- RunUMAP(sobji, dims = 1:30)
  
  sceasy::convertFormat(sobji, from = 'seurat', to = 'anndata', outFile = file.path("data", paste("svmtest",i, ".h5ad", sep="")))
  rm(sobji)
}