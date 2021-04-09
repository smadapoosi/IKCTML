rule all:
  input:
      "confusion/ukdf.csv",
      "confusion/SVC/Lake_SVC_confusion.csv",
      "confusion/SVC/Liao_SVC_confusion.csv",
      "confusion/SVC/Menon_SVC_confusion.csv",
      "confusion/SVC/Young_SVC_confusion.csv",
      "confusion/SVC/Wu_SVC_confusion.csv",
      "confusion/RandomForest/Lake_RandomForest_confusion.csv",
      "confusion/RandomForest/Liao_RandomForest_confusion.csv",
      "confusion/RandomForest/Menon_RandomForest_confusion.csv",
      "confusion/RandomForest/Young_RandomForest_confusion.csv",
      "confusion/RandomForest/Wu_RandomForest_confusion.csv",
      "confusion/MLP/Lake_MLP_confusion.csv",
      "confusion/MLP/Liao_MLP_confusion.csv",
      "confusion/MLP/Menon_MLP_confusion.csv",
      "confusion/MLP/Young_MLP_confusion.csv",
      "confusion/MLP/Wu_MLP_confusion.csv",
      "confusion/KNeighbor/Lake_KNeighbor_confusion.csv",
      "confusion/KNeighbor/Liao_KNeighbor_confusion.csv",
      "confusion/KNeighbor/Menon_KNeighbor_confusion.csv",
      "confusion/KNeighbor/Young_KNeighbor_confusion.csv",
      "confusion/KNeighbor/Wu_KNeighbor_confusion.csv",
      "confusion/XGB/Lake_XGB_confusion.csv",
      "confusion/XGB/Liao_XGB_confusion.csv",
      "confusion/XGB/Menon_XGB_confusion.csv",
      "confusion/XGB/Young_XGB_confusion.csv",
      "confusion/XGB/Wu_XGB_confusion.csv",
      "scores/Median_F1_heatmap.pdf",
      "confusion/Percent_Unknown_heatmap.pdf",
      "scores/F1_scores.csv"

rule IntegrateAll:
    input: 
      data = "data/MergedObject.RDS"
    output:
      int = "data/integratedobject.h5ad"
    log: "data/IntegrateAll.log"
    singularity: "docker://sblough18/seurat-integrate:latest"
    shell:
      "Rscript Scripts/IntegrateAll.R "
      "{input.data} "
      "&> {log} "

rule QualityControl:
    input: 
      data = "data/integratedobject.h5ad"
    output:
      tfval = "data/cellstouse.csv",
      unkn = "confusion/UnknownpercentsvmQCconfusion.csv"
    log: "data/QC.log"
    singularity: "docker://sblough18/svm:latest"
    shell:
      "python3 Scripts/QualityControl.py "
      "{input.data} "
      "&> {log} "    

rule IntegrateSome:
    input: 
      data = "data/MergedObject.RDS",
      tfval = "data/cellstouse.csv"
    output:
      int1 = "data/svmtest1.h5ad",
      int2 = "data/svmtest2.h5ad",
      int3 = "data/svmtest3.h5ad",
      int4 = "data/svmtest4.h5ad",
      int5 = "data/svmtest5.h5ad"
    log: "data/IntegrateSome.log"
    singularity: "docker://sblough18/seurat-integrate:latest"
    shell:
      "Rscript Scripts/IntegrateSome.R "
      "{input.data} "
      "{input.tfval} "
      "&> {log} "


rule PerformanceTesting:
  input:
    int1 = "data/svmtest1.h5ad",
    int2 = "data/svmtest2.h5ad",
    int3 = "data/svmtest3.h5ad",
    int4 = "data/svmtest4.h5ad",
    int5 = "data/svmtest5.h5ad"
  output:
    ukdf = "confusion/ukdf.csv",
    data1 = "confusion/SVC/Lake_SVC_confusion.csv",
    data2 = "confusion/SVC/Liao_SVC_confusion.csv",
    data3 = "confusion/SVC/Menon_SVC_confusion.csv",
    data4 = "confusion/SVC/Young_SVC_confusion.csv",
    data5 = "confusion/SVC/Wu_SVC_confusion.csv",
    data6 = "confusion/RandomForest/Lake_RandomForest_confusion.csv",
    data7 = "confusion/RandomForest/Liao_RandomForest_confusion.csv",
    data8 = "confusion/RandomForest/Menon_RandomForest_confusion.csv",
    data9 = "confusion/RandomForest/Young_RandomForest_confusion.csv",
    data10 = "confusion/RandomForest/Wu_RandomForest_confusion.csv",
    data11 = "confusion/MLP/Lake_MLP_confusion.csv",
    data12 = "confusion/MLP/Liao_MLP_confusion.csv",
    data13 = "confusion/MLP/Menon_MLP_confusion.csv",
    data14 = "confusion/MLP/Young_MLP_confusion.csv",
    data15 = "confusion/MLP/Wu_MLP_confusion.csv",
    data16 = "confusion/KNeighbor/Lake_KNeighbor_confusion.csv",
    data17 = "confusion/KNeighbor/Liao_KNeighbor_confusion.csv",
    data18 = "confusion/KNeighbor/Menon_KNeighbor_confusion.csv",
    data19 = "confusion/KNeighbor/Young_KNeighbor_confusion.csv",
    data20 = "confusion/KNeighbor/Wu_KNeighbor_confusion.csv",
    data21 = "confusion/XGB/Lake_XGB_confusion.csv",
    data22 = "confusion/XGB/Liao_XGB_confusion.csv",
    data23 = "confusion/XGB/Menon_XGB_confusion.csv",
    data24 = "confusion/XGB/Young_XGB_confusion.csv",
    data25 = "confusion/XGB/Wu_XGB_confusion.csv",
    medf1 = "scores/Median_F1_heatmap.pdf",
    unkn = "confusion/Percent_Unknown_heatmap.pdf",
    f1 = "scores/F1_scores.csv"
  log: "data/PerformanceTesting.log"
  singularity: "docker://sblough18/svm:latest"
  shell:
    "python3 Scripts/PerformanceTesting.py "
    "{input.int1} "
    "{input.int2} "
    "{input.int3} "
    "{input.int4} "
    "{input.int5} "
    "&> {log}"
