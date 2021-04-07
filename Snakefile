

"""
One rule to... rule... them all...
"""
rule all:
  input:
      "confusion/ukdf.csv"
      
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
    conf = "confusion/ukdf.csv"
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
