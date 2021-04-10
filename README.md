# Identification of Cell Types in scRNA-seq Data Using Machine Learning Algorithms

## Introduction

According to the CDC, kidney disease is the ninth leading cause of death in the United States, affecting more than 1 in 7 adults. However, advancements in RNA sequencing technologies promise to provide answers, giving revolutionary insight into the complex mechanisms of kidney disease at cell-level resolution. This project seeks to increace the accessablity of data-driven medicine for the human kidney by compiling 5 expert-annotated sn and scRNAseq datasets of healthy, adult human kidney cells and harnessing them to automatically identify unknown query data.


### Data Composition

To maximize the applicability of our model, we sought to include the most diverse a collection of samples we could obtain. We used cells from different single cell and single nucleus sequencing technologies, biopsy locations, ages, and sexes. A summary of our samples is shown in the table below.

![image](https://user-images.githubusercontent.com/77076900/114285173-18784480-9a23-11eb-9800-8bec8db82e00.png)


### Workflow

The workflow of our project is visualized below. After obtaining our data, we removed poorly annotated cells by original author notes, UMAP visualization, and SVM outlier detection. Next, we merged and standardized the samples before removing batch effects using Seurat rPCA integration. Finally, we standardized the ontology across studies by plotting the correlation between original author annotations. This processed data was then used to evaluate the efficacy of various machine learning models by predicting cell types in a single study using the other four as a naive reference.

![image](https://user-images.githubusercontent.com/77076900/114284866-e1a12f00-9a20-11eb-8ef9-3f864777b0c3.png)


### Performance

Our performance testing included 5 different models, shown below, in a rejection scheme that marked uncertian cells as unknown. We used 0.6 as our threshold for rejection; however, this threshold may be tuned for more specific applications.
- Support Vector Machine
- RandomForest
- Multi-Layer Perceptitron
- KNeighbors
- XGBoost

All algorithms showed a strong performance, as visualized in the heatmap below. The best performer was XGBoost, which achieves an average median F1 score of 0.98 across the five datasets along with a average median rejection rate of known cells just slightly above 0. 

![image](https://user-images.githubusercontent.com/77076900/114285109-94be5800-9a22-11eb-83eb-390235c1e9ff.png)

## Data Visualization

Our collected data is interactively visualized as a cellxgene object at this [link](https://nephromappercellxgene.herokuapp.com/) (Note it may require a minute or two to load). This object allows for a detailed exploration of the genes expressed by our cells along their various metadata annotations including cell type, study, sex, age, technology, and original author annotation. 

## To Reproduce Our Results

Start by downloading our dataset MergedObject.RDS from [Zenodo](https://zenodo.org/record/4671060#.YG5Dby1h0YI).

Place the dataset in a directory named data inside the root repository so that you can access it by data/ from the directory containing the Snakefile. (i.e. if root directory was named home, data should be on the path home/data/)

Next, make sure that your system has [singularity](https://sylabs.io/guides/3.0/user-guide/installation.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) installed, and load them into your environment. For our Linux server, we type the commands: 
```
module load singularity 
conda activate snakemake
```
Next, change your directory to the one containing the Snakefile.

Finally, type in this command to run the script: 
```
snakemake --use-singularity --cores <n cores>
```
Here's what happening under the hood:
1. (File 1: IntegrateAll.RDS) Integration of all five datasets
  - Input: MergedObject.RDS Seurat Object pre-QC cells
  - Output: Integrated h5ad object of pre-QC cells
2. (File 2: QualityControl.py) SVM quality control
  - Input: Integrated h5ad object of pre-QC cells from 1
  - Output: a CSV of binary cell designations from QualityControl.py
3. (File 3: IntegrateSome.RDS) Integration of the first 4 datasets and then with the last one separately
* Input: MergedObject.RDS Seurat Object pre-QC cells and CSV of binary cell designations from QualityControl.py
* Output: 5 h5ad objects, each of which will simulate a user query
4. (File 4: PerformanceTesting.py) SVM Performance testing
* Input: 5 h5ad objects from IntegrateSome.RDS
* Output: our figures and classification report (unknown percentages, f1 scores, and confusion matricies)

## To Query Our Reference

First, create a [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) object for your data. 

Second, download our integrated reference from [Zenodo]() and upload both our reference and your data to Google Drive.

Third, open our colab workflow with the link below and follow the included instructions to produce an annotated Seurat object saved to your Google Drive.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/gist/adtisch/d3f445882f32c9139a56e5772d0dd7f7/annotation-workbook.ipynb)

## Questions and Issues
If you have any questions or run into issues please leave your problem in the issues tab or contact us by email.

Maintainers: Stephen Blough <bloughst@umich.edu> & Adam Tisch <adtisch@umich.edu>

## Citations
### Included Datasets:
- [Lake, B.B. et al. A single-nucleus RNA-sequencing pipeline to decipher the molecular anatomy and pathophysiology of human kidneys. Nat Commun 10, 2832 (2019).](https://doi.org/10.1038/s41467-019-10861-2)
- [Liao, J., Yu, Z., Chen, Y. et al. Single-cell RNA sequencing of human kidney. Sci Data 7, 4 (2020).](https://doi.org/10.1038/s41597-019-0351-8)
- [Menon, R. et al. Single cell transcriptomics identifies focal segmental glomerulosclerosis remission endothelial biomarker. JCI Insight 5, e133267 (2020).](https://doi.org/10.1172/jci.insight.133267)
- [Wu, H. et al. Single-cell transcriptomics of a human kidney allograft biopsy specimen defines a diverse inflammatory response. J Am Soc Nephrol 29: 2069–2080 (2018).](https://doi.org/10.1681/asn.2018020125)
- [Young, M. D. et al. Single-cell transcriptomes from human kidneys reveal the cellular identity of renal tumors. Science 361, 594–599 (2018).](https://doi.org/10.1126/science.aat1699)

### CDC Statistics:
- [Centers for Disease Control and Prevention. Chronic Kidney Disease in the United States, 2021. Atlanta, GA: US Department of Health and Human Services, Centers for Disease Control and Prevention; 2021.](https://www.cdc.gov/kidneydisease/publications-resources/ckd-national-facts.html)

### Project Inspiration:
- [Abdelaal, T. et al. A comparison of automatic cell identification methods for single-cell RNA sequencing data. Genome Biol 20, 194 (2019).](https://doi.org/10.1186/s13059-019-1795-z)



