# Identification of Cell Types in scRNA-seq Data Using Machine Learning Algorithms

## Introduction

According to the CDC, kidney disease is the ninth leading cause of death in the United States, affecting more than 1 in 7 adults. However, advancements in RNA sequencing technologies promise to provide answers, giving revolutionary insight into the complex mechanisms of kidney disease at cell-level resolution. This project seeks to increace the accessablity of data-driven medicine for the human kidney by compiling 5 expert-annotated sn and scRNAseq datasets of healthy, adult human kidney cells and harnessing them to automatically identify unknown query data.


### Training Data Composition

To maximize the applicability of our model, we sought to include the most diverse a collection of samples we could obtain. We used cells from different single cell and single nucleus sequencing technologies, biopsy locations, ages, and sexes. A summary of our samples is shown in the table below.

| Study | Number of Cells | Number of Donors | Number of Samples | Sampling Locations | Age | Sex | Technology |
| ------|-----------------|------------------|-------------------|---------------------|-----|----|---------------------|
| Lake et al. | 13,255 | 10 | 35 | 15 Cortex; 14 Medulla; 6 Both | 12 Samples < 50; 23 Samples > 50 | 10 M; 11 F| sn Drop-Seq |
 | Liao et al. | 16,145 | 2 | 2 | Unknown | 59, 65 | 1 M; 1 F | sc 10X |
 | Menon et al. | 22,264 | 22 | 24 | Unknown | 2 Samples < 50; 13 Samples > 50 | 7 M; 6 F | sc 10X |
 | Wu et al. | 4,432 | 1  | 4 | Unknown | 21 | M | sn InDrops |
 | Young et al.  | 6,179  | 5  | 17  | 14 Cortex; 1 Both; 2 Ureter | 49; 63; 67; 70; 72 | 3 M; 2 F | sc 10X |


### Workflow

The workflow of our project is visualized below. After obtaining our data, we removed poorly annotated cells by original author notes and our own SVM outlier detection protocol. Next, we merged and standardized the samples before removing batch effects using Seurat rPCA integration. Finally, we standardized the ontology across studies by plotting the correlation between original author annotations. This processed data was then used to evaluate the efficacy of various machine learning models by predicting cell types in a single study using the other four as a naive reference.

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


## To Reproduce Our Results

Start by downloadng our dataset MergedObject.RDS from [Zenodo](https://zenodo.org/record/4671060#.YG5Dby1h0YI).

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
1. (File 1 IntegrateAll) Integration of all five datasets
  - Input: Merged RDS Seurat Object pre-QC cells
  - Output: Integrated h5ad object of pre-QC cells
2. (File 2 QualityControl) SVM quality control
* Input: Integrated h5ad object of pre-QC cells from 1
* Output: a CSV of binary cell designations from Quality control
3. (File 3 IntegrateSome) Integration of 4 datasets and then with the last one separately
* Input: Merged RDS Seurat Object pre-QC cells and CSV of binary cell designations from Quality control
* Output: 5 h5ad objects, each of which is is first integrated with four of the datsets and then again with last one, which will simulate a query dataset
4. (File 4 PerformanceTesting) SVM Performance testing
* Input: 5 h5ad objects from 3
* Output: our figures and classification report (unknown percentages, f1 scores, and confusion matricies)

## To Query Our Reference

The best performing annotation protocol from our performance testing showed an average median F1 score 

First, create a [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) object for your data. 

Second, download our integrated reference from [Zenodo]() and upload both our reference and your data to Google Drive.

Third, open our colab workflow with the link below and follow the included instructions to produce an annotated Seurat object saved to your Google Drive.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/gist/adtisch/d3f445882f32c9139a56e5772d0dd7f7/annotation-workbook.ipynb)

## Questions and Issues
If you have any questions or run into issues please leave your problem in the issues tab or contact the maintainers.

Maintainers: Stephen Blough <bloughst@umich.edu> & Adam Tisch <adtisch@umich.edu>

### Citations
Our Merged Kidney Dataset used a combination of datasets and studies from:
- Menon et al.
- Lake et al.
- Wu et al.
- Liao et al.
- Young et al.


