# Identification of Cell Types in scRNA-seq Data Using Machine Learning Algorithms

## Introduction
According to the CDC, kidney disease is the ninth leading cause of death in the US, affects more than 1 in 7 US adults, and costs Medicare about $120 B per year. This project seeks to increace the accessability of data-driven precision kidney medicine by creating a reference map of more than 62,000 kidney cells from 45 donors in 5 published studies and using it to annotate user query data. 

Our reference map is visualized as an interactive cellxgene application, linked here.

The code in this repository will both recreate our dataset and performance evaluation as well as implement our reference for user query data.

Our performance testing workflow can be replicated with our Snakefile, which intakes the merged, unintegrated datasets from each of our 5 studies and automatically reproduces our workflow, results, and figures. 

### Included Algorithms
We tested 5 algorithms for cell type classification, each following a rejection model where cells with low classification probabilities are marked as unknown:
- SVM
- RandomForest
- MLP
- KNeighbors
- XGB

## How to Use

Start by downloadng our dataset MergedObject.RDS from [zenodo](https://zenodo.org/record/4671060#.YG5Dby1h0YI).

Place the dataset in a directory named data inside the root repository so that you can access it by data/ from the directory containing the Snakefile. (i.e. if root directory was named home, data should be on the path home/data/)

Next, make sure that your system has [singularity](https://sylabs.io/guides/3.0/user-guide/installation.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) installed, and load them into your environment. For our server, we type the commands: 
```
module load singularity 
conda activate snakemake
```
Next, change your directory to the one containing the Snakefile.

Finally, type in this command to run the script: 
```
snakemake --use-singularity --cores <n cores>
```
## Results
Here's what to expect when process finishes:
1. (File 1 IntegrateAll) Integration of all five datasets
       * Input: Merged Object pre-QC cells
       * Output: Integrated h5ad object of pre-QC cells, quality control report

2. (File 2 QualityControl) SVM quality control
       * Input is Integrated object of 62000 cells
       * Output is a CSV of True or False from Quality control
3. (File 3 IntegrateSome) Integration of 4 datasets and then with the last one separately
        1. Input: Merged Object, and QC csv
        2. Output: 5 h5ad objects, each of which is different
4. (File 4 PerformanceTesting) SVM Performance testing
        1. Input is 5 h5ad objects 
        2. Output is the Figures (unknown percent and f1 scores, confusion matrix, cell by cell unknown percent)

There are also various other csvs and objects created in the script, but they are less important

## Questions and Issues
If you have any questions or run into issues please leave write your problem in the issues tab or contact the maintainers.

Maintainers: Stephen Blough <bloughst@umich.edu> & Adam Tisch <adtisch@umich.edu>

### Citations
Our Merged Kidney Dataset used a combination of datasets and studies from:
- Menon et al.
- Lake et al.
- Wu et al.
- Liao et al.
- Young et al.


