# Identification of Cell Types in scRNA-seq Data Using Machine Learning Algorithms

## Introduction
Our Scripts and Snakefile will create important files including the percent unknown and and median f1 scores for every dataset used in one of our five chosen machine learning  algorithms. Also included are confusion matrices, and a cell-by-cell breakdown of the percent unknown for each training dataset and algorithm used.

### Included Algorithms
Most of the Algorithms we opted for are based in the python sklearn:
- SVM
- RandomForest
- MLP
- KNeighbors
- XGB

## How to Use

Start by downloadng our dataset MergedObject.RDS from [zenodo](https://zenodo.org/record/4671060#.YG5Dby1h0YI).

Make sure the dataset is placed in a directory named data inside the root repository so that you can access it by data/ from the directory containing the Snakefile. (i.e. if root directory was named home, data should be on the path home/data/)

Next make sure that your system has singularity and snakemake installed on it, and load them into your environment. For our server, we type the commands: 
```
module load singularity 
conda activate snakemake
```
Also change your directory to the one containing the Snakefile.

Finally to type in the command: 
```
snakemake --use-singularity --cores <n cores>
```
## Results
Here's what to expect when process finishes:
1. (File 1 IntegrateAll) Integration of all five datasets
        a. Input: Merged Object pre-QC cells
        b. Output: Integrated h5ad object of pre-QC cells, quality control report
2. (File 2 QualityControl) SVM quality control
        a. Input is Integrated object of 62000 cells
        b. Output is a CSV of True or False from Quality control
3. (File 3 IntegrateSome) Integration of 4 datasets and then with the last one separately
        a. Input: Merged Object, and QC csv
        b. Output: 5 h5ad objects, each of which is different
4. (File 4 PerformanceTesting) SVM Performance testing
        a. Input is 5 h5ad objects 
        b. Output is the Figures (unknown percent and f1 scores, confusion matrix, cell by cell unknown percent)

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


