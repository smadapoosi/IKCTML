R Notebook
================

# On Adding New Data to Our Reference

The workflow of our identification model can be adapted to include new
datasets in our reference by the following process:

1.  Obtain well annotated kidney cell dataset in a Seurat Object,
    placing the original cell type labels in a metadata column named
    ‘Original\_Cell\_Types’, and the study the data is from in a
    metadata column named ‘Study’. Ideally, the new data will also have
    information for the other metadata columns (If you are reading this,
    I will assume you either know how to check the metadata of a Seurat
    object and add your own or will understand the process after this
    linked
    [tutorial](https://satijalab.org/seurat/articles/essential_commands.html).

2.  Merge new data with the existing data using a workflow similar to
    the CreationOfMergedObject.ipynb

3.  Plot a heatmap of the correlations of the new annotation labels
    along with the existing ones in a workflow similar to the
    correlational heatmap creation in figures.md. This will show whether
    the new cell labels should be included in an existing category or
    are part of an entirely new one. We assigned our annotations to
    categories based on both the correlations between cell types and the
    descriptions of the cells inside. In our workflow, any added cells
    must undergo this same manual assignment process.

4.  Run IntegrateAll.R and QualityControl.py. These will not be affected
    by new data.

5.  Now, in IntegrateSome.R, there is a variable called ‘iterator’ at
    line 14. The name of the new dataset needs to be added, as this will
    be the name of the train-test split object produced in which the new
    data is the testing dataset.

    Additionally, there is a for loop that currently runs 1:5 at line
    26, this needs to be changed to run for the new appropriate number
    of datasets.

6.  In PerformanceTesting.py, the variable ‘withoutnames’ needs to be
    updated with the appropriate names for testing dataset iterations.

    Next, the original labels from the new dataset need to be assigned
    to a the appropriate Leveled\_Names in the rename() function found
    in lines 34-84.

    Then, update the iterator in lines 115 & 116 to an appropriate
    number to include the new data. You will also need to add a list to
    the nested list of lists in this cell.

    Next, add a list to the nested list of lists in line 175

    Then, add your new dataset Study label to the df.columns in line 193

    Next, add a list to the nested lists in line 223

    Then, update the range at line 242

    Next, add your study to ‘datasetlist’ at line 269

7.  Finally, create the new reference. Take the final seurat object from
    IntegrateAll.R, and remove the cells identified in SVM Outlier
    Detection. An example of this can be found in lines 7-17 of
    IntegrateSome.R

8.  This is now your new reference object, which can be used in place of
    the ref.RDS downloaded from Zenodo in our Colab Annotation workflow!
