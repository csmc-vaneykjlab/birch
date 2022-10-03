# Running BIRCH using example dataset

## Table of Contents
- [Home](#home)
- [Settings](#settings)
- [Initial Analysis](#initial-analysis)
- [Diagnosis+Filtering](#diagnosis+filtering)
- [Results](#results)

### Home
When the BIRCH web-app initially loads, the home tab is displayed at first. The home tab just contains details on what BIRCH is and the workflow followed in BIRCH to go from data affected by batch-effect to batch-corrected data. 

![BIRCH settings page](../Images/home.PNG)

### Settings
From the home tab, you can navigate to the settings tab to set input parameters for batch correction. For this example, we will use a dummy dataset to demonstrate the entire process and make note of the things to look out for while dealing with your own data.  

To begin with analysis, you can either download [example_dataset.zip](https://github.com/csmc-vaneykjlab/BatchCorrectionTool/tree/main/example/example_dataset.zip) the example dataset and give it as input to the respective parameters, or you can click on "Load example data" to pre-load all the input parameters. 

Using the example dataset, once the files are uploaded, the other parameters that need to be chosen are as follows:
Argument 3 - ProteinName - this is the column that contains unique fragment/peptide/protein names.  
Argument 5 - Digestion_batch - this is the column from the annotation file that represents the field/batch that needs to be analyzed for batch-effect and corrected. 
Argument 6 - ExperimentalGroup - this column from the annotation file that represents the biological group in which we want to see and retain valid variation. 
Argument 7 - SampleName - this argument maps the sample names (columns) in normalized/unnormalized intensities file, with the sample names in the annotation file where each row is a sample. 
