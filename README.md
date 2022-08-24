# BIRCH Usage and Instructions

## Table of Contents
- [What is BIRCH](#what-is-birch)
- [Requirements](#requirements)
- [Example usage](#example-usage)
- [Using BIRCH](#using-pine-gui)
- [Cite us](#cite-us)
- [Support](#support)
- [Release notes](#release-notes)
- [License](#license)

## What is BIRCH
[BIRCH](https://batch-correction-tool.herokuapp.com/) (**B**atch-effect **I**dentification, **R**ectification, and **C**onclusive-anlaysis on **H**eterogenous data) is a web-app that can be used for reducing batch-effect in proteomics data. It generally becomes necessary to correct for batch-effect in large datasets since processing steps such as sample preparation and data acquisition tends to add noise to the data, that in-turn effects biological conclusions. This tool aims are keeping meaningful biological variation, while simultaneously reducing batch-effect due to other external factors. 

## Requirements
- Web Browser (Google Chrome, Mozilla Firefox, Microsoft Edge, Safari) with stable internet connection. 
- Input files with protein intensitites and an annotations (explained in detail in the next section). 

## Example usage
Please refer to our example (XXX link it later) for a walkthrough on how to use BIRCH.

## Using BIRCH

### Settings and Input files
When the BIRCH web-app first loads, the settings tab is seen.

![BIRCH settings page](Images/settings.PNG)

In this tab, you will have to upload 3 files:
1. Normalized data 
2. Unnormalized data 
3. Annotation file

Additionally, you will be provided with dropdown menus to select the following mandatory parameters:
1. Column representing the protein name in the normalized and unnormalized files.
2. Columns to correct for batch-effect (can choose multiple based on your experiment set-up).
3. Column containing the experimental/biological group (this is the group in which we want to retain variation coming from meaningful biological differences).
4. Column with sample names that match the headers in the normalized and unnormalized file. 

### Initial analysis tab

### Diagnosis tab

### Results tab