# Identifying risk biomarkers for future ischemic stroke using scRNA-seq mouse brain vascular transcriptomes
 
## Table of Contents
1. **Overview**
    * Introduction
    * Raw data description
2. **Installation**
    * Folder structure
    * Installing dependencies
3. **Partial least square-dimensional analysis**
    * Protocol
    * Data processing
   

## I. Overview
### Introduction
The purpose of this repo is to perform partial least square-dimensional analysis (PLS-DA) for mouse brain vascular transciptomes
to identify risk biomarkers for future ischemic stroke. 
Using raw scRNA-seq data for mouse brain transcriptomes and microarray transcriptomic profiling of ischemic/reperfusion injury 
in an in vivo wild-type mouse model.

 
### Raw data description

#### scRNA mouse brain vascular and microarray stroke/ reperfusion transcriptomes
We will draw a comparative analysis between mouse vascular transcriptomes (EBI database) and mouse ischemic stroke expression profiles (GEO database) 
to identify a correlatory biomarker associated with risk of vascular ischemic events. These files can be found in the _data_ directory. 
The normalized gene expression counts of 3186 samples are stored as sparse matrices with 3186 columns and 28338 rows, which comprise of a matrix (\*.mtx), matrix columns (\*.mtx_col), and matrix row (\*.mtx_row) file.
Files containing information about individual data points such as cell type, group (control vs. diseased) are stored as text files (\*.txt).

## II. Installation
### Folder structure
This repo is divided into three directories: _src_, _data_, and _fig_.
The _data_ directory is read-only, containing transcriptomic data.
The _results_ directory contains output figures
The _src_ directory contains scripts we used to perform our analysis. The script _main.py_ will generate one figure.

### Installing dependencies
To perform these analyses, we need to install core depencies. 
First, we will need to install the latest version of Python, GCC, and the following packages: `numpy`, `scipy`, `pandas`, and `matplotlib`.
You can use a package manager like Homebrew.
Then, we will install the `scprep` to perform RNA-seq analysis. Installation instructions are documented here:
https://scprep.readthedocs.io/en/stable/installation.html
Installation instructions are summarized below:

1. Core system dependencies
    ```
    xcode-select --install  
    brew install python3  
    brew install gcc  
    ```
2. scprep  
   Using `pip`:
   ```
   pip install --user scprep
   ```
   or from source code:
   ```
   git clone git://github.com/KrishnaswamyLab/scprep.git
   cd scprep/python
   python setup.py install --user
   ```

## III. PLS-DA
### Testing the compatibility of `scprep ` with the brain vascular transcriptome RNA-seq data (for exploratory purposes only)
We will first draw a simple gene expression histogram to gauge whether we can load the data and whether `scprep` is compatible with the data format. 
The plot delineates the variation in expression of all the genes sequenced in the 3168 samples and is saved in the _results_ directory.

### Performing PLS-DA analysis


## Citations
Dagonnier M, Donnan GA, Davis SM,Dewey HM and Howells DW (2021)Acute Stroke Biomarkers: Are WeThere Yet? Front. Neurol. 12:619721