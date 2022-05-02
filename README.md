# Identifying biomarkers for ischemic stroke severity using mouse striatum transcriptomic profiles
 
## Table of Contents
1. **Overview**
    * Introduction
    * Raw data description
2. **Installation**
    * Folder structure
    * Installing dependencies
3. **Generating Analysis**
    * _main.py_

## I. Overview
### Introduction
This repo servers two main purposes:
1. To perform principal component analysis to compare mouse stroke transcriptomic profiles in the striatum and cortex
2. To perform one-way ANOVA test with Benjamini-Hochberge multiple hypothesis correction to identify biomarkers indicating stroke severity
These analyses will be performed on raw microarray transcriptomic profiling of ischemic/reperfusion injury in an in vivo wild-type mouse model.

 
### Raw data description

#### Microarray stroke/ reperfusion transcriptomic profile 
We used mouse ischemic stroke expression profiles (GEO database) for the cortex and striatum to identify correlatory biomarkers associated stroke severity. These files can be found in the _data/I_R_txt_ directory. 
Each contain normalized gene expression counts of 16 striatum and cortex samples of 4 different experimental groups (controls, 2 h after stroke, 8 h after stroke, and 24 h after stroke). The expressions of 25,967 genes of each sample are stored in one separate \*.txt file.

## II. Installation
### Folder structure
This repo is divided into three directories: _src_, _data_, and _results_.
- The _data_ directory is read-only, containing transcriptomic data, microarray probe information, and the Human Protein Atlas dataset.
- The _results_ directory contains output figures and \*.csv files.
- The _src_ directory contains scripts used to perform our analyses. The script _process_data.py_ contains the functions used, and the script _main.py_ generates the files that will be stored in the _results_ directory.

### Installing dependencies
To perform these analyses, we need to install core depencies. 
First, we will need to install the latest version of Python, GCC, and the following packages: `numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`, `statsmodel`, and `scikit-learn`.
You can use a package manager like Homebrew.


1. Core system dependencies
    ```
    xcode-select --install  
    brew install python3  
    brew install gcc  
    ```
2. Packages  

The latest version of the packages mentioned above can be installed with the following command lines:
   ```
   pip3 install -U numpy scipy pandas  
   pip3 install -U matplotlib seaborn   
   pip3 install -U statsmodels scikit-learn
   ```
   

## III. Generating analysis
Download the striatumStrokeBiomarker repository. Inside the _src_ directory, you will find the _process_data.py_ and _main.py_ scripts. The former contains functions that the latter uses. The _main.py_ script takes the inputs described in the "Raw Data Description" section above, and if the dependencies are successfully installed, this script will create all the figures and \*.csv files included in our report.

Use the commands below:  

    git clone https://github.com/vchauu/striatumStrokeBiomarker    
    cd striatumStrokeBiomarker   
    python -m src.main     
    

## Citations
Chen MJ, Wong CH, Peng ZF, Manikandan J et al. A global transcriptomic view of the multifaceted role of glutathione peroxidase-1 in cerebral ischemic-reperfusion injury. Free Radic Biol Med 2011 Mar 15;50(6):736-48.
