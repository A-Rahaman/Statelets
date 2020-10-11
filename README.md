# Statelets

-----------------
# Project overview
----------------- 
Statelets is a novel time series motifs discovery & summarization framework. We present a probabilistic pattern mining framework for state shape representation of a dynamical system. We applied the model on nueroimaging dataset called "fBIRN" that consists of resting fMRI data from approximately 151 patient with schizophrenia (SZ) and 163 healthy control (HC) subjects. The approach runs group ICA to genertae the independent components of interest followed by sliding window correlation technique to etimate the dynamic functional networks connectivity (dFNC) of the subjects. In our experiments, we use dFNC time courses for motifs discovery and summarization.        

# Major steps.

1. Local motifs discovery from each time series
	
	- apply EMD distance as a similarity measure  
	
2. Summarization of these motifs into k number of statelets describing the most frequent trends in the data.
    
	- incorporate Kernel density estimator (KDE) using a Gaussian kernel for estimating the probability density
	- a 2D peak finding heuristic on tSNE (on a bag of shapes) to determine the dense point in the subspace

-----------------------------------------
# Our implementation of brain imaging data	
-----------------------------------------
The relevant paper for fBIRN dataset's demographics and details is (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4651841/) 
We separated the time series into two different groups SZ and HC. Our analysis follows a group wise & pairwise strategy. 
In other words, we approximate local high-density motifs for a given pair across all the subjects of a group and continue this for all the pairs (1081)  

------------
# Prerequisite 
------------
Install MATLAB version 2019a or later  

------------------------
# Running the experiments 
------------------------ 

1(a). Adjust the directories and dependencies describes within a module. Select a group SZ/HC. For a given pair 'pr' of that group, run the steps and "step1_getlocalmotifs\getmotifsPairWise.m" discovers the most recurring patterns in each time series and store corresponding information 

1(b). Run four steps methodology to create pairwise instances of shapes. Peform EMD-> tSNE -> KDE -> peak finding on the bag of local motifs extracted from all the signal 
that represents pair 'pr' in the group. So, the function generates pairwise a subset of high-density motifs for each pair


2(a). Select a group SZ/HC, collect local motifs generated in 1(c) from all the pairs using subroutine 'step4_finding dominants\Statelets_from_SZandHC.m'

2(b). Again, perform EMD-> tSNE -> KDE -> peak finding to generate Statelets for a given dynamic 

-----------
# Data usage
-----------
The fBIRN dataset was collected 10 years ago, and no data-sharing agreement was initiated at that time. According to the IRB data privacy agreement,
we are not allowed to share any subject-specific data. However, the resting-state fMRI data is available on request from fBIRN project site, 
and preprocessing pipelines is public too, as I referred to in the manuscript. So, please reach me out through mrahaman8@gatech.edu if there is anything required and missing

---------------------------------------
# General tips for running the framework
---------------------------------------
Our implementation requires storing MATLAB data from various steps and load them back to the workplace in some later steps. Therefore, it is 
necessary to set the appropriate directories before running the method. 
We ran our project on high-performance computing nodes using 'SLURM' and other workload managers.
So, the scripts are parallelized and organized in that way            
