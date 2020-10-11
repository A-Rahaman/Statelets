%{
Instructions:
Local discovery:

- Run the step 1 script 'getmotifsPairwise' to perform the local motifs discovery for allt eh time series
- It will return the high density motifs for each pair
- Other information as well. Read the instruction in aformentioned script to know what are the 
- information the subroutine provides.
    
summarization: Groupwise Statelets  

- Using 'Statelet_from_SZandHC.m' under step4 to collect all these pairwise peak motifs and create a bag of shapes for SZ and HC
- Compute EMD -> probability density (KDE)-> tSNE for both SZ and HC (script in step2) 
- Findout the peak motifs with higher probability density 
- Finally create two subsets of statelets from SZ and HC 
--------------Use these statelets for any further analysis,  --------------------------------- 

%}