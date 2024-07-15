# MonkeyAge
Epigenetic Age Estimators ECRM and ECYRM for Rhesus Macaques

Files for "A novel epigenetic clock for rhesus macaques unveils 
an association between early life adversity and epigenetic age acceleration"
by Bronk et al. 2024.


To use these scripts, first download the AgePrediction folder (including all of its files and subfolders).


The script AgePredictionScript.Rmd computes the epigenetic ages of rhesus macaques after you input a MethyLumiSet object containing their DNA methylation.

If unfamiliar with how to create a MethyLumiSet object, use MethylationProcessingScript.Rmd to turn idat files into a MethyLumiSet object. MethylationProcessingScript.Rmd also removes unwanted CpG probes and unwanted blood samples.

This repository also contains other files that need to be loaded by AgePredictionScript.Rmd or MethylationProcessingScript.Rmd.
