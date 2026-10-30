# Antarctic eDNA Analysis
Western Antarctic Peninsula marine sediment core sedaDNA metabarcoding (18S/12S/16S) biodiversity analysis

This repository accompanies the manuscript **Ancient environmental DNA reveals the distinct fingerprints of whaling and 
climate change on Antarctic marine biodiversity**. Descriptions of scripts and folders can be found below. 

File Name | Description
--- | ---
0.1.DataClean.R | Script used to clean data using negative controls, minimum read and replicate filters, and to count reads per sample
0.2.Mapping.R | Script used to prepare maps of core sites
0.3.Climate.R | Script used to prepare climate data (SAM, ENSO, sea ice concentration) for analysis
0.4.WhaleBiomass.R | Script used to extract and compare whale biomass reconstructions with metabarcoding data
0.5.CheckFamiliesRAMS.R | Script used to screen metazoan families against RAMS and build the functional groups file
1.0.Analysis.R | Main analysis script for figures and ecological analysis of data
archive | Archived code and data
cleandata | Cleaned datasets, negative control results and MCMC chains
figures | Figure outputs
mapping | Map outputs
metadata | Metadata files for samples, sites and core dates
rawdata | Raw datasets after bioinformatic analysis, climate data, biomass data and read counts
summaryfiles | Statistical model outputs and supplementary tables
taxonomy | Taxonomy files, databases and outputs
