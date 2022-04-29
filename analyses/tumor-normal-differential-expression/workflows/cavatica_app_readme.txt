#Introduction

DESeq application takes the gene expression counts data from the pediatric cancer cohorts and normal tissue from GTEx to determine the differential expression genes between all pediatric tumors and normal tissue. This data will be used in the Open Targets platform to help investigators identify potential therapeutic targets in various childhood cancers.



#Method

This app uses R for data manipulation and differential expression analysis using packages DESeq2, jsonlite, among others. 
The data set includes cancers with at least 3 cases, and DESeq analysis was performed on every combination of tumor and normal tissue. Data files for the analysis are available in the [OpenTarget Open Access](http://https://cavatica.sbgenomics.com/u/cavatica/opentarget "OpenTarget Open Access") (now Molecular Targets) project. User can also choose to upload their files in similar format.




#Result

The app delivers analysis output for each cancer cohort and normal tissue pair, with the gene symbol, logfc, pval, EFO ID, MONDO ID, mean TPM values in jsonl, json, and tsv formats.

#Contact
For additional information or questions, contact [Alvin Farrel](mailto:farrela@chop.edu) or [Sangeeta Shukla](mailto:shuklas1@chop.edu).
