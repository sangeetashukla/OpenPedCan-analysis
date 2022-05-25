# EFO and MONDO mapping

## Introduction and Purpose

The Experimental Factor Ontology (EFO) is a resource that provides systematic description of experimental variables including cancer groups available in [European Bioinformatics Institute](https://www.ebi.ac.uk/) databases and some other external projects. 
The scope of EFO is to support the annotation, analysis and visualization of data handled by many groups at the EBI and as the core ontology for [OpenTargets](https://www.opentargets.org/).

The Mondo Disease Ontology (Mondo) is another independent resource aiming to harmonize disease definitions. 

The purpose of maintaining `efo-mondo-map.tsv` is to map each cancer group found in the `histologies.tsv` to its appropriate EFO and MONDO codes.

## Format

The EFO codes may have prefix “EFO_”, “MONDO_” or “Orphanet_”, while the MONDO codes have a prefix “MONDO_”; each followed by seven digits.

## File QC
To make sure all the `cancer_group` with n>=3 in the histologies file are included in the `efo-mondo-map.tsv` file, a QC Rmd script is included in the repo. 
This notebook can be run with the command `bash run_qc.sh`
This QC step should be done for all new histologies file release before generating SNV, CNV and fusion frequencies tables.

## Note
As more cancer groups or subtypes are added to existing histologies dataset and/or in case of any ambiguity, these codes will be revised.
