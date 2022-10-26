# Molecular Subtyping NBL
**Module authors**: Aditya Lahiri, Eric Wafula, and Jo Lynne Rokita
To molecularly subtype neuroblastoma, ganglioneuroblastoma, and ganglioneuroma samples into MYCN amplified or MYCN non-amplified.

**Note**: The files in the input folder are NBL MCYN `clinical patient-status mapping files` for GMKF and TARGET samples. In V11 some values of the `pathology_free_text_diagnosis` are missing, these input files fill in some of these missing values. 

# Usage
To run the module, execute the following command from the command line in the directory `OpenPedCan-analysis/analyses/molecular-subtyping-NBL` 

`bash run_nbl_subtyping.sh `

# Module Content
`00-Analysis-RMD`: This file contains the main analysis pipeline for subtyping MYCN-NBL.
`run_nbl_subtyping.sh`: Exectues the markdown file.
`results`: Contains the alteration table `Alteration_Table.tsv` and the table `Suggested_Change.tsv` listing the biospecimens where the subtypes were reassigned based on a TPM cut off.

`plots`: Contains the TPM vs Biospecimen plots, and plots of amplification and deletions for samples status call is non-amplified but pathology_free_text_diagnosis is amplified. 