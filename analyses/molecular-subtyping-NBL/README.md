# Molecular Subtyping NBL
**Module authors**: Aditya Lahiri (@adilahiri), Eric Wafula (@ewafula), and Jo Lynne Rokita (@jharenza)

To molecularly subtype neuroblastoma, ganglioneuroblastoma, and ganglioneuroma samples into MYCN amplified or MYCN non-amplified.

**Subtyping criteria**:

**case 1**:
If `pathology_free_text_diagnosis` is **amplified** and `status` is **amplified** assign subtype as `NBL, MYCN amplified`

**case 2**:
If `pathology_free_text_diagnosis` is **non-amp** and `status` is **amplified** assign subtype as `NBL, MYCN amplified`

**case 3**:
If `pathology_free_text_diagnosis` is **non-amp** and `status` is **non-amp** assign subtype as `NBL, MYCN non-amplified`

**case 4**:
If `pathology_free_text_diagnosis` is **amplified** and status is **non-amp**
For such samples if there exists a TPM value, evaluate if the TPM is above or below the Suggested_Cutoff established in the image `results/TPM_Biospecimen_All_Samples_With_TMP.png` assign subtype as `NBL, MYCN amplified` or `NBL, MYCN non-amplified` respectively.  In case there is no TPM values then assign the subtype as `NA`. 

**Note**: The files in the input folder are NBL MCYN `clinical patient-status mapping files` for GMKF and TARGET samples. In V11 some values of the `pathology_free_text_diagnosis` are missing, these input files fill in some of these missing values. 

# Usage
To run the module, execute the following command from the command line in the directory `OpenPedCan-analysis/analyses/molecular-subtyping-NBL` 

`bash run_nbl_subtyping.sh `

# Module Content
`00-Analysis-RMD`: This file contains the main analysis pipeline for subtyping MYCN-NBL. 

`run_nbl_subtyping.sh`: Bash script to execute `00-Analysis-RMD`.

`input/gmkf_patient_clinical_mycn_status.tsv`: Has values for `pathology_free_text_diagnosis` in the GMKF cohort. These values are not present in the V11 data. 

`input/target_patient_clinical_mycn_status.tsv`: Has values for `pathology_free_text_diagnosis` in the TARGET cohort. These values are not present in the V11 data. 

`results`: Contains the alteration table `Alteration_Table.tsv` and the table `Suggested_Change.tsv` listing the biospecimens where the subtypes were reassigned based on a TPM cut off.

`plots`: Contains the TPM vs Biospecimen plots, and plots of amplification and deletions for samples status call is non-amplified but pathology_free_text_diagnosis is amplified. 










The data files being used in this analysis are:
```
histologies-base.tsv

consensus_wgs_plus_cnvkit_wxs.tsv.gz

cnv-cnvkit.seg.gz

cnv-controlfreec.tsv.gz

gene-expression-rsem-tpm-collapsed.rds
```
The module extracts the NBL samples from the `histologies-base.tsv` file by filtering for the following fields:

`pathology_diagnosis`
```
"Neuroblastoma"

"Ganglioneuroblastoma"

"Ganglioneuroblastoma, nodular"

"Ganglioneuroblastoma, intermixed"

"Ganglioneuroma, maturing subtype OR Ganglioneuroblastoma, well differentiated"

```




`sample`

