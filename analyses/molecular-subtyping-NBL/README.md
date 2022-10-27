# Molecular Subtyping NBL
**Module authors**: Aditya Lahiri ([@adilahiri](https://github.com/adilahiri)), Eric Wafula ([@ewafula](https://github.com/ewafula)), and Jo Lynne Rokita ([@jharenza](https://github.com/jharenza))

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

`input/gmkf_patient_clinical_mycn_status.tsv`: Has values for `pathology_free_text_diagnosis` in the GMKF cohort. These values are not present in the V11 data and will be made available in V12.

`input/target_patient_clinical_mycn_status.tsv`: Has values for `pathology_free_text_diagnosis` in the TARGET cohort. These values are not present in the V11 data and will be made available in V12.

`results/NBL_MYCN_Subtype.tsv`: This table is the main output of this module. For each biospecimen sample, the table  lists its DNA and RNA IDs (if available) and the subtype. 

`results/Alteration_Table.tsv`: This table is similar to the table `NBL_MYCN_Subtype.tsv`. However, this table has 
additional columns which contain information on `MYCN_TPM`,	`copy_number`,	`status`, and	`pathology_free_text_diagnosis`. Furthermore, the column `subtype` in this table provides more insights into the samples in `NBL_MYCN_Subtype.tsv` which had a subtype `NA` (samples which could not be subtyped). If samples fell into **case 4** but didn't have a TPM value, those are subtyped as `Pathology-amp,Status-non-amp,TPM-NA`.  If a sample did not fall into any of the above cases they are subtyped as `Unclassified due to insufficient info`.

`results/QC_table.tsv`: This table has matched samples from the TARGET and GMKF cohort. Samples were matched if they shared the same `Kids_First_Participant_ID` and `sample_id`. We specifically compare the `status` between each cohort and also compare their subtpyes. 

`results/Subtypes_Based_On_Cutoff.tsv`: This table contains the biospecimens which were subtyped based on their TPM values. The subtyping of these biospecimens were done using a threshold TPM value was established in `00-Analysis-RMD`. 

`plots/TPM_Biospecimen_Matching.png`: This file contains the plot of TPM vs Biospecimen_ID for those samples which had both DNA and RNA IDs. 


`plots/TPM_Biospecimen_All_Samples_With_TMP.png`:This file contains the plot of TPM vs Biospecimen_ID for those samples which had a TPM value. 

`plots/biospecimen_ID_Amp_Del.png`: The image files labeled in the following format `biospecimen_ID_Amp_Del.png` contain the chromosome 2 segment mean for that biospecimen.










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
Neuroblastoma

Ganglioneuroblastoma

Ganglioneuroblastoma, nodular

Ganglioneuroblastoma, intermixed

Ganglioneuroma, maturing subtype OR Ganglioneuroblastoma, well differentiated

```





