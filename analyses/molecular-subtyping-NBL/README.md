# Molecular Subtyping NBL
**Module authors**: Aditya Lahiri ([@adilahiri](https://github.com/adilahiri)), Eric Wafula ([@ewafula](https://github.com/ewafula)), and Jo Lynne Rokita ([@jharenza](https://github.com/jharenza))


Objective: To molecularly subtype neuroblastoma, ganglioneuroblastoma, and ganglioneuroma samples into MYCN 
amplified or MYCN non-amplified ([Ticket#417](https://github.com/PediatricOpenTargets/ticket-tracker/issues/417)).
MYCN is located on chromosone 2p and the subtyping is done by checking the CNV calls and the
clinically-reported MYCN status. We also perform QC by checking the consistency between TARGET and GMKF NBL samples.

**Subtyping criteria**:

**case 1**:
If `pathology_free_text_diagnosis` is **amplified** and `status` is **amplified** assign subtype as `NBL, MYCN amplified`

**case 2**:
If `pathology_free_text_diagnosis` is **non-amp** and `status` is **amplified** assign subtype as `NBL, MYCN amplified`

**case 3**:
If `pathology_free_text_diagnosis` is **non-amp** and `status` is **non-amp** assign subtype as `NBL, MYCN non-amplified`

**case 4**:
If `pathology_free_text_diagnosis` is **amplified** and status is **non-amp**
For such samples if there exists a TPM value, evaluate if the TPM is above or below 140.83 (`Suggested_Cutoff`) established in the image `results/TPM_Biospecimen_All_Samples_With_TMP.png` assign subtype as `NBL, MYCN amplified` or `NBL, MYCN non-amplified` respectively.  In case there is no TPM values then assign the subtype as `NA`. 

**case 5**:
If there are samples that are not yet subtyped but have a TPM value, assign them a subtype of `NBL, to be classified`.

**case 6**:
Other remaining samples are not subtyped and the subtype field is left as `NBL, to be classified`.

**Note**: The files in the input folder are NBL MCYN `clinical patient-status mapping files` for GMKF and TARGET samples. In V11 some values of the `pathology_free_text_diagnosis` are missing, these input files fill in some of these missing values. 

# Usage
To run the module, execute the following command from the command line in the directory `OpenPedCan-analysis/analyses/molecular-subtyping-NBL` 

`bash run_nbl_subtyping.sh `

# Module Content
1. `00-Preprocessing.Rmd`: This script creates a table called ``Alteration_Table.tsv` which contains the biospecimens 
for neuroblastama, ganglioneuroblastoma, and ganglioneuroma which need to be subptyped. This table is written into the intermediate folder. 

2. `01-Subtyping.Rmd`: This script loads the table `Alteration_Table.tsv` from the intermediate fold performs the subtyping. 

3. `run_nbl_subtyping.sh`: Bash script to execute `00-Preprocessing.Rmd` and `01-Subtyping.Rmd`.

4. `input/gmkf_patient_clinical_mycn_status.tsv`: Has values for `pathology_free_text_diagnosis` in the GMKF cohort. These values are
    not present in the V11 data and will be made available in V12.

5. `input/target_patient_clinical_mycn_status.tsv`: Has values for `pathology_free_text_diagnosis` in the TARGET cohort. These values are not present in the V11 data and will be made available in V12.

6. `intermediate/Alteration_Table.tsv`: Table containing NBL biospecimens that will be subtyped in `01-Subtyping.Rmd`.
    
7. `results/NBL_MYCN_Subtype.tsv`: This table is the main output of this module. For each biospecimen sample, the table  lists its DNA and RNA IDs (if available) and the subtype. 

8. `results/Alteration_Table.tsv`: This table is similar to the table `NBL_MYCN_Subtype.tsv`. However, this table has 
additional columns which contain information on `MYCN_TPM`,	`copy_number`,	`status`, and	`pathology_free_text_diagnosis`. Furthermore, the column `subtype` in this table provides more insights into the samples in `NBL_MYCN_Subtype.tsv` which had a subtype `NA` (samples which could not be subtyped). If samples fell into **case 4** but didn't have a TPM value, those are subtyped as `Pathology-amp,Status-non-amp,TPM-NA`.  If a sample did not fall into any of the above cases they are subtyped as `Unclassified due to insufficient info`.

9. `results/QC_table.tsv`: This table has matched samples from the TARGET and GMKF cohort. Samples were matched if they shared the same `Kids_First_Participant_ID` and `sample_id`. We specifically compare the `status` between each cohort and also compare their subtpyes. 

10. `results/Subtypes_Based_On_Cutoff.tsv`: This table contains the biospecimens which were subtyped based on their TPM values. The subtyping of these biospecimens were done using a threshold TPM value was established in `00-Preprocessing-Rmd`. 

11. `plots/TPM_Biospecimen_Matching.png`: This file contains the plot of TPM vs Biospecimen_ID for those samples which had both DNA and RNA IDs. 


12. `plots/TPM_Biospecimen_All_Samples_With_TMP.png`:This file contains the plot of TPM vs Biospecimen_ID for those samples which had a TPM value. 

13. `plots/biospecimen_ID_Amp_Del.png`: The image files labeled in the following format `biospecimen_ID_Amp_Del.png` contain the chromosome 2 segment mean for that biospecimen.


# Brief description of the module pipeline:
This module has only two notebooks `00-Preprocessing.Rmd` and  `01-Subtyping.Rmd`. `00-Preprocessing.Rmd` gathers the biospecimens that need to subtyped, and the 
subtyping is done in `01-Subtyping.Rmd` based on the criteria listed above. The data files being used in this analysis are: 
```
1. histologies-base.tsv

2. consensus_wgs_plus_cnvkit_wxs.tsv.gz

3. cnv-cnvkit.seg.gz

4. gene-expression-rsem-tpm-collapsed.rds
```
The module extracts the NBL tumor samples from the `histologies-base.tsv` file by filtering for the following fields:

| Histology File Column Name                                                         |Filter columnsthe following values| 
|:----------------------------------------------------------------------------|------------------------------------------:|
|Pathology_Diagnosis|`Neuroblastoma`,`Ganglioneuroblastoma`, `Ganglioneuroblastoma, nodular`, `Ganglioneuroblastoma, intermixed` , `Ganglioneuroma, maturing subtype OR Ganglioneuroblastoma, well differentiated`|                                                      
|sample_type|`Tumor`|   
|experimental_strategy|`WGS`,`WXS`,`Targeted Sequencing`, `RNA-Seq`|




We then filter the `consensus_wgs_plus_cnvkit_wxs.tsv.gz` and `gene-expression-rsem-tpm-collapsed.rds` for the gene symbol `MYCN` and then join it with the filtered histology file. We use this composite file to get the DNA and RNA biospecimen IDs for the records and then subtype them based on the above mentioned subtyping criteria. 



