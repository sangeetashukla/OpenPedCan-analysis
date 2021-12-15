# OpenPencan Tumor Mutation Burden calculation

This analysis utilizes the SNV consensus MAF file, `../../data/snv-consensus-plus-hotspots.maf.tsv.gz` from [Pediatric Open Targets, OPenPedCan-analysis](https://github.com/PediatricOpenTargets/OpenPedCan-analysis) datasets to calculate Tumor Mutation Burden (TMB) for each experimental strategy (WGS, WXS, and Targeted Sequencing) tumor sample with SNV calls for all cohorts and cancer types evaluated in the project. The Tumor Mutation Burden calculation is adapted from [snv-callers module](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers) of the OpenPBTA-analyses, and use the SNV calls Mutect2, Strelka2, Lancet, and Vardict callers.

## TMB Calculation

For each experimental strategy and TMB calculation, the intersection of the genomic regions effectively being surveyed are used. All calls whether consensus or not are used for TMB calculations. These genomic regions are used for first filtering mutations to these regions and then for using the size in bp of the genomic regions surveyed as the TMB denominator.

### All mutations TMB

For all mutation TMBs, all callers are used. For WGS samples, the size of the genome covered by the intersection of Strelka2 and Mutect2's surveyed areas which are considered representative of all callers is used for the denominator.

```
WGS_all_mutations_TMB = (total # consensus and non-consensus snvs called by all callers) / intersection_strelka_mutect_genome_size
```
For WXS samples, the size of the genome for each the WXS bed region file is used for the denominator with the associated tumor samples.
```
WXS_all_mutations_TMB = (total # consensus and non-consensus snvs called by all callers)) / wxs_genome_size
```
### Coding only TMB

Coding only TMB uses all callers as well and the intersection demoninators are calculated by using coding sequence ranges in the gtf from Gencode 27.
This file is included in the OpenPedCan data download.
SNVs outside of these coding sequences are filtered out before being summed and used for TMB calculations as follows:

```
WGS_coding_only_TMB = (total # consensus and non-consensus snvs called by all callers) / intersection_wgs_strelka_mutect_CDS_genome_size
```
For WXS samples, each the WXS bed region file is intersected with the coding sequences for filtering and for determining the denominator to be used with the with the associated tumor samples.
```
WXS_coding_only_TMB = (total # consensus and non-consensus snvs called by all callers) / intersection_wxs_CDS_genome_size
```

## General usage of scripts


#### `run_tmb_calculation.sh`
This is a bash script wrapper for setting input file paths for the main anlysis script, `01-calculate_tmb.R` and creating additional intermedaite input files `../../scratch/intersect_strelka_mutect2_vardict_WGS.bed` and `gencode.v27.primary_assembly.annotation.bed`. All file paths set in this script relative to the module directory. Therefore, this script should always run as if it were being called from the directory it lives in, the module directory (`OpenPedCan-analysis/analyses/tmb-calculation`).


#### 01-calculate_tmb.R
Uses the OpenPedCan snv consensus file to calculate TMB for all WGS, WXS, and Targeted Sequencing samples.Two TMB files are created, one including *all snv* and a *coding snvs only*, these both are made using both consensus and non-consensus mutations called by all callers.

**Argument descriptions**
```
 --consensus_maf_file : Input OpenPedCan SNV consensus MAF file.

 --bed_files : Input samples to target BED mapping file.

 --histologies_file : Input OpenPedCan histologies metadata file.

 --coding_regions : BED file for coding regions to use for coding only TMB.

 --nonsynfilter_maf : If TRUE, filter out synonymous mutations, keep
                 non-synonymous mutations, according to maftools definition.
                 Default is FALSE

 --nonsynfilter_focr : If TRUE, filter out synonymous mutations, keep non-synonymous
                 mutations, according to Friends of Cancer esearch definition.
                 Default is FALSE
```

#### `Util/split_mnv.R`
Contains a function to split multinucleotide variants (MNVs) into single nucleotide variants (SNVs).

#### `Util/tmb_functions.R`
Contains functions for calculating Tumor Mutation Burden (TMB).


## Analysis input file

### OpenPedCan data download files:
- `../../data/snv-consensus-plus-hotspots.maf.tsv.gz`
- `../../data/gencode.v27.primary_assembly.annotation.gtf.gz`
- `../../data/histologies.tsv`

### Module specific input files:
- `input/biospecimen_id_to_bed_map.txt`
- `input/S0274956_Padded_HG38.merged.bed`
- `input/S02972011_Covered_hg38_100.bed`
- `input/S04380110_Regions_hg38_100.bed`
- `input/S07604715_100bp_Padded.bed`
- `input/SeqCap_EZ_Exome_v2_Padded_HG38.merged.bed`
- `input/StrexomeLite_hg38_liftover_100bp_padded.bed`
- `input/Strexome_targets_intersect_sorted_padded100.GRCh38.bed`
- `input/TARGET_AML_NBL_WT_SeqVal79_attempt06_AllTracks_HG38_bed_expanded100.bed`
- `input/ashion_confidential_exome_v2_targets_hg38_paded100.bed`
- `input/hg38_strelka.bed`
- `input/nexterarapidcapture_exome_targetedregions_v1.2_hg38_100.bed`
- `input/truseq-exome-targeted-regions-manifest-v1-2_hg38_100.bed`
- `input/wgs_canonical_calling_regions.hg38.bed`
- `input/xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.bed`


## Analysis result files

### Output:
- `results/snv-mutation-tmb-coding.tsv`
- `results/snv-mutation-tmb-all.tsv`
