# OpenPedCan Methylation Analysis

## Purpose

Preprocess probe hybridization intensity values of selected methylated and unmethylated cytosine (CpG) loci into usable methylation measurements for the [Pediatric Open Targets, OPenPedCan-analysis](https://github.com/PediatricOpenTargets/OpenPedCan-analysis) raw DNA methylation array datasets. 

## Data description 

#### TARGET
The TARGET Illumina `Infinium HumanMethylation (450k and 27k) BeadChip` and `Roche Nimblegen HELP` methylation arrays and sample metadata for following pediatric cancers were downloaded from the [TARGET project website](https://ocg.cancer.gov/programs/target/target-methods):
 
- [Neuroblastoma (NBL)](https://target-data.nci.nih.gov/Public/NBL/methylation_array/) - `Infinium HumanMethylation450`
- [Osteosarcoma (OS)](https://target-data.nci.nih.gov/Public/OS/methylation_array/) - `Infinium HumanMethylation450`
- [Clear Cell Sarcoma of the Kidney (CCSK)](https://target-data.nci.nih.gov/Public/CCSK/methylation_array/) - `Infinium HumanMethylation450`
- [Wilms Tumor (WT)](https://target-data.nci.nih.gov/Public/WT/methylation_array/) - `Infinium HumanMethylation450`
- [Acute Myeloid Leukemia (AML)](https://target-data.nci.nih.gov/Public/AML/methylation_array/) - `Infinium HumanMethylation450 and Infinium HumanMethylation27`
- [Acute Lymphoblastic Leukemia (ALL)](https://target-data.nci.nih.gov/Public/ALL/methylation_array/Phase2/) - `Nimblegen HELP methylation`

## Analysis

- The TARGET Illumina methylation analysis results available on the [TARGET project website](https://ocg.cancer.gov/programs/target/target-methods) were preprocessed with different methylation software packages, including `GenomeStudio` (AML27k), `minfi` (AML450k), `BeadStudio` (CCSK and WT), `methylumi` (NBL), and `Lumi+BMIQ` (OS). We have re-analyzed all the cancer types with Illumina arrays using the updated version of [minfi Bioconductor package](https://academic.oup.com/bioinformatics/article/33/4/558/2666344). We utilized the `preprocessQuantile` array preprocessing method implemented in minfi for background correction, removing of outlier samples, and performing stratified quantile normalization to estimate the relative abundance of methylated and unmethylated cytosines at selected loci (Beta and M values). Using the one prepprocessing method for all datasets from the same array platform/vendor allows downstream comparisons of methylation profiles across cancer type.
- A small subset of `non-tumor (normal)` samples (12) only available in the `Neuroblastoma (NBL)` array dataset were processed seperately and could protentially be utlized as control set for both Neuroblastoma and other cancer types. 
- Unlike the probe annotations Illumina 450k arrays, the probe annotations for older Illumina 27k arrays do not specify methylation region relative to gene coding loci (TSS or body) and are therefore not reported in the summary results. Since CpG sites interogated by the Illumina 27k array probes are a subset of Illumina 450k array probes, the methylation regions of the former can be easily be inferred from the latter for common samples among the two. 
- The`Acute Lymphoblastic Leukemia (ALL)` array dataset is from a different array platform and vendor, `Nimblegen HELP (Roche)` has not been reanalyzed.
- All sample array data files have been rename using their respective `TARGET sample USI IDs` to be consistent for this analysis module.
- Reported results have been filtered to exclude array probes without annotation information.


## General usage of scripts


#### `run-preprocess-illumina-arrays.sh`
This is a bash script wrapper for setting input file paths for the main analysis script, `01-preprocess-illumina-arrays.R` All file paths set in this script relative to the module directory. Therefore, this script should always run as if it were being called from the directory it lives in, the module directory (`OpenPedCan-analysis/analyses/methylation-analysis`).

```
bash run-preprocess-illumina-arrays.sh
```

#### `01-preprocess-illumina-arrays.R`
Preprocesses raw Illumina Infinium HumanMethylation BeadArrays (27K, and 450K) intensities using `minfi Bioconductor package` into usable methylation measurements (Beta and M values) for TARGET normal and tumor samples. 

Addition data packages for `Ilumina Infinium HumanMethylation27 BeadArrays` that don't with for analyzing `minfi` will need to be installed separately as folllows:
- `BiocManager::install("IlluminaHumanMethylation27kmanifest")`
- `BiocManager::install("IlluminaHumanMethylation27kanno.ilmn12.hg19")`

**Argument descriptions:**
```
Usage: 01-preprocess-illumina-arrays.R [options]


Options:
	--base_dir=CHARACTER
		The absolute path of the base directory containing sample 
              array IDAT files.

	--metadata_file=CHARACTER
		The metedata file associated with the sample array data
              files.

	--preprocess_method=CHARACTER
		preprocesses the Illumina methylation array using the of
              the following minfi methods: preprocessQuantile, preprocessFunnorm, 
              or preprocessIllumina.
              Default is preprocessQuantile

	--snp_filter
		If TRUE, drops the probes that contain either a SNP at
              the CpG interrogation or at the single nucleotide extension.
              Default is TRUE

	-h, --help
		Show this help message and exit
```

#### `002-compare-cancers-tsne-and-umap-plots.R`
Creates comparison methylation profiles t-SNE and UMAP plots of selected cancer genes for cancer types preprocessed from Illumina Infinium HumanMethylation450 BeadArrays. Utilizes a list of cancer gene symbols (`metadata/TARGET_Methylation_GeneList.txt`) that are present in the array annotation design reported in the methylation M-values result tables (in the `results` folder) produced by the `01-preprocess-illumina-arrays.R` script.

```
Rscript --vanilla 02-compare-cancers-tsne-and-umap-plots.R
```

#### `03-summarize-methylation-values.R`
Summarizes methylation `Beta-values` and `M-values` to obtain a representative `median value` per CpG probe site for all samples from each cancer type preprocessed from Illumina Infinium HumanMethylation450 BeadArrays reported in a single combined file, `results/median-beta-values-methylation.tsv.gz` and `results/median-m-values-methylation.tsv.gz` respectively. It also annotates the CpG probe sites with the associated genes from the evidence-based annotation of the human genome (hg19/GRCh37), GENCODE version 19 (Ensembl 74), on which the array CpG site coordinates are based. 

```
Rscript --vanilla 03-summarize-methylation-values.R
```

#### `04-get-array-gencode-annotations.R`
Parses GENCODE annotation of the human genome (GRCh37), version 19 (Ensembl 74), and associated feature coordinates for genes in the Illumina Infinium HumanMethylation450 BeadArrays annotation design (`results/methylation-array-gencode.annotations.tsv.gz`) to use for plotting comparative methylation profiles with GViz R package in association with any one of the two results tables of representative median CpG probe sites methylation values. 

```
Rscript --vanilla 04-get-array-gencode-annotations.R
```

## Input datasets

#### `Methylation arrays:`
Methylation array datasets are avaliable on the CHOP HPC `Isilon` sever (location: `/mnt/isilon/opentargets/wafulae/methylation-analysis/data/`). Please contact `Avin Farrel (@afarrel)` for access. 
- `data/Normal/*.idat` - Normal samples array files
- `data/NBL/*.idat` - Neuroblastoma (NBL) tumor samples array files
- `data/OS/*.idat` - Osteosarcoma (OS) tumor samples samples array files
- `data/CCSK/*.idat` - Clear Cell Sarcoma of the Kidney (CCSK) tumor samples array
- `data/WT/*.idat` - Wilms Tumor (WT) tumor samples array
- `data/AML/*.idat` - Acute Myeloid Leukemia (AML) tumor samples array


#### `Sample and reference metadata:`
- `metadata/TARGET_Normal_MethylationArray_20160812.sdrf.txt` - metadata for Normal samples
- `metadata/TARGET_NBL_MethylationArray_20160812.sdrf.1.txt` - metadata for Neuroblastoma (NBL) tumor samples, batch 1
- `metadata/TARGET_NBL_MethylationArray_20160812.sdrf.2.txt` - metadata for Neuroblastoma (NBL) tumor samples, batch 2
- `metadata/TARGET_OS_MethylationArray_20161103.sdrf.txt` - metadata for Osteosarcoma (OS) tumor samples
- `metadata/TARGET_CCSK_MethylationArray_20160819.sdrf.txt` - metadata for Clear Cell Sarcoma of the Kidney (CCSK) tumor samples
- `metadata/TARGET_WT_MethylationArray_20160831.sdrf.txt` - metadata for Wilms Tumor (WT) tumor samples
- `metadata/TARGET_AML_MethylationArray_20160812_450k.sdrf.1.txt` - metadata for Acute Myeloid Leukemia (AML) tumor samples, 450k arrays, batch 1
- `metadata/TARGET_AML_MethylationArray_20160812_450k.sdrf.2.txt` - metadata for Acute Myeloid Leukemia (AML) tumor samples, 450k arrays, batch 2
- `metadata/TARGET_AML_MethylationArray_20160812_27k.sdrf.1.txt` - metadata for Acute Myeloid Leukemia (AML) tumor samples, 27k arrays, batch 1
- `metadata/TARGET_AML_MethylationArray_20160812_27k.sdrf.2.txt` - metadata for Acute Myeloid Leukemia (AML) tumor samples, 27k arrays, batch 2
- `metadata/TARGET_AML_MethylationArray_20160812_27k.sdrf.3.txt` - metadata for Acute Myeloid Leukemia (AML) tumor samples, 27k arrays, batch 3
- `metadata/TARGET_Methylation_GeneList.txt` - a selected set of genes expressed in a subset of the tumor types (contact @afarrel for details)
- `metadata/UCSC_hg19-GRCh37_Ensembl2RefSeq.tsv` - [UCSC hg19/GRCh37 Ensemble to RefSeq gene IDs mapping file](https://genome.ucsc.edu/cgi-bin/hgTables)
- `metadata/gencode.v19.annotation.gtf.gz` [GENCODE evidence-based annotation of the human genome (GRCh37), version 19 (Ensembl 74)](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/)

## Results
Summary result files of methylation `beta-values` and `M-values` are too large to upload to this repository and available on the CHOP HPC `Isilon` sever (location: `/mnt/isilon/opentargets/wafulae/methylation-analysis/results/`). Please contact `Avin Farrel (@afarrel)` for access. Only the median `beta-values` and `M-values` are uploaded to this repository.
- `results/Normal-beta-values-methylation.tsv.gz`
- `results/Normal-m-values-methylation.tsv.gz`
- `results/NBL-beta-values-methylation.tsv.gz`
- `results/NBL-m-values-methylation.tsv.gz`
- `results/OS-beta-values-methylation.tsv.gz`
- `results/OS-m-values-methylation.tsv.gz`
- `results/CCSK-beta-values-methylation.tsv.gz`
- `results/CCSK-m-values-methylation.tsv.gz`
- `results/WT-beta-values-methylation.tsv.gz`
- `results/WT-m-values-methylation.tsv.gz`
- `results/AML450k-beta-values-methylation.tsv.gz`
- `results/AML450k-m-values-methylation.tsv.gz`
- `results/AML27k-beta-values-methylation.tsv.gz`
- `results/AML27k-m-values-methylation.tsv.gz`
- `results/median-beta-values-methylation.tsv.gz`
- `results/median-m-values-methylation.tsv.gz`
- `results/methylation-array-gencode.annotations.tsv.gz`

## Plots
Comparison methylation t-SNE and UMAP plots among cancer types for selected genes expressed in a subset of the tumor type. ``Plots` and `M-values` gene matrices are available on the CHOP HPC `Isilon` sever (location: `/mnt/isilon/opentargets/wafulae/methylation-analysis/plots/`). Please contact `Avin Ferrel (@afarrel)` for access.
- `plots/<GeneSymbol>-m-values.csv` - M-values gene matrices 
- `plots/<GeneSymbol>-plot.png` - gene t-SNE plots
- `plots/<GeneSymbol>-plot.png` - gene UMAP plot

