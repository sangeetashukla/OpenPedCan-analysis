## Batch Correction of RNA-seq expression matrices

**Module authors:** Komal Rathi ([@komalsrathi](https://github.com/komalsrathi))

### Description

The goal of this analysis is to combine and batch-correct expression matrices from various RNA-seq datasets. 
The module generates combined uncorrected and batch-corrected matrices as well as combined metadata files that has the sample identifier, batch information and other specified columns only.
For QC, the module also generates UMAP/t-SNE clustering of most variable genes and housekeeping genes before and after correction as well as density plots with housekeeping genes.

There are some strict requirements for the input files for the module to function as expected:

1. `Expression matrix`: Input expression matrices must be collapsed to unique gene symbols.
2. `Metadata file`: The user should supply a histology/metadata file for each corresponding matrix to be merged.
3. `Sample identifier`: Each metadata file should have the same column name for the sample identifier that matches the columns in the corresponding expression matrix. For e.g. one histology cannot have `Kids_First_Biospecimen_ID` and the second cannot have `sample_id` for mapping sample names. If you have these differences, you can create temporary histology files for the purpose of running this module.
4. `Batch variables`: Each metadata file must have the same column name for the columns that are to be combined and used as a batch variable. For e.g. one histology cannot have `RNA_library` and the second cannot have `library_type` for mapping library names.
5. `Normal Tissue/Disease Groups`: Each metadata file must have the same column name for the columns that are to be used for to color the points in the clustering plots. For e.g. all GTEx tissues and TARGET diseases were mapped to a `group` column for the purpose of this module.

### Module Structure

```sh
.
├── README.md
├── code
│   ├── 01-batch-correct.R	# script to batch correct combined TPM matrices
│   ├── 02-qc-plots.R		# script to create t-SNE and density plots
├── output
│   ├── *_corrected.rds		# matrices for batch corrected datasets
│   ├── *_uncorrected.rds	# matrices for uncorrected datasets
│   ├── *_metadata.rds		# combined metadata with sample identifier, batch info and specified columns
├── plots
│   ├── *_density.pdf		# density plots of housekeeping genes before and after batch correction
│   ├── *_tsne.pdf		# t-SNE plots before and after batch correction with most variable and housekeeping genes
│   └── *_umap.pdf		# umap plots before and after batch correction with most variable and housekeeping genes
├── run_analysis.sh		# script to run full analysis
└── util
    ├── combine_mat.R		# function to combined input matrices for batch correction
    ├── combine_meta.R		# function to combined input metadata for batch correction
    ├── density_plots.R		# function to create density plots
    ├── pubTheme.R		# function for publication quality ggplot2 theme
    └── clustering_plots.R	# function to create UMAP/t-SNE plots
```

### Analysis scripts

#### code/01-batch-correct.R

This script combines various input matrices using their rownames. The input matrices MUST be collapsed to unique gene symbols. Batch correction is done on the combined matrices using `sva::ComBat` when input type is TPM/FPKM and `sva::ComBat_seq` when input type is expected counts. Normalized expression inputs for e.g. TPM/FPKM are log-transformed before batch-correction and so is the corresponding matrix resulting from `sva::ComBat`. These are back-transformed within the script for downstream analyses.


```sh
Rscript code/01-batch-correct.R --help

Options:
	--mat=MAT
		comma-separated list of expression matrices to combine (TPM, FPKM or expected counts) (.rds)

	--type=TYPE
		Type of expression data: TPM, FPKM or expected_count

	--metadata=METADATA
		comma-separated list of sample metadata files to combine (.rds)

	--id_col=ID_COL
		identifier column from metadata file that matches with expression matrix columns

	--batch_cols=BATCH_COLS
		comma-separated list of columns from meta file to be used to create batch variable

	--other_cols_to_keep=OTHER_COLS_TO_KEEP
		comma-separated list of columns to keep in the combined metadata in addition to the sample identifier and batch variables

	--output_prefix=OUTPUT_PREFIX
		prefix for output files
```

##### Example Run

```sh
# GTEx + TARGET NBL
Rscript code/01-batch-correct.R \
--mat 'input/gtex-gene-expression-rsem-tpm-collapsed.polya.rds, input/target-nbl-gene-expression-rsem-tpm-collapsed.polya.rds' \
--type 'tpm' \
--metadata 'input/gtex_metadata.rds, input/target_nbl_metadata.rds' \
--id_col 'sample_id' \
--batch_cols 'study_id, library_type' \
--other_cols_to_keep 'group' \
--output_prefix 'gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya'
```

##### Outputs

Output files in .rds format are generated under `output` folder:

```
output
├── gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya_corrected.rds		# batch corrected matrix
├── gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya_metadata.rds		# combined metadata with sample identifier, batch info and specified columns
└── gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya_uncorrected.rds	# uncorrected matrix
```

#### code/02-qc-plots.R

This script takes the output files generated in `01-batch-correct.R` and creates the following plots for uncorrected and corrected matrices: 

1. UMAP/t-SNE of most variable genes and housekeeping genes. UMAP is recommended for large datasets because t-SNE takes a lot of time.
2. density plots of house-keeping genes

```sh
Rscript code/02-qc-plots.R --help

Options:
	--uncorrected_mat=UNCORRECTED_MAT
		combined expression matrix with multiple batches (.rds)

	--corrected_mat=CORRECTED_MAT
		corrected expression matrix with multiple batches (.rds)

	--combined_metadata=COMBINED_METADATA
		combined metadata file with batch information (.rds)

	--var_prop=VAR_PROP
		proportion of most variable genes to be used

	--sample_id=SAMPLE_ID
		sample identifier column in metadata file matching column names in expression datasets

	--hk_genes=HK_GENES
		comma separated list of house keeping genes

	--clustering_type=CLUSTERING_TYPE
		type of clustering to use: umap or tsne

	--plots_prefix=PLOTS_PREFIX
		prefix for clustering plots
```

##### Example Run

```sh
Rscript code/02-qc-plots.R \
--uncorrected_mat 'output/gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya_uncorrected.rds' \
--corrected_mat 'output/gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya_corrected.rds' \
--combined_metadata 'output/gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya_metadata.rds' \
--var_prop 10 \
--sample_id 'sample_id' \
--hk_genes 'ACTB, TUBA1A, TUBB, GAPDH, LDHA, RPL19' \
--clustering_type 'umap' \
--plots_prefix 'gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya'
```

##### Outputs

UMAP/t-SNE and density plots are generated under the `plots` folder:

```
plots
├── gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya_density.pdf	# density plots of housekeeping genes before and after batch correction
├── gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya_umap.pdf		# UMAP plots before and after batch correction with most variable and housekeeping genes
```
