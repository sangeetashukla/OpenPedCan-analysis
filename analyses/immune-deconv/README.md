## Immune Deconvolution

**Module authors:** Komal S. Rathi ([@komalsrathi](https://github.com/komalsrathi))

### Description

The goal of this analysis is to use the R package `immunedeconv` to quantify and compare various immune cell types in the tumor microenvironment (TME) across various cancer and GTEx groups. 
The package `immunedeconv`, provides the following deconvolution methods: `xCell` (n = 64; immune and non-immune cell types), `CIBERSORT` (relative mode; n = 22 immune cell types); `CIBERSORT (abs.)` (absolute mode; n = 22 immune cell types), `TIMER` (n = 6), `EPIC` (n = 6), `quanTIseq` (n = 10) and `MCP-Counter` (n = 8). 

Both `CIBERSORT` and `CIBERSORT (abs.)` require two files i.e. `LM22.txt` and `CIBERSORT.R`, that are available upon request from https://cibersort.stanford.edu/. Please refer to https://icbi-lab.github.io/immunedeconv/articles/immunedeconv.html#special-case-cibersort for more details. We recommend using `xCell` instead when these files are not available to the user. 

### Method selection

1. Between sample comparisons: If the user is interested in between-sample comparisons, we recommend using `xCell` and `CIBERSORT (abs.)`. `xCell` is our method of choice because it is the most comprehensive deconvolution method (with the largest representation of cell types) and widely used in literature vs other deconvolution methods. Reference for benchmarking between xCell and other methods: `PMID: 31641033`

2. Between cell-type comparisons: If the user is interested in between-cell-type comparisons, it is recommended to use `CIBERSORT` and `CIBERSORT (abs.)`. 

For our analysis, we will pick the top two methods that represent the most cell types i.e. `xCell` and `CIBERSORT`. `xCell` and `CIBERSORT (abs.)` output enrichment scores and `CIBERSORT` outputs actual cell fractions so we will use `xCell` and `CIBERSORT (abs.)` to enable comparisons across various samples or groups. Reference: `PMID: 31510660`

### Analysis scripts

#### 01-immune-deconv.R

1. Inputs from data download

```
gene-expression-rsem-tpm-collapsed.rds
data/histologies.tsv
```

2. Function

This script deconvolutes immune cell types using the method of choice for e.g. `xCell` or `CIBERSORT (abs)`. Since, `xCell` uses the variability among the samples for a linear transformation of the output score, we split the expression matrix into individual `cohorts + cancer_group` or `cohort + gtex_group` and deconvolute them separately. Once processed, all data is combined into a single rds file which can be used to do comparisons across various groups.

3. Output: 

```
results/{deconv_method}_output.rds
````

The results in the rds file are predicted immune scores per cell type per input sample. These scores are not actual cell fractions but arbitrary scores representing enrichment of the cell types which can be compared across various cancer/gtex groups. Depending on the user requirement, the output can be used to create various visualizations. 

#### 02-summary-plots.R 

1. Input

```
results/{deconv_method}_output.rds
```

2. Function:

This script creates heatmaps from predicted immune scores.

3. Output

* `plots/{deconv_method}_heatmap_by_group.pdf`: heatmap of average immune scores per cell type per cancer group or GTEx group.

### Running the analysis

The following script will run the full analysis using two methods of choice: `xCell` and `CIBERSORT (abs)`:

```
RUN_CIBERSORT=1 bash run-immune-deconv.sh
```

As we cannot make the `LM22.txt` and `CIBERSORT.R` publicly available, so we will not run CIBERSORT via CI:

```
bash run-immune-deconv.sh
```
