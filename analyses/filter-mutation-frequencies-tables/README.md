## Filter Mutation Frequencies Tables

### Purpose
Remove `Ensembl (ESNG)` gene identifier in the OPenPedCan mutation frequency tables, including [SNV](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/snv-frequencies), [CNV](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/cnv-frequencies) and [fusion](https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/fusion-frequencies) that are not in [GENCODE v38](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/) and [Ensembl package 104](http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/).


### Analysis scripts

### `run-filter-mutation-frequencies-tables.sh`
This is a wrapper bash script for main anlysis notebook script, `filter-mutation-frequencies-tables.Rmd` that coverts JSON mutation frequencies files to JSON Line (JSONL), compresses JSONL files, and deletes intermediate JSON files. All file paths set in this script relative to the module directory. Therefore, this script should always run as if it were being called from the directory it lives in, the module directory (`OpenPedCan-analysis/analyses/filter-mutation-frequencies-tables`)


Usage:
```bash
bash run-filter-mutation-frequencies-tables.sh

```

### `filter-mutation-frequencies-tables.Rmd`
This R notebook filters SNV, CNV and fusion mutation frequencies tables to exclude Ensembl gene identifier that are not in `GENCODE v38` and `Ensembl package 104` , and lists identifiers filtered out. 

Usage:
```r
Rscript -e "rmarkdown::render('filter-mutation-frequencies-tables.Rmd', clean = TRUE)"
```

Input:
- `../snv-frequencies/results/gene-level-snv-consensus-annotated-mut-freq.tsv.gz`
- `../snv-frequencies/results/variant-level-snv-consensus-annotated-mut-freq.tsv.gz`
- `../cnv-frequencies/results/gene-level-cnv-consensus-annotated-mut-freq.tsv.gz`
- `../fusion-frequencies/results/putative-oncogene-fused-gene-freq.tsv.gz`
- `../fusion-frequencies/results/putative-oncogene-fusion-freq.tsv.gz`
- `input/input/gencode.v38.annotation.gtf.gz`
- `input/Homo_sapiens.GRCh38.104.gtf.gz` (not utilized, similar required content as gencode v38)


Results:
- `results/gene-level-snv-consensus-annotated-mut-freq.tsv.gz`
- `results/gene-level-snv-consensus-annotated-mut-freq.jsonl.gz`
- `results/variant-level-snv-consensus-annotated-mut-freq.tsv.gz`
- `results/variant-level-snv-consensus-annotated-mut-freq.jsonl.gz`
- `results/gene-level-cnv-consensus-annotated-mut-freq.tsv.gz`
- `results/gene-level-cnv-consensus-annotated-mut-freq.jsonl.gz`
- `results/putative-oncogene-fused-gene-freq.tsv.gz`
- `results/putative-oncogene-fused-gene-freq.jsonl.gz`
- `results/putative-oncogene-fusion-freq.tsv.gz`
- `results/putative-oncogene-fusion-freq.jsonl.gz`
- `filter-mutation-frequencies-tables.nb.html`

