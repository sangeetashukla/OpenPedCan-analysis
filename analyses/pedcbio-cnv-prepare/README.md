## Modify annotated CNV files for pedcBio upload

**Module author:** Run Jin ([@runjin326](https://github.com/runjin326))

Currently, annotated CNV files `consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz` and `consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz` (and the combined file `consensus_wgs_plus_cnvkit_wxs.tsv.gz`) only contain samples and segments that have CNV calls.
To capture total number of samples profiled in PedCBio, this analysis does the following:
1) for samples without CNV calls, status == `neutral` was added to the results for all available segments
2) for samples with CNV calls, for diploid segments that were not present in the annotation file, status == `neutral` was added to the results for those segments

When adding back neutral segments to X or Y segments, based on the `germline_sex_estimate` (or `reported_gender` when `germline_sex_estimate` is not available), copy number and tumor ploidy were added as followed:
1) For female, `copy_number=2`, `ploidy=2` and `status=neutral` were assigned to X segments 
2) For female, `copy_number=0`, `ploidy=0` and `status=neutral` were assigned to Y segments 
3) For male, `copy_number=1`, `ploidy=1` and `status=neutral` were assigned to both X and Y segments 

Additionally, in the notebook, there are several steps that output the number of segment entries for each sample to make sure the right amount of segments were rescued. 

Usage:
```
Rscript -e "rmarkdown::render('pedcbio_cnv_prepare.Rmd', clean = TRUE)"

```

Input:
- `../../data/consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz`
- `../../data/consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz`
- `../../data/cnv-cnvkit.seg.gz`
- `../../data/cnv-consensus.seg.gz`
- `../../data/histologies.tsv`

Output:
- `results/consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz`
- `results/consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz`
- `results/consensus_wgs_plus_cnvkit_wxs`

The output files are directly uploaded to S3 buckets for loading into PedCBio.
