## Add formatted sample id column for PedCBio upload

**Module author:** Run Jin ([@runjin326](https://github.com/runjin326))

Currently, for some of the samples, when multiple DNA or RNA specimens are associated with the same sample, there 
is no column that would distinguish between different aliquots while still tieing DNA and RNA together.
This module adds a column called `formatted_sample_id` where the base name is the sample id and additional `tiebreaks` were added when multiple RNA or DNA samples are associated with the same participant.

For PBTA samples, `sample_id` column is used as the basename
- Using `sample_id` column, we can tie all DNA and RNA samples together
- Using `formatted_sample_id` column, we can distinguish amongst multiple DNA or RNA samples 

For TARGET, TCGA, and GTEx samples, `Kids_First_Participant_ID` column is used as the basename
- Using `Kids_First_Participant_ID` column, we can tie all DNA and RNA samples together
- Using `formatted_sample_id` column, we can distinguish amongst multiple DNA or RNA samples 

Usage:
  ```
Rscript -e "rmarkdown::render('pedcbio_sample_name_col.Rmd', clean = TRUE)"

```
or
```
bash run_add_name.sh
```

Input:
- `input/cbtn_cbio_sample.csv`
- `input/oligo_nation_cbio_sample.csv`
- `input/dgd_cbio_sample.csv`
- `input/x01_fy16_nbl_maris_cbio_sample.csvz`

Output:
- `results/histologies-formatted-id-added.tsv`

The output files are directly uploaded to S3 buckets for loading into PedCBio.
