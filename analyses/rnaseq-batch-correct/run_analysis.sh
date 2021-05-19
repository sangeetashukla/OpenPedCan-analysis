# GTEx + TARGET NBL (correct by study_id and library_type)
# Rscript code/01-batch-correct.R \
# --mat 'input/gtex-gene-expression-rsem-tpm-collapsed.polya.rds, input/target-nbl-gene-expression-rsem-tpm-collapsed.polya.rds' \
# --type 'tpm' \
# --metadata 'input/gtex-histologies.tsv, input/target-nbl-histologies.tsv' \
# --id_col 'sample_id' \
# --batch_cols 'study_id, library_type' \
# --other_cols_to_keep 'group' \
# --output_prefix 'gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya'

# UMAP plots
Rscript code/02-qc-plots.R \
--uncorrected_mat 'output/gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya_uncorrected.rds' \
--corrected_mat 'output/gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya_corrected.rds' \
--combined_metadata 'output/gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya_metadata.tsv' \
--var_prop 10 \
--sample_id 'sample_id' \
--hk_genes 'https://housekeeping.unicamp.br/Housekeeping_GenesHuman.csv' \
--clustering_type 'umap' \
--plots_prefix 'gtex-target-nbl-gene-expression-rsem-tpm-collapsed.polya'

# GTEx + TARGET (correct by study_id and library_type)
# Rscript code/01-batch-correct.R \
# --mat 'input/gtex-gene-expression-rsem-tpm-collapsed.polya.rds, input/target-gene-expression-rsem-tpm-collapsed.rds' \
# --type 'tpm' \
# --metadata 'input/gtex-histologies.tsv, input/target-histologies.tsv' \
# --id_col 'sample_id' \
# --batch_cols 'study_id, library_type' \
# --other_cols_to_keep 'group' \
# --output_prefix 'gtex-target-gene-expression-rsem-tpm-collapsed'

# UMAP plots
Rscript code/02-qc-plots.R \
--uncorrected_mat 'output/gtex-target-gene-expression-rsem-tpm-collapsed_uncorrected.rds' \
--corrected_mat 'output/gtex-target-gene-expression-rsem-tpm-collapsed_corrected.rds' \
--combined_metadata 'output/gtex-target-gene-expression-rsem-tpm-collapsed_metadata.tsv' \
--var_prop 10 \
--sample_id 'sample_id' \
--hk_genes 'https://housekeeping.unicamp.br/Housekeeping_GenesHuman.csv' \
--clustering_type 'umap' \
--plots_prefix 'gtex-target-gene-expression-rsem-tpm-collapsed'

# GTEx + TARGET (only correct by library_type)
# Rscript code/01-batch-correct.R \
# --mat 'input/gtex-gene-expression-rsem-tpm-collapsed.polya.rds, input/target-gene-expression-rsem-tpm-collapsed.rds' \
# --type 'tpm' \
# --metadata 'input/gtex-histologies.tsv, input/target-histologies.tsv' \
# --id_col 'sample_id' \
# --batch_cols 'library_type' \
# --other_cols_to_keep 'group' \
# --output_prefix 'gtex-target-gene-expression-rsem-tpm-collapsed-libraryonly'

# UMAP plots
Rscript code/02-qc-plots.R \
--uncorrected_mat 'output/gtex-target-gene-expression-rsem-tpm-collapsed-libraryonly_uncorrected.rds' \
--corrected_mat 'output/gtex-target-gene-expression-rsem-tpm-collapsed-libraryonly_corrected.rds' \
--combined_metadata 'output/gtex-target-gene-expression-rsem-tpm-collapsed-libraryonly_metadata.tsv' \
--var_prop 10 \
--sample_id 'sample_id' \
--hk_genes 'https://housekeeping.unicamp.br/Housekeeping_GenesHuman.csv' \
--clustering_type 'umap' \
--plots_prefix 'gtex-target-gene-expression-rsem-tpm-collapsed-libraryonly'

# GTEx + GMKF (correct by study_id and library_type)
# Rscript code/01-batch-correct.R \
# --mat 'input/gtex-gene-expression-rsem-tpm-collapsed.polya.rds, input/kfnbl-gene-expression-rsem-tpm-collapsed.stranded.rds' \
# --type 'tpm' \
# --metadata 'input/gtex-histologies.tsv, input/kfnbl-histologies.tsv' \
# --id_col 'sample_id' \
# --batch_cols 'study_id, library_type' \
# --other_cols_to_keep 'group' \
# --output_prefix 'gtex-kfnbl-gene-expression-rsem-tpm-collapsed'

# UMAP plots
Rscript code/02-qc-plots.R \
--uncorrected_mat 'output/gtex-kfnbl-gene-expression-rsem-tpm-collapsed_uncorrected.rds' \
--corrected_mat 'output/gtex-kfnbl-gene-expression-rsem-tpm-collapsed_corrected.rds' \
--combined_metadata 'output/gtex-kfnbl-gene-expression-rsem-tpm-collapsed_metadata.tsv' \
--var_prop 10 \
--sample_id 'sample_id' \
--hk_genes 'https://housekeeping.unicamp.br/Housekeeping_GenesHuman.csv' \
--clustering_type 'umap' \
--plots_prefix 'gtex-kfnbl-gene-expression-rsem-tpm-collapsed'
