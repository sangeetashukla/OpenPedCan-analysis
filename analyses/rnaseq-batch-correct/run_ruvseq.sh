# matched samples in PBTA cohort
# 4 CNS Embryonal tumor, 4 Diffuse midline glioma and 10 High-grade glioma/astrocytoma
# sample size is okay in HGG so we will use that as a test
Rscript code/03-test_ruvseq.R \
--dataset match_pbta_hgg \
--k_value 5

# matched samples in TARGET cohort
# 24 Acute Lymphoblastic Leukemia samples have been sequenced with > 1 technology
Rscript code/03-test_ruvseq.R \
--dataset match_target_all \
--k_value 5