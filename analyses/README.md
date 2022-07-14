## Analysis Modules

This directory contains various analysis modules in the OpenPedCan project. See the README of an individual analysis modules for more information about that module.

### Modules at a glance

The table below is intended to help project organizers quickly get an idea of what files (and therefore types of data) are consumed by each analysis module, what the module does, and what output files it produces that can be consumed by other analysis modules. This is in service of documenting interdependent analyses. Note that *nearly all* modules use the harmonized clinical data file (`histologies.tsv`) even when it is not explicitly included in the table below.

<table>
<colgroup>
<col style="width: 10%" />
<col style="width: 10%" />
<col style="width: 10%" />
<col style="width: 10%" />
<col style="width: 10%" />
<col style="width: 10%" />
<col style="width: 10%" />
<col style="width: 10%" />
<col style="width: 10%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th><p>Module</p></th>
<th><p>Input Files</p></th>
<th><p>Brief Description</p></th>
<th><p>Produces files for data release?</p></th>
<th><p>Output Files Consumed by Other Analyses</p></th>
<th><p>Modules Consume Outputs</p></th>
<th><p>OT compatibility</p></th>
<th><p>Adapted for OPC?</p></th>
<th><p>Run Platform</p></th>
<th><p>Action Plan</p></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/chromosomal-instability"><code>chromosomal-instability</code></a></p></td>
<td><p><code>histologies.tsv</code> <code>sv-manta.tsv.gz</code> <code>cnv-cnvkit.seg.gz</code></p></td>
<td><p>Evaluates chromosomal instability by calculating chromosomal breakpoint densities and by creating circular plot visuals</p></td>
<td><p>No</p></td>
<td><p><code>breakpoint-data/union_of_breaks_densities.tsv</code></p></td>
<td><p><code>molecular-subtyping-EPN</code></p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>Will Adapt for OT</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/chromothripsis"><code>chromothripsis</code></a></p></td>
<td><p><code>sv-manta.tsv.gz</code> <code>cnv-consensus.seg.gz</code> <code>independent-specimens.wgs.primary-plus.tsv</code></p></td>
<td><p>chromothripsis analysis per <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1007">#1007</a></p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/cnv-chrom"><code>cnv-chrom-plot</code></a></p></td>
<td><p><code>cnv-consensus-gistic.zip</code> <code>cnv-consensus.seg</code></p></td>
<td><p>Plots genome wide visualizations relating to copy number results</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/cnv-comparison"><code>cnv-comparison</code></a></p></td>
<td><p>Earlier version of SEG files</p></td>
<td><p><em>Deprecated</em>; compared earlier version of the CNV methods.</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/cnv-frequencies"><code>cnv-frequencies</code></a></p></td>
<td><p><code>histologies.tsv</code> <code>consensus_wgs_plus_cnvkit_wxs.tsv.gz</code> <code>independent-specimens.wgswxspanel.primary.eachcohort.tsv</code> <code>independent-specimens.wgswxspanel.relapse.eachcohort.tsv</code> <code>independent-specimens.wgswxspanel.primary.tsv</code> <code>independent-specimens.wgswxspanel.relapse.tsv</code></p></td>
<td><p>Annotate CNV table with mutation frequencies</p></td>
<td><p>No</p></td>
<td><p><code>results/gene-level-cnv-consensus-annotated-mut-freq.jsonl.gz</code> <code>results/gene-level-cnv-consensus-annotated-mut-freq.tsv.gz</code></p></td>
<td><p>Upload to FNL BOX</p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/collapse-rnaseq"><code>collapse-rnaseq</code></a></p></td>
<td><p><code>gene-expression-rsem-tpm.rds</code> <code>gencode.v27.primary_assembly.annotation.gtf.gz</code></p></td>
<td><p>Collapses RSEM FPKM matrices such that gene symbols are de-duplicated.</p></td>
<td><p>Yes</p></td>
<td><p><code>results/gene-expression-rsem-fpkm-collapsed.rds</code> included in data download; too large for tracking via GitHub</p></td>
<td><p>N/A</p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>CAVATICA</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/comparative-RNASeq-analysis"><code>comparative-RNASeq-analysis</code></a></p></td>
<td><p><code>gene-expression-rsem-tpm.rds</code> <code>histologies.tsv</code> <code>mend-qc-manifest.tsv</code> <code>mend-qc-results.tar.gz</code></p></td>
<td><p><em>In progress</em>; will produce expression outlier profiles per <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/229">#229</a></p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/compare-gistic"><code>compare-gistic</code></a></p></td>
<td><p><code>cnv-consensus-gistic.zip</code> <code>analyses/run-gistic/results/cnv-consensus-hgat-gistic.zip</code> <code>analyses/run-gistic/results/cnv-consensus-lgat-gistic.zip</code> <code>analyses/run-gistic/results/cnv-consensus-medulloblastoma-gistic.zip</code></p></td>
<td><p>Comparison of the GISTIC results of the entire cohort with the GISTIC results of three individual histolgies, namely, LGAT, HGAT and medulloblastoma <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/547">#547</a></p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/copy_number_consensus_call"><code>copy_number_consensus_call</code></a></p></td>
<td><p><code>cnv-cnvkit.seg.gz</code> <code>cnv-controlfreec.tsv.gz</code> <code>sv-manta.tsv.gz</code></p></td>
<td><p>Produces consensus copy number calls per <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/128">#128</a> and a set of excluded regions where CNV calls are not made</p></td>
<td><p>Yes</p></td>
<td><p><code>results/cnv_consensus.tsv</code> <code>results/cnv-consensus.seg.gz</code> included in data download <code>ref/cnv_excluded_regions.bed</code> <code>ref/cnv_callable.bed</code></p></td>
<td><p><code>focal-cn-file-preparation</code> <code>run-gistic</code></p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>CAVATICA</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/create-subset-files"><code>create-subset-files</code></a></p></td>
<td><p>All files</p></td>
<td><p>This module contains the code to create the subset files used in continuous integration</p></td>
<td><p>No</p></td>
<td><p>All subset files for continuous integration</p></td>
<td><p>All</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>Will set up for OT ticket in</p></td>
</tr>
<tr class="odd">
<td><p>efo-mondo-mapping</p></td>
<td><p><code>histologies.tsv</code></p>
<p><code>efo-mondo-map.tsv</code></p></td>
<td><p>This module contains a file with EFO, MONDO, and NCIT codes for all <code>cancer_group</code> found in histologies.tsv and runs a script to qc in case any <code>cancer_group</code> is missed</p></td>
<td><p>Yes</p></td>
<td><p><code>efo-mondo-mapping.tsv</code></p></td>
<td><p><code>tumor-normal-differential-expression</code></p></td>
<td><p>No</p></td>
<td><p>Yes</p></td>
<td><p>N/A</p></td>
<td><p>Yes</p></td>
</tr>
<tr class="even">
<td><p>filter-mutation-frequencies-tables</p></td>
<td><p><code>gencode.v38.primary_assembly.annotation.gtf.gz</code></p>
<p><code>PMTL_v1.1.tsv</code></p>
<p><code>histologies.tsv</code></p>
<p><code>gene-level-snv-consensus-annotated-mut-freq.tsv.gz</code></p>
<p><code>snv-consensus-plus-hotspots.maf.tsv.gz</code></p>
<p><code>variant-level-snv-consensus-annotated-mut-freq.tsv.gz</code></p>
<p><code>gene-level-cnv-consensus-annotated-mut-freq.tsv.gz</code></p>
<p><code>consensus_wgs_plus_cnvkit_wxs.tsv.gz</code></p>
<p><code>putative-oncogene-fusion-freq.tsv.gz</code></p>
<p><code>fusion-putative-oncogenic.tsv</code></p>
<p><code>putative-oncogene-fused-gene-freq.tsv.gz</code></p></td>
<td><p>Remove <code>Ensembl (ESNG)</code> gene identifier in the OPenPedCan mutation frequency tables, including <a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/snv-frequencies">SNV</a>, <a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/cnv-frequencies">CNV</a> and <a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/fusion-frequencies">fusion</a> that are not in <a href="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/rellease_38/">GENCODE v38 and Ensembl package 104</a>.</p></td>
<td><p>No</p></td>
<td><p>All files from module <code>results</code> directory</p></td>
<td><p>TBD</p></td>
<td><p>No</p></td>
<td><p>Yes</p></td>
<td><p>N/A</p></td>
<td><p>Yes</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/focal-cn-file-preparation"><code>focal-cn-file-preparation</code></a></p></td>
<td><p><code>cnv-cnvkit.seg.gz</code> <code>cnv-controlfreec.tsv.gz</code> <code>gene-expression-rsem-tpm-collapsed.rds</code> <code>cnv-consensus.seg.gz</code></p></td>
<td><p>Maps from copy number variant caller segments to gene identifiers; will be updated to take into account changes that affect entire cytobands, chromosome arms <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/186">#186</a></p></td>
<td><p>Yes</p></td>
<td><p><code>results/cnvkit_annotated_cn_wxs_autosomes.tsv.gz</code> <code>results/cnvkit_annotated_cn_wxs_x_and_y.tsv.gz</code> <code>results/consensus_seg_annotated_cn_autosomes.tsv.gz</code> <code>results/consensus_seg_annotated_cn_x_and_y.tsv.gz</code> <code>results/consensus_wgs_plus_cnvkit_wxs.tsv.gz</code> included in data download <code>results/consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz</code> included in data download <code>results/consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz</code> included in data download</p></td>
<td><p><code>molecular-subtyping-ATRT</code> <code>molecular-subtyping-chordoma</code> <code>molecular-subtyping-embryonal</code> <code>molecular-subtyping-HGG</code> <code>cnv-frequencies</code></p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>CAVATICA</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/fusion_filtering"><code>fusion_filtering</code></a></p></td>
<td><p><code>fusion-arriba.tsv.gz</code> <code>fusion-starfusion.tsv.gz</code> <code>independent-specimens.rnaseq.primary.tsv</code> <code>independent-specimens.rnaseq.relapse.tsv</code></p></td>
<td><p>Standardizes, filters, and prioritizes fusion calls</p></td>
<td><p>Yes</p></td>
<td><p><code>results/fusion-putative-oncogenic.tsv</code>included in data download <code>results/fusion-recurrent-fusion-bycancergroup.tsv</code> <code>results/fusion-recurrent-fusion-bysample.tsv</code> <code>results/fusion-recurrently-fused-genes-bycancergroup.tsv</code> <code>results/fusion-recurrently-fused-genes-bysample.tsv</code></p></td>
<td><p><code>fusion-summary</code> <code>fusion-frequencies</code> <code>molecular-subtyping-HGG</code> <code>molecular-subtyping-LGAT</code> <code>oncoprint-landscape</code> <code>pedot-table-column-display-order-name</code></p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/fusion-summary"><code>fusion-summary</code></a></p></td>
<td><p><code>histologies.tsv</code> <code>fusion-putative-oncogenic.tsv</code> <code>fusion-arriba.tsv.gz</code> <code>fusion-starfusion.tsv.gz</code></p></td>
<td><p>Generate summary tables from fusion files <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/398">#398</a>; <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/623">#623</a></p></td>
<td><p>Yes</p></td>
<td><p><code>results/fusion_summary_embryonal_foi.tsv</code> <code>results/fusion_summary_ependymoma_foi.tsv</code> <code>results/fusion_summary_ewings_foi.tsv</code></p></td>
<td><p><code>molecular-subtyping-EPN</code> <code>molecular-subtyping-EWS</code> <code>molecular-subtyping-embryonal</code></p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/fusion-frequencies"><code>fusion-frequencies</code></a></p></td>
<td><p><code>fusion-putative-oncogenic.tsv</code> <code>independent-specimens.rnaseq.primary.eachcohort.tsv</code> <code>independent-specimens.rnaseq.relapse.eachcohort.tsv</code></p></td>
<td><p>Gather counts and frequencies for fusion per cancer_group and cohort</p></td>
<td><p>No</p></td>
<td><p><code>results/putative-oncogene-fused-gene-freq.jsonl.gz</code> <code>results/putative-oncogene-fused-gene-freq.tsv.gz</code> <code>results/putative-oncogene-fusion-freq.jsonl.gz</code> <code>results/putative-oncogene-fusion-freq.tsv.gz</code></p></td>
<td><p>Upload to FNL BOX</p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p>gene_match</p></td>
<td></td>
<td></td>
<td><p>Yes</p></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/gene-set-enrichment-analysis"><code>gene-set-enrichment-analysis</code></a></p></td>
<td><p><code>gene-expression-rsem-tpm-collapsed.rds</code></p></td>
<td><p>Updated gene set enrichment analysis with appropriate RNA-seq expression data</p></td>
<td><p>No</p></td>
<td><p><code>results/gsva_scores.tsv</code> combined file for all RNA library types</p></td>
<td><p><code>molecular-subtyping-ATRT</code> <code>molecular-subtyping-EPN</code></p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>Move to CAVATICA</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/hotspots-detection"><code>hotspot-detection</code></a></p></td>
<td><p><code>snv-strelka2.vep.maf.gz</code> <code>snv-mutect2.vep.maf.gz</code> <code>snv-vardict.vep.maf.gz</code> <code>snv-lancet.vep.maf.gz</code></p></td>
<td><p>Scavenges cancer any hotspot calls from each caller and merges with consensus 3/3 calls if it was missed in snv-caller workflow.</p></td>
<td><p>No</p></td>
<td><p><code>snv-hotspots-mutation.maf.tsv.gz</code></p></td>
<td><p><code>snv-caller</code></p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>CAVATICA</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/immune-deconv"><code>immune-deconv</code></a></p></td>
<td><p><code>gene-expression-rsem-tpm-collapsed.rds</code></p></td>
<td><p>Immune/Stroma characterization across PBTA part of <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/15">#15</a></p></td>
<td><p>No</p></td>
<td><p><code>results/deconv-output.RData</code></p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/independent-samples"><code>independent-samples</code></a></p></td>
<td><p><code>histologies.tsv</code></p></td>
<td><p>Generates independent specimen lists for WGS/WXS samples</p></td>
<td><p>Yes</p></td>
<td><p><code>results/independent-specimens.wgswxspanel.primary.tsv</code> included in data download <code>results/independent-specimens.wgswxspanel.relapse.tsv</code> included in data download <code>results/independent-specimens.wgswxspanel.primary.eachcohort.tsv</code> included in data download <code>results/independent-specimens.wgswxspanel.relapse.eachcohort.tsv</code> included in data download <code>results/independent-specimens.wgswxspanel.primary.prefer.wxs.tsv</code> included in data download <code>results/independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv</code> included in data download <code>results/independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv</code> included in data download <code>results/independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv</code> included in data download <code>results/independent-specimens.rnaseq.primary.tsv</code> included in data download <code>results/independent-specimens.rnaseq.relapse.tsv</code> included in data download <code>results/independent-specimens.rnaseq.primary.eachcohort.tsv</code> included in data download <code>results/independent-specimens.rnaseq.relapse.eachcohort.tsv</code> included in data download</p></td>
<td><p><code>interaction-plots</code> <code>oncoprint-landscape</code> <code>chromothripsis</code> <code>cnv-frequencies</code> <code>fusion-filtering</code> <code>fusion-frequencies</code> <code>snv-frequencies</code> <code>rna-seq-expression-summary-stats</code></p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/interaction-plots"><code>interaction-plots</code></a></p></td>
<td><p><code>independent-specimens.wgs.primary-plus.tsv</code> <code>snv-consensus-mutation.maf.tsv.gz</code></p></td>
<td><p>Creates interaction plots for mutation mutual exclusivity/co-occurrence <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/13">#13</a>; may be updated to include other data types e.g., fusions</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/long-format-table-utils"><code>long-format-table-utils</code></a></p></td>
<td><p><code>ensg-hugo-rmtl-mapping.tsv</code> <code>analyses/fusion_filtering/references/genelistreference.txt</code> <code>efo-mondo-map.tsv</code> <code>uberon-map-gtex-group.tsv</code> <code>uberon-map-gtex-subgroup.tsv</code></p></td>
<td><p>Functions and scripts for handling long-format tables</p></td>
<td><p>No</p></td>
<td><p><code>annotator/annotation-data/ensg-gene-full-name-refseq-protein.tsv</code> <code>annotator/annotation-data/oncokb-cancer-gene-list.tsv</code></p></td>
<td><p><code>snv-frequencies</code> <code>cnv-frequencies</code> <code>fusion-frequencies</code> <code>rna-seq-expression-summary-stats</code></p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/molecular-subtyping-ATRT"><code>molecular-subtyping-ATRT</code></a></p></td>
<td><p><code>analyses/gene-set-enrichment-analysis/results/gsva_scores.tsv</code> <code>gene-expression-rsem-tpm-collapsed.rds</code> <code>analyses/focal-cn-file-preparation/results/consensus_seg_annotated_cn_autosomes.tsv.gz</code> <code>snv-consensus-mutation-tmb-all.tsv</code> <code>cnv-consensus-gistic.zip</code></p></td>
<td><p><em>Deprecated</em>; Summarizing data into tabular format in order to molecularly subtype ATRT samples <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/244">#244</a>; this analysis did not work</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/molecular-subtyping-CRANIO"><code>molecular-subtyping-CRANIO</code></a></p></td>
<td><p><code>histologies-base.tsv</code> <code>snv-consensus-plus-hotspots.maf.tsv.gz</code></p></td>
<td><p>Molecular subtyping of craniopharyngiomas samples <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/810">#810</a></p></td>
<td><p>No</p></td>
<td><p><code>results/CRANIO_molecular_subtype.tsv</code></p></td>
<td><p><code>molecular-subtyping-pathology</code></p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>Prepare for scaling</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/molecular-subtyping-EPN"><code>molecular-subtyping-EPN</code></a></p></td>
<td><p><code>histologies-base.tsv</code> <code>gene-expression-rsem-tpm-collapsed.rds</code> <code>analyses/chromosomal-instability/breakpoint-data/union_of_breaks_densities.tsv</code> <code>analyses/fusion-summary/results/fusion_summary_ependymoma_foi.tsv</code> <code>analyses/gene-set-enrichment-analysis/results/gsva_scores.tsv</code></p></td>
<td><p>molecular subtyping of ependymoma tumors</p></td>
<td><p>No</p></td>
<td><p><code>results/EPN_all_data_withsubgroup.tsv</code></p></td>
<td><p><code>molecular-subtyping-pathology</code></p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>Will Adapt for OT</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/molecular-subtyping-EWS"><code>molecular-subtyping-EWS</code></a></p></td>
<td><p><code>histologies-base.tsv</code> <code>analyses/fusion-summary/results/fusion_summary_ewings_foi.tsv</code></p></td>
<td><p>Reclassification of tumors based on the presence of defining fusions for Ewing Sarcoma per <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/623">#623</a></p></td>
<td><p>No</p></td>
<td><p><code>results/EWS_samples.tsv</code></p></td>
<td><p><code>molecular-subtyping-pathology</code></p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>Will Adapt for OT</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-HGG"><code>molecular-subtyping-HGG</code></a></p></td>
<td><p><code>histologies-base.tsv</code> <code>snv-consensus-plus-hotspots.maf.tsv.gz</code> <code>consensus_wgs_plus_cnvkit_wxs.tsv.gz</code> <code>fusion-putative-oncogenic.tsv</code> <code>cnv-consensus-gistic.zip</code> <code>gene-expression-rsem-tpm-collapsed.rds</code> <code>tp53_altered_status.tsv</code></p></td>
<td><p>Molecular subtyping of high-grade glioma samples <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/249">#249</a></p></td>
<td><p>No</p></td>
<td><p><code>results/HGG_molecular_subtype.tsv</code></p></td>
<td><p><code>molecular-subtyping-pathology</code></p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-LGAT"><code>molecular-subtyping-LGAT</code></a></p></td>
<td><p><code>histologies-base.tsv</code> <code>snv-consensus-plus-hotspots.maf.tsv.gz</code> <code>fusion-putative-oncogenic.tsv</code> <code>analyses/fusion_filtering/results/fusion-recurrently-fused-genes-bysample.tsv</code></p></td>
<td><p>Molecular subtyping of Low-grade astrocytic tumor samples <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/631">#631</a></p></td>
<td><p>No</p></td>
<td><p><code>results/lgat_subtyping.tsv</code></p></td>
<td><p><code>molecular-subtyping-pathology</code></p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-MB"><code>molecular-subtyping-MB</code></a></p></td>
<td><p><code>histologies.tsv</code> <code>gene-expression-rsem-tpm-collapsed.rds</code></p></td>
<td><p>Molecular classification of Medulloblastoma subtypes part of <a href="https://github.com/PediatricOpenTargets/ticket-tracker/issues/116">#116</a></p></td>
<td><p>No</p></td>
<td><p><code>results/MB_molecular_subtype.tsv</code></p></td>
<td><p><code>molecular-subtyping-pathology</code></p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/molecular-subtyping-SHH-tp53"><code>molecular-subtyping-SHH-tp53</code></a></p></td>
<td><p><code>histologies</code> <code>snv-consensus-plus-hotspots.maf.tsv.gz</code></p></td>
<td><p><em>Deprecated</em>; Identify the SHH-classified medulloblastoma samples that have TP53 mutations <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/247">#247</a></p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/molecular-subtyping-chordoma"><code>molecular-subtyping-chordoma</code></a></p></td>
<td><p><code>analyses/focal-cn-file-preparation/results/consensus_seg_annotated_cn_autosomes.tsv.gz</code> <code>gene-expression-rsem-fpkm-collapsed.stranded.rds</code></p></td>
<td><p>identifying poorly-differentiated chordoma samples per <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/250">#250</a></p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>Will Adapt for OT</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/molecular-subtyping-embryonal"><code>molecular-subtyping-embryonal</code></a></p></td>
<td><p><code>histologies-base.tsv</code> <code>analyses/fusion-summary/fusion_summary_embryonal_foi.tsv</code> <code>sv-manta.tsv.gz</code> <code>consensus_wgs_plus_cnvkit_wxs.tsv.gz</code> <code>analyses/focal-cn-file-preparation/cnvkit_annotated_cn_x_and_y.tsv.gz</code> <code>analyses/focal-cn-file-preparation/controlfreec_annotated_cn_x_and_y.tsv.gz</code> <code>gene-expression-rsem-tpm-collapsed.rds</code></p></td>
<td><p>Molecular subtyping of non-medulloblastoma, non-ATRT embryonal tumors <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/251">#251</a></p></td>
<td><p>No</p></td>
<td><p><code>results/embryonal_tumor_molecular_subtypes.tsv</code></p></td>
<td><p><code>molecular-subtyping-pathology</code></p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>Will Adapt for OT</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-integrate"><code>molecular-subtyping-integrate</code></a></p></td>
<td><p><code>histologies-base.tsv</code> <code>results/compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv</code></p></td>
<td><p>Add molecular subtype information to base histology</p></td>
<td><p>No</p></td>
<td><p><code>results/histologies.tsv</code></p></td>
<td><p>data release</p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/molecular-subtyping-neurocytoma"><code>molecular-subtyping-neurocytoma</code></a></p></td>
<td><p><code>histologies-base.tsv</code></p></td>
<td><p>Molecular subtyping of Neurocytoma samples <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/805">#805</a></p></td>
<td><p>No</p></td>
<td><p><code>results/neurocytoma_subtyping.tsv</code></p></td>
<td><p><code>molecular-subtyping-pathology</code></p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>Will Adapt for OT</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-pathology"><code>molecular-subtyping-pathology</code></a></p></td>
<td><p><code>analyses/molecular-subtyping-CRANIO/results/CRANIO_molecular_subtype.tsv</code> <code>analyses/molecular-subtyping-EPN/results/CRANIO_molecular_subtype.tsv</code> <code>analyses/molecular-subtyping-MB/results/MB_molecular_subtype.tsv</code> <code>analyses/molecular-subtyping-neurocytoma/results/neurocytoma_subtyping.tsv</code> <code>analyses/molecular-subtyping-EWS/results/EWS_samples.tsv</code> <code>analyses/molecular-subtyping-HGG/results/HGG_molecular_subtype.tsv</code> <code>analyses/molecular-subtyping-LGAT/results/lgat_subtyping.tsv</code> <code>analyses/molecular-subtyping-embryonal/results/embryonal_tumor_molecular_subtypes.tsv</code></p></td>
<td><p>Compile output from other molecular subtyping modules and incorporate pathology feedback <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/645">#645</a></p></td>
<td><p>No</p></td>
<td><p><code>results/compiled_molecular_subtyping_with_clinical_feedback.tsv</code> <code>results/compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv</code></p></td>
<td><p><code>molecular-subtyping-integrate</code></p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/mutational-signatures"><code>mutational-signatures</code></a></p></td>
<td><p><code>snv-consensus-plus-hotspots.maf.tsv.gz</code></p></td>
<td><p>Performs COSMIC and Alexandrov et al. mutational signature analysis using the consensus SNV data</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/mutect2-vs-strelka2"><code>mutect2-vs-strelka2</code></a></p></td>
<td><p><code>snv-mutect2.vep.maf.gz</code> <code>snv-strelka2.vep.maf.gz</code></p></td>
<td><p><em>Deprecated</em>; comparison of only two SNV callers, subsumed by <code>snv-callers</code></p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/oncoprint-landscape"><code>oncoprint-landscape</code></a></p></td>
<td><p><code>snv-consensus-plus-hotspots.maf.tsv.gz</code> <code>fusion-putative-oncogenic.tsv</code> <code>analyses/focal-cn-file-preparation/results/controlfreec_annotated_cn_autosomes.tsv.gz</code> <code>independent-specimens.*</code></p></td>
<td><p>Combines mutation, copy number, and fusion data into an OncoPrint plot <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/6">#6</a>; will need to be updated as all data types are refined</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/pedcbio-cnv-prepare"><code>pedcbio-cnv-prepare</code></a></p></td>
<td><p><code>consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz</code> <code>consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz</code></p></td>
<td><p>Generate annotated CNV files that are similar to seg files for PedCBio uploads to include all samples with neutral CNV calls</p></td>
<td><p>Yes</p></td>
<td><p>Upload to PedCBio S3 bucket for ingestion</p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
<td></td>
</tr>
<tr class="odd">
<td><p>pedcbio-sample-name</p></td>
<td></td>
<td></td>
<td><p>No</p></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/pedot-table-column-display-order-name"><code>pedot-table-column-display-order-name</code></a></p></td>
<td><p><code>analyses/snv-frequencies/results/gene-level-snv-consensus-annotated-mut-freq.tsv</code> <code>analyses/snv-frequencies/results/variant-level-snv-consensus-annotated-mut-freq.tsv.gz</code> <code>analyses/cnv-frequencies/results/gene-level-cnv-consensus-annotated-mut-freq.tsv.gz</code> <code>analyses/fusion-frequencies/results/putative-oncogene-fused-gene-freq.tsv.gz</code> <code>analyses/fusion-frequencies/results/putative-oncogene-fusion-freq.tsv.gz</code> <code>analyses/rna-seq-expression-summary-stats/results/long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv.gz</code> <code>analyses/rna-seq-expression-summary-stats/results/long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv.gz</code></p></td>
<td><p>Generate and validate an Excel spreadsheet for Pediatric Open Targets PedOT website table display orders and names</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>Upload to FNL BOX</p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/rna-seq-composition"><code>rna-seq-composition</code></a></p></td>
<td><p><code>gene-expression-rsem-tpm.rds</code> <code>histologies.tsv</code> <code>mend-qc-results.tar.gz</code> <code>mend-qc-manifest.tsv</code> <code>star-log-manifest.tsv</code> <code>star-log-final.tar.gz</code></p></td>
<td><p>Analyzes the fraction of read types that comprise each RNA-Seq sample; flags samples with unusual composition</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/rna-seq-protocol-ruvseq"><code>rna-seq-protocol-ruvseq</code></a></p></td>
<td><p><code>gene-counts-rsem-expected_count-collapsed.rds</code></p></td>
<td><p>Evaluate the use of empirical negative control genes for batch correction</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>Github</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/rna-seq-protocol-dge"><code>rna-seq-protocol-dge</code></a></p></td>
<td><p><code>gene-counts-rsem-expected_count-collapsed.rds</code></p></td>
<td><p><em>In progress</em> <a href="https://github.com/PediatricOpenTargets/ticket-tracker/issues/17">PediatricOpenTargets/ticket-tracker#17</a>; check if the DGE analysis between poly-A and stranded RNA-seq data follow a null-p-value distribution; determine stably expressed genes between poly-A and stranded samples.</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/rna-seq-expression-summary-stats"><code>rna-seq-expression-summary-stats</code></a></p></td>
<td><p><code>gene-expression-rsem-tpm-collapsed.rds</code> <code>histologies.tsv</code></p></td>
<td><p>Calculate TPM summary statistics within each cancer group and cohort. <a href="https://github.com/PediatricOpenTargets/ticket-tracker/issues/51">PediatricOpenTargets/ticket-tracker#51</a>.</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>Upload to FNL Box</p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/run-gistic"><code>run-gistic</code></a></p></td>
<td><p><code>histologies.tsv</code> <code>cnv-consensus.seg.gz</code></p></td>
<td><p>Runs GISTIC 2.0 on SEG files</p></td>
<td><p>Yes</p></td>
<td><p><code>cnv-consensus-gistic.zip</code> included in data download</p></td>
<td><p><code>cnv-chrom-plot</code> <code>compare-gistic</code> <code>molecular-subtyping-ATRT</code> <code>molecular-subtyping-HGG</code></p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>Move to CAVATICA</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/sample-distribution-analysis"><code>sample-distribution-analysis</code></a></p></td>
<td><p><code>histologies.tsv</code></p></td>
<td><p>Produces plots and tables that illustrate the distribution of different histologies in the PBTA data</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/selection-strategy-comparison"><code>selection-strategy-comparison</code></a></p></td>
<td><p><code>gene-expression-rsem-tpm-collapsed.rds</code></p></td>
<td><p><em>Deprecated</em>; Comparison of RNA-seq data from different selection strategies</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/sex-prediction-from-RNASeq"><code>sex-prediction-from-RNASeq</code></a></p></td>
<td><p><code>gene-expression-kallisto.stranded.rds</code> <code>histologies.tsv</code></p></td>
<td><p>predicts genetic sex using RNA-seq data <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/6">#84</a></p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers"><code>snv-callers</code></a></p></td>
<td><p><code>snv-lancet.vep.maf.gz</code> <code>snv-mutect2.vep.maf.gz</code> <code>snv-strelka2.vep.maf.gz</code> <code>snv-vardict.vep.maf.gz</code> <code>tcga-snv-lancet.vep.maf.gz</code> <code>tcga-snv-mutect2.vep.maf.gz</code> <code>tcga-snv-strelka2.vep.maf.gz</code></p></td>
<td><p>Generates consensus SNV and indel calls for PBTA and TCGA data; calculates tumor mutation burden using the consensus calls</p></td>
<td><p>Yes</p></td>
<td><p><code>results/consensus/snv-consensus-plus-hotspots.maf.tsv</code> included in data download; too large for tracking via GitHub <code>results/consensus/snv-consensus-mutation-tmb-all.tsv</code> <code>results/consensus/snv-consensus-mutation-tmb-coding.tsv</code>included in data download; too large for tracking via GitHub <code>results/consensus/tcga-snv-consensus-mutation.maf.tsv.gz</code> <code>results/consensus/tcga-snv-mutation-tmb.tsv</code> <code>results/consensus/tcga-snv-mutation-tmb-coding.tsv</code></p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/snv-frequencies"><code>snv-frequencies</code></a></p></td>
<td><p><code>histologies.tsv</code> <code>snv-consensus-plus-hotspots.maf.tsv.gz</code> <code>independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv</code> <code>independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv</code> <code>independent-specimens.wgswxspanel.primary.prefer.wxs.tsv</code> <code>independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv</code></p></td>
<td><p>Annotate SNV table with mutation frequencies</p></td>
<td><p>No</p></td>
<td><p><code>results/gene-level-snv-consensus-annotated-mut-freq.jsonl.gz</code> <code>results/gene-level-snv-consensus-annotated-mut-freq.tsv.gz</code> <code>variant-level-snv-consensus-annotated-mut-freq.jsonl.gz</code> <code>variant-level-snv-consensus-annotated-mut-freq.tsv.gz</code></p></td>
<td><p>Upload to FNL BOX</p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/snv-frequencies"><code>oncoprint-landscape</code></a></p></td>
<td><p><code>histologies.tsv</code> <code>snv-consensus-plus-hotspots.maf.tsv.gz</code> <code>independent-specimens.wgswxspanel.primary.eachcohort.prefer.wxs.tsv</code> <code>independent-specimens.wgswxspanel.relapse.eachcohort.prefer.wxs.tsv</code> <code>independent-specimens.wgswxspanel.primary.prefer.wxs.tsv</code> <code>independent-specimens.wgswxspanel.relapse.prefer.wxs.tsv</code></p></td>
<td><p>Annotate SNV table with mutation frequencies</p></td>
<td><p>No</p></td>
<td><p><code>oncoprint-landscape</code> <code>tp53_nf1_score</code> <code>molecular-subtyping-CRANIO</code> <code>molecular-subtyping-HGG</code> <code>molecular-subtyping-LGAT</code> <code>molecular-subtyping-SHH-tp53</code> <code>mutational-signatures</code></p></td>
<td><p>N/A</p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/ssgsea-hallmark"><code>ssgsea-hallmark</code></a></p></td>
<td><p><code>gene-counts-rsem-expected_count-collapsed.rds</code></p></td>
<td><p><em>Deprecated</em>; performs GSVA using Hallmark gene sets</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/survival-analysis"><code>survival-analysis</code></a></p></td>
<td><p>TBD</p></td>
<td><p><em>In progress</em>; will eventually contain functions for various types of survival analysis <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/18">#18</a></p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/telomerase-activity-prediction"><code>telomerase-activity-prediction</code></a></p></td>
<td><p><code>gene-expression-rsem-tpm-collapsed.rds</code> <code>gene-counts-rsem-expected_count-collapsed.rds</code></p></td>
<td><p>Quantify telomerase activity across pediatric brain tumors part of <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/148">#148</a></p></td>
<td><p>No</p></td>
<td><p><code>results/TelomeraseScores_PTBAPolya_counts</code> <code>results/TelomeraseScores_PTBAPolya_FPKM.txt</code> <code>results/TelomeraseScores_PTBAStranded_counts.txt</code> <code>results/TelomeraseScores_PTBAStranded_FPKM.txt</code></p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p>tmb-calculation</p></td>
<td><p><code>gencode.v27.primary_assembly.annotation.bed</code></p>
<p><code>intersect_strelka_mutect2_vardict_WGS.bed</code></p>
<p><code>snv-consensus-plus-hotspots.maf.tsv.gz</code></p>
<p><code>biospecimen_id_to_bed_map.tsv</code></p>
<p><code>histologies-base.tsv</code></p>
<p><code>hg38_strelka.bed</code></p>
<p><code>wgs_canonical_calling_regions.hg38.bed</code></p>
<p><code>gencode.v27.primary_assembly.annotation.gtf.gz</code></p></td>
<td><p>This analysis utilizes the SNV consensus MAF file, <code>../../data/snv-consensus-plus-hotspots.maf.tsv.gz</code> from <a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis">Pediatric Open Targets, OPenPedCan-analysis</a> datasets to calculate Tumor Mutation Burden TMB for each experimental strategy WGS, WXS, and Targeted Sequencing tumor sample with SNV calls for all cohorts and cancer types evaluated in the project. The Tumor Mutation Burden calculation is adapted from <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers">snv-callers module</a> of the OpenPBTA-analyses, and use the SNV calls Mutect2, Strelka2, Lancet, and Vardict callers.</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/tmb-compare"><code>tmb-compare</code></a></p></td>
<td><p><code>snv-consensus-mutation-tmb-coding.tsv</code></p></td>
<td><p>Compares PBTA tumor mutation burden to adult TCGA data. The D3B TMB calculations <code>TMB_d3b_code</code> and its comparison notebook <code>compare-tmb-calculations.Rmd</code> are <em>deprecated</em>.</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/PediatricOpenTargets/OpenPedCan-analysis/tree/dev/analyses/tp53_nf1_score"><code>tp53_nf1_score</code></a></p></td>
<td><p><code>snv-consensus-plus-hotspots.maf.tsv</code> <code>gene-expression-rsem-tpm-collapsed.rds</code> <code>consensus_wgs_plus_cnvkit_wxs.tsv.gz</code></p></td>
<td><p>Applies <em>TP53</em> inactivation, <em>NF1</em> inactivation, and Ras activation classifiers to RNA-seq data <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/165">#165</a></p></td>
<td><p>No</p></td>
<td><p><code>TP53_NF1_snv_alteration.tsv</code> <code>gene-expression-rsem-tpm-collapsed_classifier_scores.tsv</code> <code>loss_overlap_domains_tp53.tsv</code> <code>poly-A_TP53.png</code> <code>stranded_TP53.png</code> <code>sv_overlap_tp53.tsv</code> <code>tp53_altered_status.tsv</code></p></td>
<td><p><code>molecular-subtyping-HGG</code></p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/transcriptomic-dimension-reduction"><code>transcriptomic-dimension-reduction</code></a></p></td>
<td><p><code>gene-expression-rsem-tpm.rds</code> <code>gene-expression-kallisto.rds</code></p></td>
<td><p>Dimension reduction and visualization of RNA-seq data part of <a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/9">#9</a></p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>No</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p><a href="https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/tcga-capture-kit-investigation"><code>tcga-capture-kit-investigation</code></a></p></td>
<td><p><code>snv-lancet.vep.maf.gz</code> <code>snv-mutect2.vep.maf.gz</code> <code>snv-strelka2.vep.maf.gz</code> <code>tcga-snv-lancet.vep.maf.gz</code> <code>tcga-snv-mutect2.vep.maf.gz</code> <code>tcga-snv-strelka2.vep.maf.gz</code> <code>histologies.tsv</code> <code>tcga-manifest.tsv</code> <code>WGS.hg38.lancet.unpadded.bed</code> <code>WGS.hg38.strelka2.unpadded.bed</code> <code>WGS.hg38.mutect2.vardict.unpadded.bed</code></p></td>
<td><p>Investigation of the TMB discrepancy between PBTA and TCGA data</p></td>
<td><p>No</p></td>
<td><p><code>results/*.bed</code></p></td>
<td><p>N/A</p></td>
<td><p>Yes</p></td>
<td><p>No</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="even">
<td><p><code>tumor-gtex-plots</code></p></td>
<td><p><code>gene-expression-rsem-tpm-collapsed.rds</code> <code>histologies.tsv</code></p></td>
<td><p><em>In progress</em> <a href="https://github.com/PediatricOpenTargets/ticket-tracker/issues/38">PediatricOpenTargets/ticket-tracker#38</a>; tumor vs normal and tumor only expression plots</p></td>
<td><p>No</p></td>
<td><p><code>results/pan_cancer_plots_cancer_group_level.{tsv, jsonl.gz}</code> <code>results/pan_cancer_plots_cohort_cancer_group_level.{tsv, jsonl.gz}</code> <code>results/tumor_normal_gtex_plots_cancer_group_level.{tsv, jsonl.gz}</code> <code>results/tumor_normal_gtex_plots_cohort_cancer_group_level.{tsv, jsonl.gz}</code> <code>results/metadata.tsv</code> <code>plots/*.png</code></p></td>
<td><p>N/A</p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>GitHub</p></td>
<td><p>N/A</p></td>
</tr>
<tr class="odd">
<td><p>tumor-normal-differential-expression</p></td>
<td><p><code>histologies.tsv</code></p>
<p><code>gene-counts-rsem-expected_count-collapsed.rds</code></p>
<p><code>independent-specimens.rnaseq.primary.tsv</code></p>
<p><code>independent-specimens.rnaseq.primary.eachcohort.tsv</code></p>
<p><code>gene-expression-rsem-tpm-collapsed.rds</code></p>
<p><code>ensg-hugo-pmtl-mapping.tsv</code></p>
<p><code>efo-mondo-map.tsv</code></p>
<p><code>uberon-map-gtex-subgroup.tsv</code></p></td>
<td><p>This module takes as input histologies and the RNA-Seq expression matrices data, and performs differential expression analysis for all combinations of GTEx subgroup normal and cancer histology type tumor.</p></td>
<td><p>No</p></td>
<td><p>N/A</p></td>
<td><p>N/A</p></td>
<td><p>Yes</p></td>
<td><p>Yes</p></td>
<td><p>HPC</p>
<p>CAVATICA user can create application for personal analysis purpose using scripts provided in the module</p></td>
<td><p>N/A</p></td>
</tr>
</tbody>
</table>
