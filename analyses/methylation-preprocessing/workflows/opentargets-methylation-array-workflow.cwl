cwlVersion: v1.2
class: Workflow
id: opentargets-methylation-array-workflow
label: OpenTargets Methylation Array Workflow
doc: |
  Preprocess probe hybridization intensity values of selected methylated and
  unmethylated cytosine (CpG) loci into usable methylation measurements and copy number for the
  [Pediatric Open Targets, OPenPedCan-analysis](https://github.com/PediatricOpenTargets/OpenPedCan-analysis)
  raw DNA methylation array datasets.

  This workflow is used to efficiently handle a large number of IDAT files. The workflow
  will:
  1. run the preprocess-illumina-arrays.R module on a given directory with IDAT files
  2. returns methylation beta-values, m-values, and cn-values matrices, and array probe annotation table
requirements:
- class: InlineJavascriptRequirement
- class: MultipleInputFeatureRequirement
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
inputs:
  # Required Inputs
  input_idats_dir: { type: Directory, doc: "Directory containing the IDATs to process." }
  output_basename: { type: 'string', doc: "String to use as basename for output filenames." }

  # Optional Arguments
  controls_present: { type: 'boolean?', doc: "If set, preprocesses the Illumina methylation array dataset assuming presence of either normal and tumor samples or samples of mutiple cancer groups or both." }
  snp_filter: { type: 'boolean?', doc: "If set, drops the probes that contain either a SNP at the CpG interrogation or at the single nucleotide extension." }

  # Resource Control
  preprocess_methylation_array_cores: { type: 'int?', doc: "Cores to give to the preprocess_methylation_arrays step" }
  preprocess_methylation_array_ram: { type: 'int?', doc: "GB of RAM to give to the preprocess_methylation_arrays step" }
outputs:
  methylation_beta_values: { type: File, outputSource: preprocess_methylation_arrays/beta_values, doc: "" }
  methylation_m_values: { type: File, outputSource: preprocess_methylation_arrays/m_values, doc: "" }
  methylation_cn_values: { type: File, outputSource: preprocess_methylation_arrays/cn_values, doc: "" }

steps:
  preprocess_methylation_arrays:
    run: ../tools/preprocess_illumina_arrays.cwl
    in:
      input_idats_dir: input_idats_dir
      controls_present: controls_present
      snp_filter: snp_filter
      cores: preprocess_methylation_array_cores
      ram: preprocess_methylation_array_ram
    out: [beta_values, m_values, cn_values]

$namespaces:
  sbg: https://sevenbridges.com
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 2
