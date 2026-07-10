Analysis code and aggregate result tables for a transcriptomic study of primary DCIS,
identifying co-expression network and isoform-usage features in the primary lesion that
distinguish patients who subsequently develop an ipsilateral recurrence from those who do not.

## Data availability
Raw and individual-level data are **controlled-access** (UK Sloane Project) and are **not**
included in this repository. Only analysis code and aggregate, gene-level summary tables are provided.

## Repository structure
- `scripts/`
  - `01_WGCNA_pipeline_HPC.Rmd` — primary WGCNA (batch-corrected matrix; 98 QC-passing samples)
  - `01b_WGCNA_uncorrected_HPC.Rmd` — sensitivity: uncorrected matrix, batch modelled at module–trait testing
  - `02_HER2_validation_HPC.Rmd` — HER2 module–trait ROC validation
  - `03_paired_recurrence_HPC.Rmd` — paired primary-vs-recurrence DE (DESeq2)
  - `04b_isoform_switch_QC98.Rmd` — isoform switching (satuRn), 98 QC-passing samples
  - `04c_qc98_compare.R`, `04e_candidate_switchplots.R` — isoform helpers / candidate switch plots
  - `05_network_visualisation_HPC.R` — co-expression network figures
  - `config.R` — paths and helper functions
  - `run_in_container.sh`, `run_01b.sh`, `dcis_isoform.def` — run wrappers and container definition
- `results/` — aggregate result tables (WGCNA, uncorrected-WGCNA sensitivity, enrichment, QC-98 isoform, paired and case/control DEGs)

## Environment
- WGCNA / paired pipelines: R 4.5.1 / Bioconductor 3.22 (KCL CREATE HPC)
- Isoform pipeline: Singularity container, R 4.4.1 / Bioconductor 3.19 (IsoformSwitchAnalyzeR 2.4.0, satuRn 1.12.0)
- Full session details: `scripts/sessionInfo_HPC.txt`, `scripts/sessionInfo_isoform_container.txt`
