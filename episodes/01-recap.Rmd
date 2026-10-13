---
title: "Recap from Day 1"
teaching: 10
exercises: 0
---

:::::::::::::::::::::::::::::::::::::: questions

- What did we produce on Day 1, and where does Day 2 pick up?
- Which `nf-core/taxprofiler` and `nf-core/funcprofiler` outputs will we analyse today?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Review the main stages from FASTQ files to nf-core pipeline outputs
- Identify the files from Day 1 that are the inputs for today's analysis

::::::::::::::::::::::::::::::::::::::::::::::::

## From reads to abundance tables

On Day 1 we went from raw shotgun metagenomic reads to standardised abundance tables:

| Stage | What happened | Tool | Output |
|---|---|---|---|
| 1. Data | Downloaded paired-end FASTQ files and sample metadata | ENA browser, `wget`, MD5 check | Raw reads, `metadata.tsv` |
| 2. Pipeline setup | Built sample and database sheets and a reproducible launch command | Nextflow, nf-core | `samplesheet.csv`, `databases.csv`, `params.json` |
| 3. Pre-processing | Read QC, adapter trimming, host-read removal | nf-core/taxprofiler (fastp, Bowtie2, ...) | MultiQC report |
| 4. Taxonomic profiling | The same reads were classified by four profilers | MetaPhlAn, mOTUs, Kraken2/Bracken, sylph | One profile per tool and sample |
| 5. Standardisation | Profiles were converted to a common format and merged | Taxpasta | One merged table per profiler |
| 6. Import into R | Tables and metadata were loaded into a single object | `mia`, `TreeSummarizedExperiment` | `TreeSummarizedExperiment` object |

<!-- FIXME: replace with a workflow figure (e.g. episodes/fig/workflow.png) shared with the Day 1 lesson -->

## The dataset

We are using samples from the benchmarking study by Sun *et al.* (2021), [*Challenges in benchmarking metagenomic profilers*](https://doi.org/10.1038/s41592-021-01141-3). These are mock and simulated communities, so we know the **true** composition of each sample. On Day 1, we used the a subset of the original data as input for the pipelines, which does not reflect the complete taxonomic and functional profiles of the samples. Therefore, we also provide you with the outputs from the full metagen

<!-- FIXME: add the number of samples, community types and sequencing depth once the subset is fixed -->

## What we do today

1. **Concepts:** why the four profilers give different answers for the same reads.
2. **Practice:** compare the profilers against the ground truth, and see how the choice of profiler changes alpha and beta diversity.
3. **Differential abundance:** test which taxa differ between groups of samples.

::::::::::::::::::::::::::::::::::::: discussion

### Quick check

On Day 1, every sample was run through four profilers. Before we look at the results, do you expect them to agree? On what: the list of species present, their abundances, or both?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: keypoints

- nf-core/taxprofiler ran four profilers on the same reads, and Taxpasta standardised their outputs into comparable tables.
- Our samples have a known composition (the ground truth), so we can judge how accurate each profiler is.
- Today's analysis starts from the merged Taxpasta tables and the sample metadata.

::::::::::::::::::::::::::::::::::::::::::::::::
