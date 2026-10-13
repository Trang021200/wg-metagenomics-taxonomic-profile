---
title: Setup
---

## Tutorial overview

This is **Day 2** of the Melbourne Bioinformatics WGS metagenomics workshop series. Day 2 is about analysing the taxonomic profiling outputs from **nf-core/taxprofiler** in R. On Day 1, we went from raw shotgun metagenomic reads to standardised abundance tables. Today we import those tables into R. We will compute alpha and beta diversity based on the output of each profiler, and then compare the statistics among profilers against the ground truth.

We will begin with a short recap of Day 1. This is followed by the concepts behind four taxonomic profilers (MetaPhlAn, mOTUs, Kraken2/Bracken and sylph): the reference databases they use, how they assign reads, and the types of abundance they report. Participants will then compare the profilers' outputs against the known composition of the samples, and see how the choice of profiler affects alpha and beta diversity. Finally, we will introduce differential abundance testing to identify taxa that differ between groups of samples.


::::callout
# LEARNING OBJECTIVES

By the end of this tutorial, you will be able to:

- Explain how reference databases, profiling methods and abundance types cause taxonomic profilers to give different results.
- Choose an appropriate profiler and output for a research question, and use the other outputs as complementary information.
- Compare taxonomic profiles against a ground truth using distance metrics and precision–recall analysis.
- Calculate and compare alpha and beta diversity in R using `mia`, `scater` and `vegan`.
- Perform and interpret a differential abundance analysis with ANCOM-BC2.

::::

<br>

### Timing

The anticipated workshop duration when delivered to a group of participants is **4 hours**, plus breaks.

<br>

### Skill level

This intermediate workshop assumes a basic knowledge of R and familiarity with the content of Day 1. No background in statistics beyond basic hypothesis testing (p-values) is assumed. We introduce the diversity and distance metrics as we go.

## Prior to the workshop

:::: checklist

## REQUIRED KNOWLEDGE

- This workshop assumes participants have a **basic understanding of R**, including variables, data frames, installing and loading packages, and running code in RStudio, or have previously attended an Introduction to R workshop. Please review introductory materials [here](https://melbournebioinformatics.github.io/intro-to-r/).
- Familiarity with the content of **Day 1**: *Taxonomic and functional profiling from shotgun metagenomics data with nf-core*. We will recap the key points at the start of the session, but you will get the most out of Day 2 if you know what the nf-core/taxprofiler outputs are.
- Installation of R (version 4.6.0 or later) and RStudio
::::



:::: checklist

## REQUIRED SOFTWARE

Attendees are required to bring their own laptop computers (Windows, macOS or Linux) with at least 8 GB of RAM and about 5 GB of free disk space. **Please ensure you have installed:**

- [Chrome](https://www.google.com/chrome/) or [Firefox](https://www.mozilla.org/en-US/firefox/)
- [R](https://cran.ms.unimelb.edu.au/), **version 4.6.0 or later** (download and install the latest version of R using the UniMelb mirror)
- [RStudio Desktop](https://posit.co/download/rstudio-desktop/) (install R **before** RStudio, so RStudio can find it)
- **R packages** required for this workshop (see below)

### Installing required R packages

Please copy and run the code below in the RStudio **Console** before the workshop. Installing the Bioconductor packages can take 15–30 minutes. When asked *"Update all/some/none? [a/s/n]"*, type `a` and press Enter. If asked whether to install from source, type `no`.

```r
# Packages from CRAN
cran_packages <- c(
  "dplyr",       # Data manipulation
  "tidyr",       # Reshaping tables
  "readr",       # Reading tables
  "ggplot2",     # Plotting
  "vegan"        # PERMANOVA (adonis2) and dispersion tests (betadisper)
)

# Packages from Bioconductor
bioc_packages <- c(
  "TreeSummarizedExperiment",  # Data container for microbiome data
  "SummarizedExperiment",      # Accessing assays and sample metadata
  "mia",                       # Diversity, transformations, ordination, PERMANOVA
  "miaViz",                    # Microbiome plotting
  "scater",                    # plotColData(), plotReducedDim()
  "ANCOMBC"                    # Differential abundance testing
)

# 1. Install CRAN packages
install.packages(cran_packages)

# 2. Install Bioconductor manager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 3. Install Bioconductor packages
BiocManager::install(bioc_packages)
```

### Checking your installation

Run the code below. Every package should print `TRUE`.

```r
sapply(c(cran_packages, bioc_packages), requireNamespace, quietly = TRUE)

R.version.string     # should be R version 4.6.x or later
BiocManager::version()
```

::::

:::::::::::::::: spoiler

### Why R 4.6?

The current Bioconductor release (3.23) requires R 4.6. Most of the packages we use (`mia`, `miaViz`, `ANCOMBC`, `scater`) come from Bioconductor. On an older version of R, `BiocManager` will install outdated versions of these packages, or will not install them at all. Check your version by typing `R.version.string` in the R console.

::::::::::::::::::::::::


::::::::::::::::::::::::::::::::::::::: caution

### Having trouble?

If any package returns `FALSE`, try installing it again on its own, e.g. `BiocManager::install("ANCOMBC")`, and read the error message carefully. The most common causes are an outdated R version and, on Linux, missing system libraries. If you are still stuck, email us at [bioinformatics-training@unimelb.edu.au](mailto:bioinformatics-training@unimelb.edu.au) with a screenshot of the error **before** the workshop.

:::::::::::::::::::::::::::::::::::::::::::::::

<br>

:::: checklist

## REQUIRED DATA

<!--
FIXME: update the link to the data once setup is finalised
HOW WILL THE OUTPUT FROM DAY1 GET DOWNLOADED?
-->

**Please download the data required for this workshop:**

1. Create a new folder for the workshop, e.g. `metagenomics_day2` on your Desktop.
2. Download the data and unzip it into that folder.
3. In RStudio, go to **File → New Project → Existing Directory** and select the folder. This sets your working directory so the file paths in the lesson work as written.

We will use the outputs of **nf-core/taxprofiler** and **nf-core/funcprofiler** generated on Day 1, which reflect the taxonomic profiles of a subset of the original data. We also provide the full taxonomic profiles, which were generated using the same **nf-core** pipelines used in Day 1. All statistical inferences and biological interpretations in this workshop will be done for the complete profiles.

You do **not** need to re-run the Day 1 pipelines. We provide the processed tables so that everyone starts Day 2 from the same point.

The data folder contains:

| File | Description |
|---|---|
| `taxpasta_*.tsv` | Standardised, merged abundance tables from Taxpasta, one for each profiler (MetaPhlAn, mOTUs, Kraken2/Bracken, sylph) |
| `*_genefamilies.tsv` | Gene family abundances in reads per kilobase (RPK), stratified by contributing species |
| `*_pathabundance.tsv` | Metabolic pathway abundances in RPK, stratified by species contribution |
| `*_pathcoverage.tsv` | Pathway coverage scores (0–1), indicating the fraction of reactions detected per pathway |
| `*_profile.txt` | Species-level taxonomic abundance profile used as input to HUMAnN |
| `ground_truth.tsv` | Known composition of each sample, used as the reference for the distance and precision/recall analyses |
| `metadata.tsv` | Sample metadata |

The data come from the benchmarking study by Sun *et al.* (2021), [*Challenges in benchmarking metagenomic profilers*](https://doi.org/10.1038/s41592-021-01141-3), *Nature Methods* 18, 618–626. The raw data are available from the authors' [Figshare project](https://figshare.com/projects/Pitfalls_and_Opportunities_in_Benchmarking_Metagenomic_Classifiers/79916).

::::

::: discussion

# Attended Day 1?

If you ran nf-core/taxprofiler yourself on Day 1, you are welcome to use your own outputs. They are in `results/taxpasta/` in your pipeline output directory. To make it easier to follow along with the instructors, we recommend using the provided files during the workshop.

:::
