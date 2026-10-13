---
title: "Wrap up"
teaching: 10
exercises: 0
---

:::::::::::::::::::::::::::::::::::::: questions

- What biological questions can we now answer with the tools we learnt today?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Connect the methods covered in this workshop to the biological questions they can answer
- Identify the next steps and resources for your own analyses

::::::::::::::::::::::::::::::::::::::::::::::::

## What we covered

| Biological question | Approach | Key functions |
|---|---|---|
| Who is in my samples, and how confident am I? | Compare profilers; check agreement between tools | nf-core/taxprofiler, Taxpasta |
| How many species are there, and how even is the community? | Alpha diversity | `mia::addAlpha()`, `scater::plotColData()`, `wilcox.test()` |
| Do communities differ between groups? | Beta diversity, PCoA, PERMANOVA | `mia::getDissimilarity()`, `mia::addMDS()`, `vegan::adonis2()`, `vegan::betadisper()` |
| Which taxa differ between groups? | Differential abundance | `ANCOMBC::ancombc2()` |

## Lessons for your own analyses

- **Choose the profiler and abundance type deliberately**, based on your question, and report the tool, version, database and abundance type.
- **Treat the profiler as a source of variation.** Don't pool data processed with different tools or databases without accounting for it.
- **Use agreement between profilers as supporting evidence**, especially for rare taxa.
- **Pair statistics with plots**, and report effect sizes with confidence intervals.

## Where to go next

- [Orchestrating Microbiome Analysis with Bioconductor (OMA)](https://microbiome.github.io/OMA/): the online book behind `mia` and `miaViz`
- [nf-core/taxprofiler](https://nf-co.re/taxprofiler) and [nf-core/funcprofiler](https://nf-co.re/funcprofiler) documentation
- Sun *et al.* (2021), [*Challenges in benchmarking metagenomic profilers*](https://doi.org/10.1038/s41592-021-01141-3), *Nature Methods*
- [Melbourne Bioinformatics training](https://www.eventbrite.com.au/d/australia--melbourne/melbourne-bioinformatics/): upcoming workshops

<!-- FIXME: add functional profiling follow-up (HUMAnN / eggNOG outputs from Day 1) if covered in a later session -->

::::::::::::::::::::::::::::::::::::: discussion

### Q&A

What questions do you have about applying these methods to your own data?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: keypoints

- Profiling, diversity analysis and differential abundance testing each answer a different biological question.
- The choice of profiler, database and abundance type shapes every downstream result.

::::::::::::::::::::::::::::::::::::::::::::::::
