---
title: "Comparing taxonomic profilers: concepts"
teaching: 35
exercises: 10
---

:::::::::::::::::::::::::::::::::::::: questions

- Why do different taxonomic profilers give different answers for the same reads?
- What does "abundance" actually mean in each profiler's output?
- How to choose the right profiler and database for the research question?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Describe the two main reference taxonomies used in metagenomics: NCBI RefSeq and GTDB
- Distinguish between sequence abundance, relative taxonomic abundance, and absolute abundance
- Explain Average Nucleotide Identity (ANI) and containment ANI
- Describe the database and approach of MetaPhlAn, mOTUs, Kraken2/Bracken and sylph
- Weigh the pros and cons of each approach for a given scenario
- Choose an appropriate output for a research question, and use the other outputs as complementary information

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

A taxonomic profiler answers two questions: *who is there?* and *how much of each?* Each profiler's answer depends on three choices:

1. the **reference database** it compares reads against;
2. the **method** it uses to assign reads, or parts of reads, to taxa;
3. the **type of abundance** it reports.

Two profilers that differ in any of these can give different results for the same sample, and both can be "correct". Before we compare the outputs in R (next episode), we need to understand these differences.

## Reference databases: the two pillars

Most metagenomic databases are built on one of two taxonomies:

| | NCBI Taxonomy / RefSeq | GTDB |
|---|---|---|
| What it is | NCBI's curated collection of reference sequences, labelled with the NCBI Taxonomy | The Genome Taxonomy Database: a taxonomy of bacteria and archaea built from genome phylogeny |
| How taxa are defined | Historical, literature-based classification. Ranks are not standardised across lineages | Phylogeny of conserved marker proteins, with ranks normalised by relative evolutionary divergence. Species are clusters of genomes at about 95% ANI |
| Coverage | Bacteria, archaea, viruses, fungi, other eukaryotes, plasmids, the human genome | Bacteria and archaea only, including many uncultured genomes (MAGs) |
| Names | Familiar, e.g. *Shigella flexneri*, phylum Firmicutes / Bacillota | Can differ from NCBI. For example, *Shigella* species fall within *Escherichia coli*, and genera are split with suffixes such as `Clostridium_A` |

<!-- FIXME: double-check the examples and add the GTDB release used by the workshop's sylph database -->

::::::::::::::::::::::::::::::::::::: callout

### Names are not interchangeable

A profiler using a GTDB database and one using an NCBI database can detect the same organism under **different names**. Before comparing the profilers' outputs, we have to map them to a common taxonomy. Otherwise a correct detection looks like a false positive plus a false negative.

::::::::::::::::::::::::::::::::::::::::::::::::

## What does "abundance" mean?

Sun *et al.* (2021) showed that much of the disagreement between profilers comes from reporting different kinds of abundance:

- **Sequence abundance**: the fraction of *sequenced reads (or bases)* that come from a taxon.
- **Taxonomic abundance**: the fraction of *cells (genomes)* that belong to a taxon.
- **Absolute abundance**: the actual amount of a taxon, e.g. cells per gram of stool. None of the profilers we use can estimate this from sequencing alone. It requires extra measurements such as spike-ins, qPCR or flow cytometry.

Sequence and taxonomic abundance are related through the genome size. A taxon with a genome twice as large contributes roughly twice as many reads per cell:

$$\text{sequence abundance}_i \propto \text{taxonomic abundance}_i \times \text{genome size}_i$$

::::::::::::::::::::::::::::::::::::: challenge

### Challenge 1: Same cells, different numbers

A community has equal numbers of cells of species A (genome size 2 Mb) and species B (genome size 6 Mb). Assuming every read is classified, what is each species' taxonomic abundance and sequence abundance?

:::::::::::::::::::::::: solution

Taxonomic abundance: A = 50%, B = 50%.

Sequence abundance: B contributes three times as many bases per cell as A, so A = 2 / (2 + 6) = 25% and B = 75%.

The same community gives very different numbers depending on what is reported.

:::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::

## ANI and containment ANI

**Average nucleotide identity (ANI)** measures how similar two genomes are. Genomes with roughly ≥ 95% ANI are usually considered the same species, which is the species definition GTDB uses.

**Containment ANI** estimates the ANI between a *reference genome* and the *reads in a metagenome*. It is based on the fraction of the reference genome's k-mers found (contained) in the sample. If the reference is present, nearly all its k-mers should be found. A lower containment ANI suggests a related strain or species rather than the reference itself. At low sequencing depth, many k-mers are missed simply by chance. sylph corrects for this and reports an **adjusted ANI**.

## The four profilers

| | MetaPhlAn | mOTUs | Kraken2 + Bracken | sylph |
|---|---|---|---|---|
| Approach | Maps reads to **clade-specific marker genes** | Maps reads to **universal single-copy marker genes** (10 genes) | **Exact k-mer matching**, each read assigned to the lowest common ancestor (LCA) of matching taxa. Bracken then re-estimates abundance at a chosen rank | **k-mer sketching**. Estimates containment ANI and coverage for each reference genome |
| Database | ChocoPhlAn marker database, organised into species-level genome bins (SGBs) | mOTUs database: marker-gene clusters from reference genomes *and* metagenome-assembled genomes | User-chosen, often built from RefSeq (e.g. the Kraken2 *Standard* database). Bracken files are built per database and read length | Pre-built genome databases, e.g. GTDB |
| Reports | Taxonomic abundance (relative) | Relative abundance from single-copy markers (≈ taxonomic). Read counts optional | Reads per taxon (sequence abundance) | **Both** sequence and taxonomic abundance, plus ANI and coverage |
| Unknown taxa | Unknown SGBs (uSGBs) from MAGs | Reports an *unassigned* (`-1`) fraction for species not in the database | Reads stay unclassified, or are assigned to a higher rank | Genomes below the ANI threshold are not reported |

<!-- FIXME: confirm tool/database versions used on Day 1 and which abundance column Taxpasta keeps for each profiler -->

### Pros and cons

- **MetaPhlAn:** fast, with few false positives thanks to the unique marker genes. Covers bacteria, archaea and some eukaryotes and viruses. Can miss taxa that have few markers or are present at low coverage. Gives no read counts by default.
- **mOTUs:** estimates the proportion of *unknown* species, and works well for poorly characterised environments thanks to its MAG-derived markers. Only covers prokaryotes.
- **Kraken2/Bracken:** classifies every read, can include viruses and eukaryotes, and its read counts suit count-based statistics. Prone to false positives, especially with small or incomplete databases, so it often needs abundance or read-count filtering. Its abundance is sequence abundance.
- **sylph:** very fast, accurate at species level, reports both abundance types, and its ANI output helps tell strains apart. Only as good as its genome database (GTDB covers bacteria and archaea only), and uses GTDB names.

::::::::::::::::::::::::::::::::::::: challenge

### Challenge 2: Pick a profiler

Which profiler (or output) would you choose for each scenario, and why?

1. You want to know what fraction of the community belongs to species that are missing from reference databases.
2. You plan to use a count-based differential abundance method that needs read counts.
3. Your samples may contain DNA viruses and fungi as well as bacteria.
4. You want to compare cell proportions with flow-cytometry estimates.

:::::::::::::::::::::::: solution

1. **mOTUs.** Its `-1` fraction estimates the unassigned species. sylph with a MAG-rich database is also an option.
2. **Kraken2/Bracken** read counts, or mOTUs run with read counts enabled.
3. **Kraken2** with a database that includes viral and fungal genomes. MetaPhlAn also covers some eukaryotes and viruses.
4. **Taxonomic abundance:** sylph's `Taxonomic_abundance`, MetaPhlAn or mOTUs. Not Bracken's sequence abundance.

These are discussion points rather than single right answers.

:::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::

## Using outputs as complementary information

There is rarely one "best" profiler. A practical strategy:

1. **Choose a primary output** that matches your question: taxonomic abundance for questions about cell proportions, read counts for count-based statistical models.
2. **Use the other profilers as a sanity check.** Taxa detected by several tools, built on different databases and methods, are more trustworthy. Taxa detected by only one tool, at low abundance, deserve scrutiny.
3. **Report what you used:** tool, version, database and abundance type. Results from different profilers, or different abundance types, should not be mixed in one analysis.

::::::::::::::::::::::::::::::::::::: keypoints

- Profilers differ in their database (NCBI/RefSeq vs GTDB), method (marker genes, k-mers, sketching) and type of abundance, and all three affect the results.
- Sequence abundance (fraction of reads) and taxonomic abundance (fraction of cells) differ systematically because genome sizes differ.
- Containment ANI estimates how much of a reference genome in "contained" in the sample.
- Choose the profiler and abundance type that match your research question, and use agreement between tools as supporting evidence.

::::::::::::::::::::::::::::::::::::::::::::::::
