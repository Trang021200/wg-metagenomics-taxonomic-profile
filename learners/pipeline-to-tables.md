---
title: 'From pipeline output to analysis tables'
---

The tables in the workshop `data/` folder (`metaphlan.tsv`, `motus.tsv`, `bracken.tsv`, `sylph_taxonomic.tsv`, `sylph_sequence.tsv` and `metadata.tsv`) were made from the per-sample output files of nf-core/taxprofiler. This page shows how to do the same with your own pipeline results, so that you can run the Day 2 analysis on your own samples.

Each profiler writes **one file per sample, in its own format**. For the analysis we want **one table per profiler**, all with the same layout:

| taxonomy_id | name | sample_1 | sample_2 | ... |
|---|---|---|---|---|
| 29391 | Gemella morbillorum | 17.88 | 0.91 | ... |

with one row per species and one column per sample. Each value stays in the profiler's own unit (percent, or read counts). Converting to relative abundance is part of the analysis in [Community analysis](../episodes/03-community-analysis.Rmd).

::::::::::::::::::::::::::::::::::::: callout

### Where the files are

Inside the nf-core/taxprofiler output folder (`--outdir`), each profiler has its own folder, with one subfolder per database and one file per sample in it:

| Profiler | Files | What we keep |
|---|---|---|
| MetaPhlAn | `metaphlan/<database>/<sample>*.metaphlan_profile.txt` | species rows, NCBI taxonomy ID, relative abundance (%) |
| mOTUs | `motus/<database>/<sample>*.out` | species, NCBI taxonomy ID, read counts |
| Kraken2 + Bracken | `bracken/<database>/<sample>*.bracken.tsv` | species, NCBI taxonomy ID, re-estimated reads (`new_est_reads`) |
| sylph | `sylph/<database>/<sample>*.sylphmpa` | species, GTDB lineage, taxonomic and sequence abundance (%) |
| fastp | `fastp/<sample>*.fastp.json` | read pairs before and after quality control |
| Kraken2 | `kraken2/<database>/<sample>*.kraken2.report.txt` | number of unclassified reads |

::::::::::::::::::::::::::::::::::::::::::::::::

## Settings and helper functions

You need `readr`, `dplyr` and `tidyr` (installed in [Setup](../learners/setup.md)) and `jsonlite` for the fastp reports (`install.packages("jsonlite")` if you don't have it). Change `results` to your pipeline output folder and `samples` to your sample IDs, as in your samplesheet.

```r
library(readr)
library(dplyr)
library(tidyr)

results <- "taxprofiler_results"   # the pipeline's output folder (--outdir)
samples <- c("Gut_sample1", "Gut_sample2", "Gut_sample3", "Gut_sample4", "Gut_sample5")
out_dir <- "my_tables"             # a new folder, so the workshop data/ is not overwritten
dir.create(out_dir, showWarnings = FALSE)
```

Three small functions do the work that is the same for every profiler:

```r
# Find one output file per sample for a profiler, e.g. all MetaPhlAn profiles.
find_files <- function(folder, pattern) {
  files <- list.files(file.path(results, folder), pattern = pattern,
                      recursive = TRUE, full.names = TRUE)
  if (length(files) != length(samples)) {
    stop("Expected ", length(samples), " files in ", folder, ", found ", length(files))
  }
  files
}

# Which sample does a file belong to? File names start with the sample ID,
# followed by "_" or "." (e.g. Gut_sample1_db1.bracken.tsv).
sample_of <- function(file) {
  for (s in samples) {
    if (grepl(paste0("^", s, "[._]"), basename(file))) return(s)
  }
  stop("No sample ID found in ", basename(file))
}

# Long table (sample_id, taxonomy_id, name, value) -> one row per taxon,
# one column per sample, most abundant taxa first.
to_wide <- function(long) {
  long |>
    group_by(taxonomy_id) |>
    mutate(name = first(name)) |>
    group_by(taxonomy_id, name, sample_id) |>
    summarise(value = sum(value), .groups = "drop") |>
    pivot_wider(names_from = sample_id, values_from = value, values_fill = 0) |>
    select(taxonomy_id, name, all_of(samples)) |>
    arrange(desc(rowSums(across(all_of(samples)))))
}
```

For each profiler, we write a function that reads **one** sample's file into a long table with four columns (`sample_id`, `taxonomy_id`, `name`, `value`), apply it to all the files with `lapply()`, stack the results with `bind_rows()`, and reshape with `to_wide()`. A species missing from a sample gets a 0.

## MetaPhlAn

A MetaPhlAn profile starts with a few lines beginning with `#`, then has one row per clade at every rank, from kingdom (`k__`) down to species (`s__`) and SGB (`t__`). The `clade_name` and `NCBI_tax_id` columns are full paths separated by `|`, e.g. `k__Bacteria|...|s__Gemella_morbillorum`. We keep only the rows that end at a species, and take the last taxonomy ID in the path. Some MetaPhlAn species are SGBs with no NCBI species (e.g. *Rothia* SGB49305), so their ID is empty. We give those an ID of their own, `MetaPhlAn:<name>`.

```r
read_metaphlan <- function(file) {
  read_tsv(file, comment = "#", quote = "",
           col_names = c("clade_name", "ncbi_tax_id", "relative_abundance", "additional_species"),
           col_types = "ccdc") |>
    filter(grepl("(^|\\|)s__[^|]*$", clade_name)) |>     # species rows only
    mutate(species     = sub(".*s__", "", clade_name),
           taxonomy_id = sub(".*\\|", "", ncbi_tax_id),    # last ID in the path
           taxonomy_id = ifelse(taxonomy_id == "", paste0("MetaPhlAn:", species), taxonomy_id),
           name        = gsub("_", " ", species),
           sample_id   = sample_of(file)) |>
    select(sample_id, taxonomy_id, name, value = relative_abundance)
}

metaphlan <- find_files("metaphlan", "metaphlan_profile\\.txt$") |>
  lapply(read_metaphlan) |>
  bind_rows() |>
  to_wide()
```

## mOTUs

A mOTUs profile has three columns: the mOTU name with its ID in brackets (`Gemella morbillorum [ref_mOTU_v31_...]`), the NCBI taxonomy ID, and the read count. Several mOTUs can belong to the same NCBI species, so `to_wide()` adds up their counts. mOTUs without an NCBI species (`NA`) keep their mOTU ID, and the `unassigned` row (reads from species missing from the mOTUs database) gets the ID `-1`.

```r
read_motus <- function(file) {
  read_tsv(file, comment = "#", quote = "",
           col_names = c("label", "ncbi_tax_id", "count"), col_types = "ccd") |>
    filter(count > 0) |>
    mutate(motu_id     = sub(".*\\[(.*)\\]$", "\\1", label),
           name        = sub("\\s*\\[[^]]*\\]$", "", label),
           taxonomy_id = case_when(
             startsWith(label, "unassigned") ~ "-1",
             is.na(ncbi_tax_id) ~ paste0("mOTUs:", motu_id),   # no NCBI species
             TRUE ~ ncbi_tax_id),
           sample_id   = sample_of(file)) |>
    select(sample_id, taxonomy_id, name, value = count)
}

motus <- find_files("motus", "\\.out$") |>
  lapply(read_motus) |>
  bind_rows() |>
  to_wide()
```

## Kraken2 + Bracken

Bracken's species report is already one row per species, with a header. `new_est_reads` is Bracken's estimate of the reads per species, after redistributing the reads that Kraken2 could only assign to a higher rank.

```r
read_bracken <- function(file) {
  read_tsv(file, quote = "", col_types = cols(taxonomy_id = "c", .default = col_guess())) |>
    transmute(sample_id = sample_of(file), taxonomy_id, name, value = new_est_reads)
}

bracken <- find_files("bracken", "\\.bracken\\.tsv$") |>
  lapply(read_bracken) |>
  bind_rows() |>
  to_wide()
```

## sylph

sylph-tax writes a MetaPhlAn-style table with GTDB names: the lineage (`d__Bacteria|...|s__Gemella morbillorum`), the taxonomic abundance, the sequence abundance, and the ANI for genome-level (`t__`) rows. One file gives two tables, one per abundance type. GTDB has no taxonomy IDs, so the species name is used as the ID, and the full lineage is kept in a `gtdb_lineage` column.

```r
read_sylph <- function(file) {
  read_tsv(file, comment = "#", quote = "", col_types = cols(.default = "c")) |>
    select(clade_name = 1, taxonomic = 2, sequence = 3) |>
    filter(grepl("(^|\\|)s__[^|]*$", clade_name)) |>
    mutate(across(c(taxonomic, sequence), as.numeric),
           name         = sub(".*s__", "", clade_name),
           taxonomy_id  = name,                          # GTDB has no taxonomy IDs
           gtdb_lineage = gsub("|", ";", clade_name, fixed = TRUE),
           sample_id    = sample_of(file))
}

sylph_long <- find_files("sylph", "\\.sylphmpa$") |>
  lapply(read_sylph) |>
  bind_rows()
lineages <- distinct(sylph_long, taxonomy_id, gtdb_lineage)

sylph_taxonomic <- sylph_long |>
  mutate(value = taxonomic) |>
  to_wide() |>
  left_join(lineages, by = "taxonomy_id") |>
  relocate(gtdb_lineage, .after = name)

sylph_sequence <- sylph_long |>
  mutate(value = sequence) |>
  to_wide() |>
  left_join(lineages, by = "taxonomy_id") |>
  relocate(gtdb_lineage, .after = name)
```

## Sample metadata

The workshop `metadata.tsv` adds two quality numbers from the pipeline to each sample: the number of read pairs (from the fastp report) and the percentage of reads Kraken2 left unclassified (from the Kraken2 report). For your own study, also add the information about each sample that your question is about, such as the group, time point or patient.

```r
read_fastp <- function(file) {
  qc <- jsonlite::fromJSON(file)$summary
  tibble(sample_id      = sample_of(file),
         read_pairs_raw = qc$before_filtering$total_reads / 2,
         read_pairs_qc  = qc$after_filtering$total_reads / 2)
}

read_kraken2_report <- function(file) {
  report <- read_tsv(file, col_names = FALSE, quote = "", col_types = cols(.default = "c"))
  unclassified <- sum(as.numeric(report$X2[report$X5 == "0"]))   # taxid 0
  classified   <- sum(as.numeric(report$X2[report$X5 == "1"]))   # taxid 1 = root
  tibble(sample_id = sample_of(file),
         kraken2_unclassified_pct = round(100 * unclassified / (unclassified + classified), 2))
}

metadata <- bind_rows(lapply(find_files("fastp", "\\.fastp\\.json$"), read_fastp)) |>
  left_join(bind_rows(lapply(find_files("kraken2", "kraken2\\.report\\.txt$"),
                             read_kraken2_report)),
            by = "sample_id")
```

## Save the tables

```r
write_tsv(metaphlan,       file.path(out_dir, "metaphlan.tsv"))
write_tsv(motus,           file.path(out_dir, "motus.tsv"))
write_tsv(bracken,         file.path(out_dir, "bracken.tsv"))
write_tsv(sylph_taxonomic, file.path(out_dir, "sylph_taxonomic.tsv"))
write_tsv(sylph_sequence,  file.path(out_dir, "sylph_sequence.tsv"))
write_tsv(metadata,        file.path(out_dir, "metadata.tsv"))
```

These tables have the same layout as the workshop files, so the code in [Community analysis](../episodes/03-community-analysis.Rmd) works on them once you change `data/` to `my_tables/` in the file paths. Your own samples have no ground truth, so skip the parts that compare against `truth_*.tsv`.

::::::::::::::::::::::::::::::::::::: callout

### More than one database per profiler

If you ran a profiler with two databases, its folder has a subfolder for each, and `find_files()` stops because it finds two files per sample. Point it at one database, e.g. `find_files("metaphlan/<database>", ...)`.

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: callout

### What about Taxpasta?

nf-core/taxprofiler can also standardise the profiles itself with Taxpasta (the `taxpasta/` folder). Without a taxonomy, Taxpasta's tables contain taxonomy IDs only, with all ranks mixed in one table. To get names and ranks, run the pipeline with `--taxpasta_taxonomy_dir` (an NCBI taxdump folder), `--taxpasta_add_name`, `--taxpasta_add_rank` and `--taxpasta_add_lineage`. Taxpasta 0.7.0 does not support sylph, so you still need the code above for sylph.

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: instructor

The workshop files were made with `data-prep/prepare_day2_data.py`, which applies the same rules in Python and also builds the ground-truth tables and the NCBI–GTDB mapping. The R code on this page gives the same five profiler tables and the same read-pair and Kraken2 columns of `metadata.tsv`.

::::::::::::::::::::::::::::::::::::::::::::::::
