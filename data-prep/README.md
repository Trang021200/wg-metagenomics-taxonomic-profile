# Data preparation for Day 2

`prepare_day2_data.py` builds the files in `episodes/data/` from the nf-core/taxprofiler
results of the five Sun et al. (2021) simulated gut samples. See
`instructors/instructor-notes.md` for the command, inputs and caveats.

| File | What it is |
|---|---|
| `prepare_day2_data.py` | Reads the per-sample MetaPhlAn, mOTUs, Bracken, Kraken2, sylph and fastp outputs and writes species-level tables, the truth tables, the NCBI↔GTDB mapping and the sample metadata. Python 3 standard library only. |
| `species_map_gtdb_r220.tsv` | Name-based mapping of the 95 truth species to GTDB R220 (type-strain genomes, built from GTDB's per-genome taxonomy history). The script refines 12 split species with sylph's genome-level evidence. |
