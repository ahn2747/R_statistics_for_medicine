# Adding a new cancer / an updated .sav

1. Put `<CANCER>.sav` in `database/`, add the code to `cfg$cancers`, and add `cfg$schema$<CANCER> = list(...)` with **only the fields that differ** from `schema_default`. Examples:
   - `stage = "ajcc_pathologic_tumor_stage", stage_format = "ajcc_text"`
   - `event_value = "Dead"`
   - `time_unit = "months"`
   - `neoadjuvant = NA` if the column doesn't exist
2. Add `cfg$value_labels$<CANCER>` and `cfg$label_rules$<CANCER>` only if the codes or cutoffs differ (for example a different CEA cutoff or age group).
3. Put the gene CSVs in `database/gene_files/` as `<CANCER>_<split>_<GENE>.csv` (format in docs/data_flow.md).
4. Run 01 and read everything it prints:
   - schema errors name the exact `cfg$schema` entry to fix
   - warnings flag the event rate (<5% or >70%), the time unit, possible reversed status coding (from missing `last_contact`), and sex values
   - the label-consistency results
5. **Updated .sav for an existing cancer:** 01 prints the manifest diff (n, survival counts, columns, genes, neoadjuvant IDs). Check that it matches what you meant to change. If the md5 is the same but the counts changed, a code change caused it: investigate before deleting the manifest entry.
6. Run 02b (check `merge_QC*.csv`), then 03 and 04. Check the survival summary for skipped genes, EPV, PH flags and failures.
7. For a new cancer, set `primary_gene`, or leave it as is if the same gene is primary. A primary gene missing from a cancer produces a warning.
