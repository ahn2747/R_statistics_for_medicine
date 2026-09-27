# Data flow

## `config.R` key groups
`cfg` contains:
- data: cancers, paths, `gene_pattern`, `drop_columns`, `exclude_neoadjuvant`, `primary_gene`
- schema: `schema_default` plus per-cancer overrides in `cfg$schema$<C>`
- labels: per-cancer overrides in `cfg$value_labels$<C>` and `cfg$label_rules$<C>`
- Table 1: `table1_vars`
- covariates, strata, reference levels and exclusions: see docs/covariates.md
- survival: `km_times`, `ph_split_months`, `rmst_tau`
- figures: `group_colors`
- GSEA: `gsea_genes`, `gsea_design_covariates`, `gsea_sensitivity_covariates`, `gsea_extra_collections`, `gsea_highlight`, `gsea_c8_pattern`, `gsea_enrichment_top`, `gsea_size`, `n_cores`, `seed`

## Column names
- Schema fields are renamed to canonical names by `apply_schema_names()` (mapping in CLAUDE.md).
- All other columns go through `standardize_names()`: gene columns become `<gene>_expression` / `<gene>_group`, and the rest go through `janitor::clean_names`. A non-schema column whose cleaned name collides with a canonical one gets the suffix `_src`.

## Clinical data flow
`validate_schema()` runs first, then `load_clinical()` and `recode_clinical(df, cancer)`.
- `load_clinical()` reads the file, renames columns, runs `clean_tcga()`, then:
  - converts status to 0/1 using `event_value`
  - converts stage to 1–4 (`stage_format` "numeric" or "ajcc_text")
  - computes `os_months` from `time_unit`
- `recode_clinical()`:
  - creates `pathologic_stage_12_34` if it's missing
  - runs `check_label_consistency()`
  - applies `labels_for(cancer)`
  - adds the `stage` alias
  - turns YES/NO columns into factors
  - adds `neoadjuvant_flag` and `tss` (characters 6–7 of `sample_id`)
- If `<C>.sav` doesn't exist, `load_gdc_clinical()` builds the canonical columns from GDC (untested; see docs/notes.md).

## Gene CSVs
`database/gene_files/<CANCER>_<split>_<GENE>.csv` with columns `Patient, Days, Status, Expression, Group`:
- Only Expression and Group are merged. Status is used only for QC.
- A gene that already exists in the .sav is QC-only.

## 05 inputs
- `data/processed/<C>_counts.rds`: raw STAR unstranded counts from 02a (symbol × 12-character patient ID, tumor samples only, first aliquot)
- `<C>_merged.rds`
- `database/TCGA_<C>_RNAseq_Expression.csv` (log2 TPM+1, used for QC only)
- the MSigDB `.chip` file
