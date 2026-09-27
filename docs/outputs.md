# Outputs

**File naming:** every table and figure is `<prefix>_<analysis>[_<GENE>].<ext>`. The prefix is the cancer code for 01–05 (e.g. `COAD_table1_MS4A1.docx`, `READ_km_TIMP1.pdf`, `COAD_volcano.pdf` under `gsea/<GENE>/`) and the GSE ID for 06 (`GSE39582_km_OS.pdf`). The file lists below omit the prefix.

## Output helpers (`R/utils.R`)
- All of these take an optional `subdir` (e.g. `"gsea/MS4A1"`) under `output/<type>/<C>/`.
- `save_table()`, `save_df_table()` and `save_fig()` add the prefix themselves (`out_file(prefix, name)`; `prefix` defaults to `cancer`). Code that builds a path by hand (03's Table 1 .docx, 04's `zph_` PDF) must call `out_file()` too.
- 06 writes with `cancer = "GEO"`, `subdir = <GSE>` and `prefix = <GSE>`, i.e. `output/<type>/GEO/<GSE>/<GSE>_*`.
- `save_table()` (CSV, UTF-8 BOM) and `save_df_table()` (.docx + .csv)
- `save_fig()` (cairo PDF + 300-dpi LZW TIFF)
- `fmt_hr()` and `fmt_p()` (2 decimals; `<0.001`)
- `move_stale_outputs()`: gene-specific files that this run didn't recreate are moved to `output/*/<C>/_stale/`, never deleted. Files without the `<prefix>_` start (the old naming) are moved too. It matches `<prefix>_` + the families each script owns:
  - 03: `table1_`, `tss_by_group_`
  - 04: `km_`, `forest_multi_`, `zph_`, `cox_multi_`
  - 06 (under `GEO/<GSE>/`): `km_`, `forest`, `cox_`, `zph_`
- `move_stale_dirs()`: gene subfolders under `gsea/` that aren't in the current run are moved to `gsea/_stale/`.

## 02a outputs (`data/processed/`, not under `output/`, no prefix rule)
- `<C>_counts.rds`, `<C>_se.rda` and `database/TCGA_<C>_RNAseq_Expression.csv` (the existing data files; see docs/data_flow.md)
- `<C>_gdc_provenance.json`: GDC provenance for the Methods.
  - `gdc_data_release`: the release at download (`release` and `release_source`), `metadata_se` and the release current when the JSON was made. A `release_source` of `queried_at_<date>, not at download` means the download-time release was unavailable.
  - `query`: project / data.category / data.type / workflow.type, and the number of files
  - `assays`: counts and TPM assays used
  - `reference`: gene model (from the cached TSV headers) and genome build with its source
  - `samples`: tumor filter, tumor n, duplicate patients (kept and dropped aliquots) and final patient n
  - `genes`: gene counts and the final matrix size
  - `cache`: GDC cache folder, file mtime min/max as download-time evidence, and md5 agreement with GDC
  - `outputs_md5`, and `outputs_counts_identical` (`--provenance-only` rebuilds the counts from `_se.rda` and compares them)
  - `software`: R, TCGAbiolinks and SummarizedExperiment versions
- `<C>_gdc_files.csv`: one row per GDC file, with file_id, file_name, barcode, md5sum, size, version, created/updated_datetime, data_release, workflow version, cache mtime/md5, header line, and the `in_query`, `in_cache`, `md5_match`, `in_se` and `used` flags. `used` marks the aliquots in the final matrix.
- Run `Rscript 02a_download_tcga.R --provenance-only [C]` to regenerate these two files from `<C>_se.rda` and the GDC cache without downloading. It stops if `<C>_se.rda` is missing and writes no data files.

## 03 outputs
- `table1_<GENE>.{docx,csv}`, `table1_overall`, `tss_by_group_<GENE>.csv`

## 04 outputs (`output/tables/<C>/`)
- `survival_summary.{docx,csv}`: one row per gene, with the primary gene first and a `Role` column. It includes time-split HRs, RMST, the exploratory flag and a note.
- `survival_summary_raw.csv`
- `km_summary`, `cox_uni`, `cox_multi_<GENE>`, `ph_tests.csv`
- Figures: `km_`, `forest_multi_`, `forest_genes`, and `zph_<GENE>.pdf` (only when PH is violated).

## 05 outputs
- `output/tables/<C>/gsea/<GENE>/`: `de_results.csv`, `gsea_<Hallmark|Reactome|KEGG|GOBP|C8>.csv`, `gsea_QC.csv`, `gsea_summary.{docx,csv}`
- `output/figures/<C>/gsea/<GENE>/`: `volcano`, `pca`, `nes_<collection>`, `enrichment_<PATHWAY>` (selection rules in docs/methods.md)
- `data/processed/<C>_dds_<GENE>.rds`

## 06 outputs
- `output/tables/GEO/<GSE>/`:
  - `probe_QC.csv`: every probe matching `cfg$primary_gene` (mean, median, IQR, SD over tumors, multi-mapping flag, pairwise r, `chosen`, rule)
  - `pheno_QC.csv`: every mapped field with its raw key, n, missing and values, plus usable n/events per endpoint
  - `survival_summary.{docx,csv}` and `survival_summary_raw.csv`: one row per endpoint × population × model. Includes EPV, PH p for the gene term, the interaction p (+ extra-covariate model), log-rank p, RMST, and the direction compared with the TCGA multivariable HR.
  - `km_summary.{docx,csv}`: median and 3-/5-year survival by group, plus log-rank p
  - `cox_<EP>.{docx,csv}` and `cox_<EP>_<subgroup>`: every covariate row of every model
  - `ph_tests.csv`: cox.zph for every model
- `output/figures/GEO/<GSE>/`: `km_<EP>`, `km_<EP>_<subgroup>`, `forest` (gene HR across all models) and `forest_multi_<EP>` (all covariates of the most-adjusted model)
- `inspect/` (from `test_geo_inspect.R`) is left alone.
- Stale prefixes: `km_`, `forest`, `cox_`, `zph_`.
