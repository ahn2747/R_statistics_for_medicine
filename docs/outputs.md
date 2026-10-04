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
  - 03: `table1_`, `tss_by_group_`, `correlation_`
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
- `correlation_<GENE>.csv`: Pearson r/p/n (long: var1, var2, r, p, n; every ordered pair plus the diagonal) between the gene, `cfg$cor_vars$genes` (log2 expression) and `cfg$cor_vars$clinical` (raw values), in the gene's Table 1 patients. Feeds Table 2·3 of `<GENE>_tables.docx`.

## 04 outputs (`output/tables/<C>/`)
- `survival_summary.{docx,csv}`: one row per gene, with the primary gene first and a `Role` column. It includes time-split HRs, RMST, the exploratory flag and a note.
- `survival_summary_raw.csv`
- `km_summary`, `cox_uni`, `cox_multi_<GENE>`, `ph_tests.csv`
- Figures: `km_`, `forest_multi_`, `forest_genes`, and `zph_<GENE>.pdf` (only when PH is violated).

## 05 outputs
- `output/tables/<C>/gsea/<GENE>/`: `de_results.csv`, `gsea_Hallmark.csv` (one `gsea_<collection>.csv` per `cfg$gsea_collections` entry), `gsea_QC.csv`, `gsea_summary.{docx,csv}`
- `output/figures/<C>/gsea/<GENE>/`: `volcano`, `pca`, `nes_<collection>`, `enrichment_<PATHWAY>` (selection rules in docs/methods.md)
- `data/processed/<C>_dds_<GENE>.rds`

## 06 outputs
- `output/tables/GEO/<GSE>/`:
  - `probe_QC.csv`: every probe matching `cfg$primary_gene` (mean, median, IQR, SD over tumors, multi-mapping flag, pairwise r, `chosen`, rule)
  - `pheno_QC.csv`: every mapped field with its raw key, n, missing and values, plus usable n/events per endpoint
  - `survival_summary.{docx,csv}` and `survival_summary_raw.csv`: one row per endpoint × population × model. Includes EPV, PH p for the gene term, the interaction p (+ extra-covariate model), log-rank p, RMST, and the direction compared with the TCGA multivariable HR. The raw file also records `gene` and `cox_gene_term`, which 90_export.R checks before copying GEO results.
  - `km_summary.{docx,csv}`: median and 3-/5-year survival by group, plus log-rank p
  - `cox_<EP>.{docx,csv}` and `cox_<EP>_<subgroup>`: every covariate row of every model, including one univariable model per clinical covariate (Model = "Univariable")
  - `ph_tests.csv`: cox.zph for every model
- `output/figures/GEO/<GSE>/`: `km_<EP>`, `km_<EP>_<subgroup>`, `forest` (gene HR across all models) and `forest_multi_<EP>` (all covariates of the most-adjusted model)
- `inspect/` (from `test_geo_inspect.R`) is left alone.
- Stale prefixes: `km_`, `forest`, `cox_`, `zph_`.

## Export to the student folder (`90_export.R`, settings in `cfg$export`)
`Rscript 90_export.R [GENE] [--dry-run]` copies one gene's finished outputs (GENE defaults to `cfg$primary_gene`) into `<dest_root>/<folder_fmt % GENE>/`, e.g. `0.학생연구자료/MAD2L1(COAD&READ)/`. Shared helpers are in `R/export.R`.
- **Layout** (`cfg$export$layout`; patterns use `{G}` = gene, `{coll}` = each `cfg$gsea_collections` name, `*` = wildcard, without the `<prefix>_` and extension):
  - `firstline/`: 03/04 CSVs (`table1_{G}`, `cox_uni`, `cox_multi_{G}`, `survival_summary`, `km_summary`, `ph_tests`, `correlation_{G}`)
  - `figure/`: `km_{G}`, `forest_multi_{G}` TIFFs
  - `gsea/data/`, `gsea/figure/`: from `gsea/<GENE>/` (`de_results` only if `include_de = TRUE`)
  - `geo/<GSE>_<validates>/data|figure/`: 06 outputs, only when the GEO results belong to this gene
  - Only `fig_ext` (TIFF) and `table_ext` (CSV) files are copied. `table1_overall`, `forest_genes`, PDFs and the per-table .docx files are not. `~$*` lock files, `_stale/` and `inspect/` are excluded (`cfg$export$exclude`).
- **Checks** (stop with the script to rerun): the gene is in every `<C>_merged.rds`; every TCGA pattern has a file (`table1_`/`correlation_` → 03, the rest → 04). Outputs older than `<C>_merged.rds` give a warning. A missing `gsea/<GENE>/` skips GSEA with a warning (and notes a `gsea/_stale/<GENE>/` copy, which is never used). GEO is copied only if `survival_summary_raw.csv`'s `gene` column (or, for old files, the chosen probe's symbol in `probe_QC.csv`) equals GENE; otherwise it is skipped with a "rerun 06 with MEDIN_GENE=<GENE>" note.
- **Copy rule:** names are kept. If the destination exists with the same content md5 (for .docx, the md5 of the zip entries without `docProps/core.xml`) the file is `same` and skipped; otherwise the old file is moved to `_old/<YYYYMMDD_HHMM>/<relative path>` (never deleted) and replaced. Files in the destination that aren't in the plan (manuscripts, manual files) are not touched. Copy failures are listed at the end and the exit code is 1. `--dry-run` prints the plan (new / same / replace / skip + reason) and changes nothing.
- **Generated files:**
  - `<GENE>_key_results.csv`: one row per cohort × endpoint × population × model (TCGA: univariable, multivariable, strata, time-split if present, from 04's `survival_summary_raw.csv`; GEO: every row of 06's raw file) with gene_term, n, events, HR, lo, hi, p, logrank_p, rmst_diff, ph_flag, exploratory; plus one `GSEA <collection>` row per cancer with the number of pathways with padj < `gsea_fdr` and the top `gsea_top` up/down by NES.
  - `<GENE>_tables.docx`: `cfg$export$template` (read only) filled from the pipeline CSVs. Captions are the template's, with its gene (auto-detected from the Table 1 caption, or `template_gene`) replaced and `caption_fixes` applied. Font and size are read from the template's tables; tables use three rules (header top/bottom, body bottom). Table 1 = `table1_<GENE>.csv` rows (`cfg$table1_vars`, with Missing rows) for every cancer side by side, High | Low | P-value (`table1_alt_first`), counts only (`table1_counts_only`). Table 2·3 = `correlation_<GENE>.csv` as an R/P matrix. Figure 1·2 = `km_<GENE>.tiff`. Table 4·5 = one row per multivariable term, "<Variable> (<level> vs <reference>)", with the univariable HR from `cox_uni` on the same scale (gene term per `cox_gene_term`) and a note with the model's n, events, EPV and the gene HR scale. Template Table 1 variables with no pipeline match are printed (`template_aliases` maps renamed ones).
  - `_export_manifest.csv`: source (or `generated`), dest, status, md5, size, source_mtime, exported_at, git_head, git_dirty, gene, cox_gene_term. Rewritten only when something changed (the previous one goes to `_old/`).
- `run_gene.R <GENE> [--skip-download] [--from=<step>]` runs `csv_download.py` → 02b → 03 → 04 → 05 (each cancer) → 06 → 90 with `MEDIN_GENE=<GENE>`, stopping at the first non-zero exit. Steps: `download`, `02b`, `03`, `04`, `05`, `06`, `90`. Timings go to `output/run_logs/<GENE>_<timestamp>.log` (git-ignored). It never commits; it lists the changed output files at the end.
