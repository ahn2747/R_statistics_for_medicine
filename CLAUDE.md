# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project
An R pipeline for a medical paper using TCGA data, currently colorectal cancer (COAD and READ).
- **COAD and READ are always analyzed separately. Never pool them.**
- The code is config-driven so more TCGA cancers can be added later with no code changes.
- The primary hypothesis gene is `cfg$primary_gene` (MS4A1). All other genes are exploratory.

## Commands
R 4.6.1 is not on PATH. Use `"C:/Program Files/R/R-4.6.1/bin/Rscript.exe" <script>` and run from the project root, because scripts use relative paths. Run in this order:

```
00_setup.R --core       # CRAN packages only (enough for 01–04); without --core it also installs Bioconductor
01_import_clinical.R    # schema check → database/<C>.sav → data/processed/<C>_clinical{,_all}.rds + db_manifest.json
02a_download_tcga.R     # GDC STAR-Counts via TCGAbiolinks → D:/GDCdata (slow; needed for GSEA / genes_from_gdc)
02b_merge_genes.R       # database/gene_files/*.csv → data/processed/<C>_merged.{rds,csv,sav} + merge_QC
03_table1.R             # table1_<GENE>.{docx,csv}, table1_overall, tss_by_group_<GENE>.csv
04_survival.R           # KM, Cox, strata(tss), cox.zph, time-split Cox, RMST, BH q, forest plots
05_gsea.R               # DESeq2 High vs Low (~ tss + group) + fgsea on MSigDB; needs 02a outputs and Bioconductor packages
```
- Not written yet: `06_external_geo.R` (GSE39582). `07_timer` and `08_proteomics` are optional.
- 05 is slow: DESeq2 runs twice, plus apeglm and GO:BP fgsea. By default it runs only `cfg$primary_gene`; `cfg$gsea_genes = "all"` runs every gene. Run it in the background.
- **There is no test suite.** Correctness comes from built-in stops:
  - `validate_schema()` and `check_label_consistency()` in 01
  - the manifest check in 01
  - the Status cross-check in 02b
- To test changes without touching real outputs:
  1. Source `config.R` and `R/utils.R`.
  2. Point `cfg$processed_dir`, `cfg$output_dir`, `cfg$manifest` and `cfg$gene_dir` at scratch folders.
  3. Evaluate the script body without its `source()` lines.
- For ad-hoc R containing non-ASCII text (≤, –, Korean), write a temporary .R file. Inline `Rscript -e` with Unicode has segfaulted.

## Architecture
- **`config.R`** holds `cfg`, the only place for project settings:
  - data: cancers, paths, `gene_pattern`, `drop_columns`, `exclude_neoadjuvant`, `primary_gene`
  - schema: `schema_default` plus per-cancer overrides in `cfg$schema$<C>`
  - labels: per-cancer overrides in `cfg$value_labels$<C>` and `cfg$label_rules$<C>`
  - Table 1: `table1_vars`
  - survival: `cox_covariates`, `stage_collapsed`, `epv_min`, `min_events`, `tss_min_n`, `km_times`, `ph_split_months`, `rmst_tau`
  - figures: `group_colors`
  - GSEA: `gsea_genes`, `gsea_design_covariates`, `gsea_extra_collections`, `gsea_highlight`, `gsea_size`, `n_cores`, `seed`
- **`R/utils.R`** holds all shared logic. Every script starts with `source("config.R"); source("R/utils.R")`, loops over `cfg$cancers`, and finds genes with `detect_genes()` (stems that have both `_expression` and `_group` columns). Don't hardcode cancer codes, genes or patient IDs.
- **Schema → canonical columns:** source column names live only in the schema. `apply_schema_names()` renames the schema columns to the canonical names in `canon`:

  | Schema field | Canonical column |
  |---|---|
  | id | `sample_id` |
  | time | `surv_time` |
  | status | `status` |
  | stage | `pathologic_stage` |
  | age | `age` |
  | sex | `gender` |
  | neoadjuvant | `neoadjuvant` |
  | last_contact | `last_contact` |

  - Downstream code uses only canonical names.
  - All other columns go through `standardize_names()`: gene columns become `<gene>_expression` / `<gene>_group`, and the rest go through `janitor::clean_names`. A non-schema column whose cleaned name collides with a canonical one gets the suffix `_src`.
- **Clinical data flow:** `validate_schema()` runs first, then `load_clinical()` and `recode_clinical(df, cancer)`.
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
  - If `<C>.sav` doesn't exist, `load_gdc_clinical()` builds the canonical columns from GDC. This path is untested.
- **`data/processed/db_manifest.json`** stores, for each cancer: the .sav md5, n, survival n/events (with and without neoadjuvant patients), columns, genes, and neoadjuvant IDs.
  - Same md5 but different n or survival counts: 01 stops, because a code change altered the results. If the change was intended, delete that cancer's entry and rerun.
  - Changed md5: 01 prints a diff and continues.
  - The manifest replaced the old hardcoded `expected_surv`.
- **Neoadjuvant exclusion:** 01 saves both `_clinical_all.rds` and `_clinical.rds`. 02b merges onto `_all`, because the groups were defined on the full cohort, and applies `cfg$exclude_neoadjuvant` at the end.
- **Gene CSVs** (`<CANCER>_<split>_<GENE>.csv` with columns `Patient, Days, Status, Expression, Group`):
  - Only Expression and Group are merged. Status is used only for QC.
  - A gene that already exists in the .sav is QC-only.
- **Output helpers:**
  - All of these take an optional `subdir` (e.g. `"gsea/MS4A1"`) under `output/<type>/<C>/`.
  - `save_table()` (CSV, UTF-8 BOM) and `save_df_table()` (.docx + .csv)
  - `save_fig()` (cairo PDF + 300-dpi LZW TIFF)
  - `move_stale_dirs()`: gene subfolders under `gsea/` that aren't in the current run are moved to `gsea/_stale/`.
  - `fmt_hr()` and `fmt_p()` (2 decimals; `<0.001`)
  - `move_stale_outputs()`: gene-specific files that this run didn't recreate are moved to `output/*/<C>/_stale/`, never deleted. It uses the prefixes each script owns:
    - 03: `table1_`, `tss_by_group_`
    - 04: `km_`, `forest_multi_`, `zph_`, `cox_multi_`
- **04 outputs** (`output/tables/<C>/`):
  - `survival_summary.{docx,csv}`: one row per gene, with the primary gene first and a `Role` column. It includes time-split HRs, RMST, the exploratory flag and a note.
  - `survival_summary_raw.csv`
  - `km_summary`, `cox_uni`, `cox_multi_<GENE>`, `ph_tests.csv`
  - Figures: `km_`, `forest_multi_`, `forest_genes`, and `zph_<GENE>.pdf` (only when PH is violated).
- **05 inputs:**
  - `data/processed/<C>_counts.rds`: raw STAR unstranded counts from 02a (symbol × 12-character patient ID, tumor samples only, first aliquot)
  - `<C>_merged.rds`
  - `database/TCGA_<C>_RNAseq_Expression.csv` (log2 TPM+1, used for QC only)
  - the MSigDB `.chip` file
- **05 outputs:**
  - `output/tables/<C>/gsea/<GENE>/`: `de_results.csv`, `gsea_<Hallmark|Reactome|KEGG|GOBP|C8>.csv`, `gsea_QC.csv`, `gsea_summary.{docx,csv}`
  - `output/figures/<C>/gsea/<GENE>/`: `volcano`, `pca`, `nes_<collection>`, `enrichment_<PATHWAY>`
  - `data/processed/<C>_dds_<GENE>.rds`
- **msigdbr 26.x** uses `collection` / `subcollection` (not `category`), and its data download on first use. KEGG uses `CP:KEGG_LEGACY`, falling back to `CP:KEGG_MEDICUS`.

## Data rules agreed with the user (don't change silently)
- **Survival** comes from the schema's time and status (currently `Days`/`Status`, 1 = death). `os_months = days / 30.44`. Reference counts, now stored in the manifest:

  | | COAD n (events) | READ n (events) |
  |---|---|---|
  | all patients | 439 (97) | 157 (25) |
  | neoadjuvant excluded (default) | 436 (96) | 156 (25) |
- **`database/` is read-only input.** The user edits the .sav files between sessions, so re-inspect them rather than assuming their structure.
- **Expression groups:** `*_group` values from the .sav are never recomputed. `median_split()` is only for new genes that come without a group. The median QC allows a difference of 1 patient.
- **Indicator columns** (`kras_/braf_gene_analysis_indicator`, `mismatch_rep_proteins_tested_by_ihc`) mean a test was **performed**. They are not mutation or MMR results. Never derive mutation/MSI/MMR variables from them. They are also excluded from Table 1.
- **Codes and consistency rules** (`label_rules_default`; a mismatch stops the run):
  - `cea_g`: 0 = ≤5 ng/mL, 1 = >5; checked against `cea_level_pretreatment > 5`
  - `age_g`: ≤65 / ≥66; checked against `age > 65`
  - `pathologic_stage_12_34`: 0 = I–II, 1 = III–IV; checked against `stage > 2`
- **Cox models:**
  - The multivariable model is the gene group + age (per 10 years) + gender + stage, as a complete-case analysis.
  - EPV = events / number of parameters. If EPV < `epv_min`, `stage_12_34` replaces stage I–IV. If EPV is still low, the result is `exploratory = TRUE`; all READ models are (EPV about 5.8).
  - Sensitivity analysis: `strata(tss)`, with sites under 10 patients grouped as "Other". The expression group is strongly associated with TSS for many genes (`tss_by_group_*.csv`), so report this analysis.
- **PH violations:**
  - Time-split Cox (0–24 / >24 months, adjusted) is run only when the gene **group term** has cox.zph p < 0.05. Currently that is READ TIMP1 and READ APC.
  - COAD TIMP1 violates PH only through its stage term.
  - RMST (High − Low, tau = 60) is reported for every gene.
- **Multiple testing:** BH q-values are computed across the exploratory genes only. The primary gene reports its raw p, and its q is shown as "–".
- **Reference results** (checked against independent code):
  - COAD MS4A1: log-rank p 0.018; multivariable HR 0.55 (0.35–0.85), p 0.007; strata(TSS) HR 0.50; RMST difference 4.88 months
  - READ MS4A1: log-rank p 0.008; RMST difference 9.11 months
- **GSEA (05):**
  - Use the existing `<gene>_group` from the merged data (the same patients and grouping as survival). Never re-split.
  - Main design: `~ tss + group`, with sites under `tss_min_n` patients grouped as "Other" (groups are unbalanced across sites). Contrast: High vs Low.
  - Pre-filter: keep genes with count ≥ 10 in at least as many samples as the smaller group.
  - Rank by the DESeq2 Wald stat, remap symbols with the `.chip` file, and average duplicates.
  - `fgseaMultilevel` with minSize 15, maxSize 500, eps 0 and a fixed seed.
  - GO:BP runs `collapsePathways` on the top 300 significant pathways.
  - `robust_no_covariate` = same NES sign and padj < 0.05 in the `~ group` run.
  - QC: Spearman correlation between GDC log2 TPM and the .sav expression (warns if r < 0.8). The gene itself must rank in the top 1% up in High, otherwise a warning is printed.
- **Table 1:** Wilcoxon for continuous variables; chi-square, switching to Fisher's exact when any expected cell is < 5. Missing values are shown as a row and excluded from the tests.

## Adding a new cancer / an updated .sav
1. Put `<CANCER>.sav` in `database/`, add the code to `cfg$cancers`, and add `cfg$schema$<CANCER> = list(...)` with **only the fields that differ** from `schema_default`. Examples:
   - `stage = "ajcc_pathologic_tumor_stage", stage_format = "ajcc_text"`
   - `event_value = "Dead"`
   - `time_unit = "months"`
   - `neoadjuvant = NA` if the column doesn't exist
2. Add `cfg$value_labels$<CANCER>` and `cfg$label_rules$<CANCER>` only if the codes or cutoffs differ (for example a different CEA cutoff or age group).
3. Put the gene CSVs in `database/gene_files/` as `<CANCER>_<split>_<GENE>.csv`.
4. Run 01 and read everything it prints:
   - schema errors name the exact `cfg$schema` entry to fix
   - warnings flag the event rate (<5% or >70%), the time unit, possible reversed status coding (from missing `last_contact`), and sex values
   - the label-consistency results
5. **Updated .sav for an existing cancer:** 01 prints the manifest diff (n, survival counts, columns, genes, neoadjuvant IDs). Check that it matches what you meant to change. If the md5 is the same but the counts changed, a code change caused it: investigate before deleting the manifest entry.
6. Run 02b (check `merge_QC*.csv`), then 03 and 04. Check the survival summary for skipped genes, EPV, PH flags and failures.
7. For a new cancer, set `primary_gene`, or leave it as is if the same gene is primary. A primary gene missing from a cancer produces a warning.

## Conventions
- Code comments and console messages are in Korean. Table and figure labels are in English.
- On data inconsistencies, `stop()` and name the patient IDs, plus the `cfg` entry to fix when there is one. Don't silently drop or coerce.
- Keep flextable captions as plain text; markdown `**` shows up literally in the .docx.
- survminer 0.5.2 with ggplot2 4.x prints "Ignoring unknown labels". It is harmless and wrapped in `suppressMessages()`. Combine the KM plot and risk table with patchwork (`p$plot / p$table`).
- Build PDFs with `cairo_pdf` so that — ≤ – render correctly.
- Write CIs that can be negative (RMST) as `(lo, hi)`, not `lo–hi`.
- After changing `R/utils.R`, rerun 01 → 02b → 03 → 04. The manifest will stop 01 if the counts shifted.
