# CLAUDE.md

## Project
An R pipeline for a medical paper using TCGA data, currently colorectal cancer (COAD and READ).
- **COAD and READ are always analyzed separately. Never pool them.**
- The code is config-driven so more TCGA cancers can be added later with no code changes.
- The primary hypothesis gene is `cfg$primary_gene`; all other genes are exploratory. Several genes are analyzed, so never assume a specific gene or its results.

## Commands
R 4.6.1 is not on PATH. Use `"C:/Program Files/R/R-4.6.1/bin/Rscript.exe" <script>` and run from the project root, because scripts use relative paths. Run in this order:

```
00_setup.R --core       # CRAN packages only (enough for 01–04); without --core it also installs Bioconductor
01_import_clinical.R    # schema check → database/<C>.sav → data/processed/<C>_clinical{,_all}.rds + db_manifest.json
02a_download_tcga.R     # GDC STAR-Counts via TCGAbiolinks → D:/GDCdata (slow; needed for GSEA / genes_from_gdc) + <C>_gdc_provenance.json, <C>_gdc_files.csv
                        #   --provenance-only [C]: rebuild only the provenance files from <C>_se.rda + GDC cache (no download)
02b_merge_genes.R       # database/gene_files/*.csv → data/processed/<C>_merged.{rds,csv,sav} + merge_QC
03_table1.R             # table1_<GENE>.{docx,csv}, table1_overall, tss_by_group_<GENE>.csv
04_survival.R           # KM, Cox, strata(tss), cox.zph, time-split Cox, RMST, BH q, forest plots
05_gsea.R [COAD]        # DESeq2 High vs Low (~ tss + group) + fgsea on MSigDB Hallmark only (cfg$gsea_collections); needs 02a outputs; run one cancer per process
06_external_geo.R [GSE] # GEO validation of cfg$primary_gene (cfg$geo_datasets; needs GEOquery, i.e. 00_setup.R without --core; reads 04's <C>_survival_summary_raw.csv)
11_sig_scores.R [SIG ...] [C ...]   # signature track (separate from genes): z/PC1/ssGSEA/singscore per cancer → data/processed/<C>_sigscores_<SIG>.csv + QC under signature/<SIG>/ (cfg$signatures)
12_sig_survival.R [SIG ...] [C ...] # OS per score method (KM/RMST on median split, Cox per 1 SD, EPV, strata(tss), PH); primary = cfg$signature_primary_method; reads 11's score files + _clinical(_all).rds
13_sig_geo.R [GSE ...] [SIG ...]    # GEO: all-gene probe collapse (geo_probe_rule) → scores within cohort → OS/RFS × all/pMMR, 06's MMR models; HR direction vs 12 (needs 12 for the comparison)
20_paper_prep.R [C ...] [--refresh-gdc]   # paper track (lung): score check, GDC smoking (cached), ESTIMATE purity → data/processed/<C>_paper.rds (cfg$paper)
21_paper_analysis.R [C ...]          # items 1–8, 10 per score (z_sav main, ssgsea) → output/*/<C>/paper/
22_paper_geo.R [GSE ...]             # item 9: GEO validation (cfg$paper$geo, separate from 06/13) → output/*/GEO/paper/
90_export.R [GENE] [--dry-run]   # copy one gene's 03–06 outputs to cfg$export$dest_root + <G>_tables.docx, <G>_key_results.csv, _export_manifest.csv
run_gene.R <GENE> [--skip-download] [--from=<step>]   # csv_download.py → 02b → 03 → 04 → 05 ×cancers → 06 → 90, with MEDIN_GENE=<GENE>
```
- **Cancer sets:** `MEDIN_CANCER_SET` (default `colorectal` = COAD, READ; `lung` = LUAD, LUSC) selects `cfg$cancers`, the Cox covariates and the GEO datasets (`.cancer_sets` in config.R). Run each set separately, e.g. `MEDIN_CANCER_SET=lung Rscript 01_import_clinical.R`; never mix sets in one run (BH set, export folder and template are per set).
- **Python:** not on PATH; `cfg$python` falls back to the `py` launcher. `py csv_download.py GENE [GENE ...] [--cancers COAD READ]` downloads OncoLnc CSVs into `database/gene_files/` (an existing file with different content stops it; it never overwrites).
- **Primary gene override:** `cfg$primary_gene` reads `MEDIN_GENE` (default MAD2L1). `run_gene.R` sets it for child processes. Changing it changes 04's BH set (other genes' q) and which gene 05/06 run, so export each gene before moving on; `run_gene.R` never commits.
- **90:** settings in `cfg$export` (docs/outputs.md "Export"). `use_custom_filename = TRUE` + `custom_filename` renames the destination folder, generated files and the `{G}` position of copied file names (only there; results and content stay the real gene). Same-content files are skipped, changed ones moved to `_old/<timestamp>/` first; GEO is copied only if 06's `gene` column matches. The template `table양식.docx` is read only.
- `07_timer` and `08_proteomics` are optional and not written yet.
- **Signature track (11–13, `R/signature.R`):** never add signatures as fake genes to `merged.rds`; it must not affect `detect_genes()`, 04's BH set, 05, 90 or run_gene.R. Score files are the interface; settings in `cfg$signatures` / `cfg$signature_*` (docs/covariates.md, docs/methods.md). Needs GSVA, singscore (Bioconductor) and psych.
- **Paper track (20–22, `R/paper.R`, `cfg$paper`):** reads merged.rds and 11's scores but writes only `<C>_paper.rds`, `<C>_gdc_clinical.csv`, `<C>_purity_estimate.csv` and `paper/` folders. Needs tidyestimate (CRAN) and network on the first 20 run (GDC clinical). Run 12 again after a .sav change before comparing 21's ssGSEA results with it.
- **06:** the first run downloads the series matrix + GPL annotation into `cfg$geo_dir` (`data/geo`, git-ignored via `data/`) and caches `<GSE>_<GPL>_eset.rds`. Later runs take ~1 min. A GSE argument overrides `cfg$geo_datasets`.
- **05:** run one cancer per process: `Rscript 05_gsea.R COAD`, then `Rscript 05_gsea.R READ`. A command-line argument overrides `cfg$cancers`. Keep `n_cores` at 2: each Windows SnowParam worker holds its own copy of the data. By default only `cfg$primary_gene` runs; `cfg$gsea_genes = "all"` runs every gene.
- **There is no test suite.** Correctness comes from built-in stops: `validate_schema()` and `check_label_consistency()` in 01, the manifest check in 01, and the Status cross-check in 02b.
- To test changes without touching real outputs:
  1. Source `config.R` and `R/utils.R`.
  2. Point `cfg$processed_dir`, `cfg$output_dir`, `cfg$manifest` and `cfg$gene_dir` at scratch folders.
  3. Evaluate the script body without its `source()` lines.
- For ad-hoc R containing non-ASCII text (≤, –, Korean), write a temporary .R file. Inline `Rscript -e` with Unicode has segfaulted.

## Architecture
- **`config.R`** holds `cfg`, the only place for project settings.
- **`R/utils.R`** holds all shared logic. Every script starts with `source("config.R"); source("R/utils.R")`, loops over `cfg$cancers`, and finds genes with `detect_genes()` (stems that have both `_expression` and `_group` columns). Don't hardcode cancer codes, genes or patient IDs.
- **`R/survival.R`** holds the Cox/KM/forest/EPV/time-split/RMST helpers shared by 04 and 06 (both also source it). Time/event columns are arguments that default to 04's `os_months`/`status`. After changing it, rerun 04 and check that `git status` shows no changed CSVs under `output/tables/<C>/`.
- **`R/geo.R`** holds the GEO helpers shared by 06 and 13 (eset download/cache, log2 check, `geo_symbols()`/`probe_rule_col()`/`probe_qc()`, characteristics parsing, field mapping, `recode_geo()`, `pheno_qc()`, `geo_load_cohort()`, `geo_endpoint_data()`). It was split out of 06 as a pure refactor; after changing it, rerun 06 with the gene of the committed GEO outputs (`MEDIN_GENE=<gene in survival_summary_raw.csv>`) and check that no CSV under `output/tables/GEO/` changed.
- **`R/export.R`** holds the 90_export.R helpers (file matching, md5 copy plan, GEO gene check, key results, `build_tables_docx()`). `cfg$export` and `cfg$cor_vars` (03 correlation) are in `config.R`.
- **GEO datasets** (06) are config-only: `cfg$geo_datasets$<GSE>` maps canonical names to raw `characteristics_ch1` keys, and sets the platform, sample types to exclude, endpoints, factor levels, extra covariate and subgroup. Don't hardcode GSE IDs, probe IDs or pheno keys.
- **Covariate names are canonical column names** (`age_g`, not `ageG`). 04/06 stop on a configured covariate that isn't in the data (`covariate_terms()`); a GEO dataset lists the ones it lacks in `cox_covariates_drop`.
- **No hardcoded covariates or levels.** Every covariate, stratum, reference level, adjustment and exclusion is read from `config.R`. Nothing is hardcoded in 01–06, `R/utils.R` or `R/survival.R`: an audit grep for `"Low"`, `"High"`, `"Male"`, `/ 10`, `"tss"`, `"Other"`, `"stage"`, `"yes"` and `~ group` finds only structural schema field names and the KM/log-rank formula `~ group`. Keep it that way; add new settings to `config.R` and to the table in docs/covariates.md.
- **Canonical columns:** source column names live only in the schema. `apply_schema_names()` renames the schema fields to the canonical names in `canon` (id → `sample_id`, time → `surv_time`, status → `status`, stage → `pathologic_stage`, age → `age`, sex → `gender`, neoadjuvant → `neoadjuvant`, last_contact → `last_contact`). Downstream code uses only canonical names.
- **`data/processed/db_manifest.json`** stores, for each cancer: the .sav md5, n, survival n/events (with and without neoadjuvant patients), columns, genes, and neoadjuvant IDs.
  - Same md5 but different n or survival counts: 01 stops, because a code change altered the results. If the change was intended, delete that cancer's entry and rerun.
  - Changed md5: 01 prints a diff and continues.
- **Neoadjuvant exclusion:** 01 saves both `_clinical_all.rds` and `_clinical.rds`. 02b merges onto `_all`, because the groups were defined on the full cohort, and applies `cfg$exclude_neoadjuvant` at the end.

## Data rules agreed with the user (don't change silently)
- **`database/` is read-only input.** The user edits the .sav files between sessions, so re-inspect them rather than assuming their structure.
- **Expression groups:** `*_group` values from the .sav are never recomputed. `median_split()` is only for new genes that come without a group. The median QC allows a difference of 1 patient.
- **Indicator columns** (`kras_/braf_gene_analysis_indicator`, `mismatch_rep_proteins_tested_by_ihc`) mean a test was **performed**. They are not mutation or MMR results. Never derive mutation/MSI/MMR variables from them. They are also excluded from Table 1.
- **MSI/MMR status: pending.** Neither .sav has an MSI or MMR result column, so MSI is not a covariate in any model yet. If the user supplies one, add it through the schema or a data column, then add it to `cox_covariates` (check EPV; READ is already below 10) and `table1_vars`, consider it for `gsea_design_covariates`, and rerun 01 → 05.
- **Codes and consistency rules** (`cfg$label_rules_default` in config.R; a mismatch stops the run):
  - `cea_g`: 0 = ≤5 ng/mL, 1 = >5; checked against `cea_level_pretreatment > 5`
  - `age_g`: ≤65 / ≥66; checked against `age > 65`
  - `pathologic_stage_12_34`: 0 = I–II, 1 = III–IV; checked against `stage > 2`
- **Survival** comes from the schema's time and status (currently `Days`/`Status`, 1 = death). `os_months = days / 30.44`. Reference counts, stored in the manifest:

  | | COAD n (events) | READ n (events) |
  |---|---|---|
  | all patients | 439 (97) | 157 (25) |
  | neoadjuvant excluded (default) | 436 (96) | 156 (25) |

## Regression checks (gene-independent)
Gene-level results (HR, p, RMST, GSEA counts, GEO HRs) are not fixed here: genes and covariates change between runs. Check regressions against the committed outputs instead.
- **Refactors with no intended change in results:** rerun the affected scripts and check that `git status` / `git diff --stat output/tables/` shows no changed CSVs. Any change means the refactor altered results.
- **Intended changes** (covariates, genes, settings): commit the new outputs with the config change, so they become the next baseline.
- **Cohort-level invariants** (do not depend on the gene):
  - TCGA survival counts: the table above, enforced by the manifest check in 01.
  - Each gene's High/Low sizes: `merge_QC` in 02b and the median QC (±1 patient).
  - 05: the PCA group panel colours must match each gene's High/Low n; every cancer runs in its own process with no failures or warnings.
  - 06 GSE39582 (colon only, so it validates COAD; READ has no external cohort yet): 585 samples → 566 tumors (19 "Non Tumoral" excluded). Stage 0 (4 tumors) → NA. RFS excludes stage IV (half of them have `rfs.delay` = 0). OS n 562 (191 events); RFS n 497 (140 events); MMR known for 519 (dMMR 75).

## Conventions
- Code comments and console messages are in Korean. Table and figure labels are in English.
- On data inconsistencies, `stop()` and name the patient IDs, plus the `cfg` entry to fix when there is one. Don't silently drop or coerce.
- Output files are named `<cancer>_<analysis>[_<GENE>]`, or `<GSE>_…` for 06. Write them through `save_table()`/`save_df_table()`/`save_fig()`, which add the prefix; hand-built paths must use `out_file()`. File names in the docs omit the prefix.
- Call flextable functions as `flextable::<fn>`: survminer attaches ggpubr, whose `font()` masks `flextable::font()`.
- Keep flextable captions as plain text; markdown `**` shows up literally in the .docx.
- Combine the KM plot and risk table with patchwork (`p$plot / p$table`); wrap survminer plotting in `suppressMessages()` (see docs/notes.md).
- Build PDFs with `cairo_pdf` so that — ≤ – render correctly.
- Write CIs that can be negative (RMST) as `(lo, hi)`, not `lo–hi`.
- After changing `R/utils.R` or covariates in `config.R`, rerun 01 → 02b → 03 → 04. The manifest will stop 01 if the counts shifted.

## Reference docs (read only when the task needs them)
- docs/data_flow.md — `config.R` key groups, clinical load/recode steps, name standardization, gene CSVs, 05 inputs. Read before touching import/merge code.
- docs/covariates.md — full `cfg` covariate/strata/level table, what is unadjusted by design, how to change covariates. Read before changing any model setting.
- docs/methods.md — Cox/EPV/strata/PH/RMST/BH, GSEA design and parameters, Table 1 tests, 05 figure selection, 06 GEO validation. Read before changing 03–06 statistics or writing the paper's Methods.
- docs/outputs.md — every output file, output helpers, `_stale` behavior and prefixes. Read before adding or renaming outputs.
- docs/adding_cancer.md — procedure for a new cancer or an updated .sav. Read before adding a cancer or when the user supplies a new .sav.
- docs/notes.md — library quirks and history (PCA bug, msigdbr 26.x, survminer warning, 05 memory, untested GDC clinical path).
