# CLAUDE.md

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
05_gsea.R [COAD]        # DESeq2 High vs Low (~ tss + group) + fgsea on MSigDB; needs 02a outputs; run one cancer per process
06_external_geo.R [GSE] # GEO validation of cfg$primary_gene (cfg$geo_datasets; needs GEOquery, i.e. 00_setup.R without --core; reads 04's survival_summary_raw.csv)
```
- `07_timer` and `08_proteomics` are optional and not written yet.
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
- **GEO datasets** (06) are config-only: `cfg$geo_datasets$<GSE>` maps canonical names to raw `characteristics_ch1` keys, and sets the platform, sample types to exclude, endpoints, factor levels, extra covariate and subgroup. Don't hardcode GSE IDs, probe IDs or pheno keys.
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

## Reference results (regression checks; checked against independent code)
- COAD MS4A1: log-rank p 0.018; multivariable HR 0.55 (0.35–0.85), p 0.007; strata(TSS) HR 0.50; RMST difference 4.88 months
- READ MS4A1: log-rank p 0.008; RMST difference 9.11 months
- 05 GSEA (full run, each cancer in its own process, no failures or warnings):
  - MS4A1 ranks 1st of all tested genes in both cancers (COAD 17,748, log2FC 4.09; READ 18,089, log2FC 3.64).
  - Significant pathways (padj < 0.05), with the number still significant in the no-covariate run in brackets:

    | | Hallmark | Reactome | KEGG | GO:BP | C8 |
    |---|---|---|---|---|---|
    | COAD | 30 (29) | 412 (373) | 72 (65) | 1288 (1243) | 483 (471) |
    | READ | 28 (27) | 336 (315) | 54 (49) | 790 (742) | 384 (354) |

  - Up in High: allograft rejection, IFN-γ response, inflammatory response.
  - Up in Low: MYC targets, oxidative phosphorylation, E2F, G2M.
  - The PCA group panel is coloured correctly (COAD 218/218, READ 79/78).
- 06 GEO, GSE39582 (colon only, so it validates COAD; READ has no external cohort yet):
  - 585 samples → 566 tumors (19 "Non Tumoral" excluded). Stage 0 (4 tumors) → NA. RFS excludes stage IV (half of them have `rfs.delay` = 0).
  - MS4A1 has 4 GPL570 probes (r 0.78–0.95). `max_mean` picks `228592_at` (mean 4.29, IQR 1.57); median split 283/283.
  - OS n 562 (191 events); RFS n 497 (140 events); MMR known for 519 (dMMR 75).
  - **Does not replicate TCGA.** OS: univariable HR 1.04 (0.78–1.38), p 0.80; multivariable 1.14 (0.85–1.52), p 0.38; + MMR 1.22 (0.90–1.66), p 0.19; pMMR multivariable 1.27 (0.92–1.74). RFS: multivariable 1.11 (0.80–1.55). No probe (or the probe mean) gives HR < 1.
  - PH is violated for the gene term (OS univariable zph p 0.005). Time-split OS: 0–24 mo HR 2.50 (1.54–4.07), >24 mo HR 0.68 (0.46–0.99). OS RMST (60 mo) High − Low −2.58 (−5.57, 0.41).
  - Group × MMR interaction: OS p 0.36, RFS p 0.20.

## Conventions
- Code comments and console messages are in Korean. Table and figure labels are in English.
- On data inconsistencies, `stop()` and name the patient IDs, plus the `cfg` entry to fix when there is one. Don't silently drop or coerce.
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
