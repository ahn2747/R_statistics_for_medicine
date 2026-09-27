# Statistical methods

Rules agreed with the user; don't change silently. Setting names are in docs/covariates.md. Written so it can feed the paper's Methods section.

## Site (TSS) handling
- Sites (`tss`) with fewer than `collapse_small_levels["tss"]` (10) patients are grouped as "Other", both for the Cox `strata(tss)` sensitivity analysis and the DESeq2 design.

## Table 1
- Wilcoxon for continuous variables; chi-square, switching to Fisher's exact when any expected cell is < 5. Missing values are shown as a row and excluded from the tests.

## Survival (04)
- **Cox models:**
  - The multivariable model is the gene group + `cox_covariates` (age per 10 years, gender, stage), as a complete-case analysis.
  - EPV = events / number of parameters. If EPV < `epv_min`, `stage_12_34` replaces stage I–IV. If EPV is still low, the result is `exploratory = TRUE`; all READ models are (EPV about 5.8).
  - Sensitivity analysis: `strata(tss)`. The expression group is strongly associated with TSS for many genes (`tss_by_group_*.csv`), so report this analysis.
- **PH violations:**
  - Time-split Cox (0–24 / >24 months, adjusted) is run only when the gene **group term** has cox.zph p < 0.05. Currently that is READ TIMP1 and READ APC.
  - COAD TIMP1 violates PH only through its stage term.
  - RMST (High − Low, tau = 60) is reported for every gene.
- **Multiple testing:** BH q-values are computed across the exploratory genes only. The primary gene reports its raw p, and its q is shown as "–".

## GSEA (05)
- Use the existing `<gene>_group` from the merged data (the same patients and grouping as survival). Never re-split.
- Main design: `~ tss + group` (TSS collapsed as above; groups are unbalanced across sites). Contrast: `group_levels[2]` vs `group_levels[1]` (High vs Low).
- Pre-filter: keep genes with count ≥ 10 in at least as many samples as the smaller group.
- Rank by the DESeq2 Wald stat, remap symbols with the `.chip` file, and average duplicates.
- **Gene set database: MSigDB v2026.1.Hs** (human). Both the gene sets and the symbol remapping use this release:
  - Gene sets come from msigdbr 26.1.1, which downloads `msigdb.2026.1.zip` (Zenodo record 18968178) into `tools::R_user_dir("msigdbr", "cache")`. Upgrading msigdbr can change the release, so check the cache file names (`msigdb.<release>.Hs.*.rds`) after any upgrade.
  - The remapping file is `database/Human_Gene_Symbol_with_Remapping_MSigDB.v2026.1.Hs.chip`. Keep its release the same as msigdbr's.
  - Collections and set counts as loaded (before the size filter): Hallmark H 50, Reactome C2:CP:REACTOME 1,839, KEGG **C2:CP:KEGG_LEGACY** 186 (not MEDICUS), GO:BP C5:GO:BP 7,538, C8 866.
- Software: R 4.6.1, DESeq2 1.52.0, apeglm 1.34.0, fgsea 1.38.0, msigdbr 26.1.1, TCGAbiolinks 2.40.0.
- `fgseaMultilevel` with minSize 15, maxSize 500, eps 0 and a fixed seed.
- GO:BP runs `collapsePathways` on the top 300 significant pathways.
- `robust_no_covariate` = same NES sign and padj < 0.05 in the `~ group` run.
- QC: Spearman correlation between GDC log2 TPM and the .sav expression (warns if r < 0.8). The gene itself must rank in the top 1% up in High, otherwise a warning is printed.

## 05 figure and summary selection
- Enrichment plots cover:
  - the top `gsea_enrichment_top` Hallmark pathways by padj (significant only)
  - every `gsea_highlight` pathway
  - the top C8 cell-type sets matching `gsea_c8_pattern` (default B cells / plasma cells, significant only)
- `gsea_summary` lists the Hallmark top 5 up and top 5 down by NES, and the top 10 by padj for the other collections.
- apeglm optimizer warnings are suppressed, and the number suppressed is printed.

## External validation in GEO (06)
- Validates only `cfg$primary_gene`, in each `cfg$geo_datasets` entry. GSE39582 (Marisa et al. 2013; CIT; Affymetrix HG-U133 Plus 2.0, RMA + ComBat as deposited) is colon cancer only, so it validates COAD.
- **Expression:** the series matrix as deposited. If max > 100 it is log2-transformed, otherwise used as is (GSE39582 is already log2). Every probe annotated to the gene in the GPL `Gene Symbol` column (splitting on `///`) is reported. One probe is chosen by `geo_probe_rule` (default highest mean over tumors), with no hardcoded probe IDs.
- **Samples:** non-tumor samples are excluded first. High/Low = above/at or below the median of the chosen probe over **all tumors** (`median_split()`), before any endpoint-specific exclusion.
- **Endpoints:** OS is primary and RFS secondary. RFS excludes stage IV, which has no disease-free interval. Patients with only one of time/event are excluded from that endpoint and listed by ID.
- **Models (per endpoint, all tumors):**
  - univariable (group)
  - univariable continuous (per 1 log2 unit)
  - multivariable with `cox_covariates`, using the same EPV rule as 04 (`build_multi()`)
  - the same multivariable model restricted to patients with known extra covariate
  - multivariable + `extra_covariate` (MMR status), with the group × MMR likelihood-ratio interaction test
  - Comparing the last two on the same patients isolates the effect of MMR adjustment. This answers the MSI/immune-confounding question that TCGA can't address yet (no MSI result).
- **Subgroup:** KM, univariable and multivariable models are repeated in pMMR.
- **PH, time-split and RMST:** these follow the 04 rules. cox.zph runs for every model. The time-split Cox (0–24 / >24 months) runs on the multivariable model only when the gene term violates PH. RMST (High − Low, tau 60) is always reported. Time-0 patients are kept (the survSplit origin is shifted below 0).
- **Direction check:** each GEO HR is compared with the sign of the TCGA multivariable HR from `output/tables/<validates>/<validates>_survival_summary_raw.csv`.
