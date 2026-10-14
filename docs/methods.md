# Statistical methods

Rules agreed with the user; don't change silently. Setting names are in docs/covariates.md. Written so it can feed the paper's Methods section.

## Site (TSS) handling
- Sites (`tss`) with fewer than `collapse_small_levels["tss"]` (10) patients are grouped as "Other", both for the Cox `strata(tss)` sensitivity analysis and the DESeq2 design.

## Table 1
- Wilcoxon for continuous variables; chi-square, switching to Fisher's exact when any expected cell is < 5. Missing values are shown as a row and excluded from the tests.
- Correlation (03, `correlation_<GENE>.csv`): Pearson r between log2(expression + 1) of the gene and `cfg$cor_vars$genes`, and the raw values of `cfg$cor_vars$clinical` (age, pretreatment CEA), using pairwise complete observations in the gene's Table 1 patients. Exploratory; no multiple-comparison adjustment.

## Survival (04)
- **Cox models:**
  - The multivariable model is the gene term + `cox_covariates` (age group ≤65/≥66, sex, stage I–II/III–IV, pretreatment CEA ≤5/>5 ng/mL), as a complete-case analysis. Before 2026-10-04 the config listed `ageG`, which isn't a canonical column, so age was silently left out of every univariable and multivariable model; 04/06 now stop on a covariate that isn't in the data.
  - **Gene term** (`cfg$cox_gene_term`): `"group"` (default) = High vs Low; `"continuous"` = `<gene>_expression_log2`, HR per 1 log2 unit. It applies to the multivariable, `strata(tss)` and time-split models, to the univariable HR in `survival_summary` and `forest_genes`, and to the PH check that triggers the time-split model. KM, log-rank and RMST always compare High vs Low. The `cox_uni` table always lists both the group and continuous univariable rows. `survival_summary_raw.csv` records the setting in `cox_gene_term`.
  - EPV = events / number of parameters. If EPV < `epv_min`, `stage_12_34` replaces stage I–IV. If EPV is still low, the result is `exploratory = TRUE`; all READ models are (EPV about 5.8).
  - Sensitivity analysis: `strata(tss)`. The expression group is strongly associated with TSS for many genes (`tss_by_group_*.csv`), so report this analysis.
- **PH violations:**
  - Time-split Cox (0–24 / >24 months, adjusted) is run only when the **gene term** (`cox_gene_term`) has cox.zph p < 0.05. With the group term, that was READ TIMP1 and READ APC.
  - COAD TIMP1 violates PH only through its stage term.
  - RMST (High − Low, tau = 60) is reported for every gene.
- **Multiple testing:** BH q-values are computed across the exploratory genes only. The primary gene reports its raw p, and its q is shown as "–".

## RNA-seq data source (02a)
Values are from `data/processed/<C>_gdc_provenance.json` (see docs/outputs.md). Rerun `02a --provenance-only` and update this section if the data are re-downloaded.
- GDC Data Release 46.0 (August 10, 2026), recorded in `metadata(se)` by GDCprepare at download. Files were downloaded on 2026-09-27 (cache mtimes: COAD 11:11–11:18, READ 11:23–11:25 KST). All cached files (COAD 524, READ 177) match the GDC md5 checksums.
- Gene model: GENCODE v36 (TSV header). The genome build is not stored in the se; TCGAbiolinks reports "hg38" for this harmonized query.
- Assays: `unstranded` counts → DESeq2; `tpm_unstrand` → log2(TPM + 1) for QC and `genes_from_gdc`.
- Samples:
  - COAD: 524 files → 481 primary tumor (01) samples → 458 patients. 13 patients had several tumor aliquots; 23 aliquots were dropped, keeping the first per patient.
  - READ: 177 files → 166 primary tumor samples → 166 patients (no duplicates).
- Genes: 60,660 features → 59,427 unique symbols. For duplicated symbols, the one with the highest mean TPM is kept.
- Software: TCGAbiolinks 2.40.0, SummarizedExperiment 1.42.0, R 4.6.1.

Methods text:

> RNA-seq gene expression data (STAR - Counts workflow; data category "Transcriptome Profiling", data type "Gene Expression Quantification") for TCGA-COAD and TCGA-READ were obtained from the NCI Genomic Data Commons (GDC Data Release 46.0, August 10, 2026) using TCGAbiolinks (v2.40.0) on September 27, 2026. Gene-level quantification was based on the GENCODE v36 annotation on the GRCh38 reference, as provided by the GDC harmonized pipeline. Only primary tumor samples (TCGA sample type code 01) were retained; when a patient had more than one primary tumor aliquot, the first aliquot was kept (COAD, 13 patients). Ensembl gene identifiers were collapsed to HGNC symbols by retaining, for each symbol, the gene with the highest mean TPM. Unstranded read counts were used for differential expression analysis, and log2(TPM + 1) values were used for quality control. The final matrices comprised 59,427 genes for 458 COAD and 166 READ patients.

## GSEA (05)
- Use the existing `<gene>_group` from the merged data (the same patients and grouping as survival). Never re-split.
- Main design: `~ tss + group` (TSS collapsed as above; groups are unbalanced across sites). Contrast: `group_levels[2]` vs `group_levels[1]` (High vs Low).
- Pre-filter: keep genes with count ≥ 10 in at least as many samples as the smaller group.
- Rank by the DESeq2 Wald stat, remap symbols with the `.chip` file, and average duplicates.
- **Gene set database: MSigDB v2026.1.Hs** (human). Both the gene sets and the symbol remapping use this release:
  - Gene sets come from msigdbr 26.1.1, which downloads `msigdb.2026.1.zip` (Zenodo record 18968178) into `tools::R_user_dir("msigdbr", "cache")`. Upgrading msigdbr can change the release, so check the cache file names (`msigdb.<release>.Hs.*.rds`) after any upgrade.
  - The remapping file is `database/Human_Gene_Symbol_with_Remapping_MSigDB.v2026.1.Hs.chip`. Keep its release the same as msigdbr's.
  - **Collection: Hallmark only** (H, 50 sets as loaded, before the size filter). Collections are set in `cfg$gsea_collections`. Reactome, KEGG, GO:BP and C8 are no longer run; example specs for them are commented out in config.R.
- Software: R 4.6.1, DESeq2 1.52.0, apeglm 1.34.0, fgsea 1.38.0, msigdbr 26.1.1, TCGAbiolinks 2.40.0.
- `fgseaMultilevel` with minSize 15, maxSize 500, eps 0 and a fixed seed.
- If GO:BP is added back to `cfg$gsea_collections`, it runs `collapsePathways` on the top 300 significant pathways.
- `robust_no_covariate` = same NES sign and padj < 0.05 in the `~ group` run.
- QC: Spearman correlation between GDC log2 TPM and the .sav expression (warns if r < 0.8). The gene itself must rank in the top 1% up in High, otherwise a warning is printed.

## 05 figure and summary selection
- Enrichment plots cover:
  - the top `gsea_enrichment_top` Hallmark pathways by padj (significant only)
  - every `gsea_highlight` pathway
  - only if C8 is added back to `cfg$gsea_collections`: the top C8 cell-type sets matching `gsea_c8_pattern` (significant only)
- `gsea_summary` lists the significant Hallmark top 5 up and top 5 down by NES. Any other collection added to `cfg$gsea_collections` gets its top 10 by padj.
- apeglm optimizer warnings are suppressed, and the number suppressed is printed.

## External validation in GEO (06)
- Validates only `cfg$primary_gene`, in each `cfg$geo_datasets` entry. GSE39582 (Marisa et al. 2013; CIT; Affymetrix HG-U133 Plus 2.0, RMA + ComBat as deposited) is colon cancer only, so it validates COAD.
- **Expression:** the series matrix as deposited. If max > 100 it is log2-transformed, otherwise used as is (GSE39582 is already log2). Every probe annotated to the gene in the GPL `Gene Symbol` column (splitting on `///`) is reported. One probe is chosen by `geo_probe_rule` (default highest mean over tumors), with no hardcoded probe IDs.
- **Samples:** non-tumor samples are excluded first. High/Low = above/at or below the median of the chosen probe over **all tumors** (`median_split()`), before any endpoint-specific exclusion.
- **Endpoints:** OS is primary and RFS secondary. RFS excludes stage IV, which has no disease-free interval. Patients with only one of time/event are excluded from that endpoint and listed by ID.
- **Models (per endpoint, all tumors):**
  - univariable (group)
  - univariable continuous (per 1 log2 unit)
  - univariable, one model per clinical covariate (`cox_uni_covariates` minus `cox_covariates_drop`, plus `extra_covariate` in all tumors; each on its own complete cases, as in 04). These appear only in the `cox_<EP>` tables, not in the summary, forest plot or key results
  - multivariable with the gene term (`cox_gene_term`, as in 04) + `cox_covariates` minus `cox_covariates_drop` (GSE39582: no CEA), using the same EPV rule as 04 (`build_multi()`). `age_g` is derived from GEO age with the TCGA rule (age > 65 → ≥66).
  - the same multivariable model restricted to patients with known extra covariate
  - multivariable + `extra_covariate` (MMR status), with the gene term × MMR likelihood-ratio interaction test
  - The TCGA direction comparison uses 04's `multi_hr`. It is skipped with a warning if 04 ran with a different `cox_gene_term`.
  - Comparing the last two on the same patients isolates the effect of MMR adjustment. This answers the MSI/immune-confounding question that TCGA can't address yet (no MSI result).
- **Subgroup:** KM, univariable and multivariable models are repeated in pMMR.
- **PH, time-split and RMST:** these follow the 04 rules. cox.zph runs for every model. The time-split Cox (0–24 / >24 months) runs on the multivariable model only when the gene term violates PH. RMST (High − Low, tau 60) is always reported. Time-0 patients are kept (the survSplit origin is shifted below 0).
- **Direction check:** each GEO HR is compared with the sign of the TCGA multivariable HR from `output/tables/<validates>/<validates>_survival_summary_raw.csv`.

## Signature scores (11; separate track)
- **Definition:** `cfg$signatures$<SIG>` (up and/or down gene sets, required `source`) from external literature/MSigDB, not derived from this TCGA data. Exact symbol matching, no alias conversion. The track never touches `merged.rds`, `detect_genes()`, 04's BH set, 05, 90 or run_gene.R.
- **Input:** 02a log2(TPM + 1) (`database/TCGA_<C>_RNAseq_Expression.csv`), each cancer separately. Score samples = `_clinical_all.rds` patients with expression (before the neoadjuvant exclusion; COAD 454, READ 165).
- **Filter:** genes with log2(TPM + 1) ≥ 1 in ≥ 20% of score samples and non-zero variance (`signature_expr_filter$tcga`; COAD/READ ≈ 19,400 genes). The filtered matrix is the common background for all four scores and the random gene sets.
- **Matching:** used (matched and passing the filter) / defined must be ≥ `signature_min_match` (0.8) per direction and overall, with ≥ `signature_min_genes` (3) per direction; otherwise 11 stops and lists the genes.
- **Scores** (per set, then score = score(up) − score(down)):
  - z: per-gene z over score samples, averaged
  - PC1: `prcomp(scale. = TRUE)`, sign set so that PC1 correlates positively with the z-score
  - ssGSEA: GSVA `ssgseaParam` (alpha 0.25, `normalize = FALSE`), up and down in one call
  - singscore: `rankGenes` + `simpleScore(centerScore = TRUE)` per set (equal to `simpleScore(up, down)$TotalScore`)
  - Main method fixed before results: `signature_primary_method` = ssGSEA.
- **QC (no outcome used):** Cronbach's alpha (`psych::alpha`, down genes sign-flipped) with per-gene r.drop and alpha if dropped; PC1 variance explained per set vs `signature_random_n` (1000) random gene sets of the same size from the filtered background (seed `cfg$seed`), as a percentile and empirical p; Spearman correlations between the four scores.

## Signature survival (12; TCGA OS)
- **Cohort:** scores joined with `_clinical.rds` (the analytic file with `cfg$exclude_neoadjuvant` applied; 12 stops if it doesn't match the current setting), patients with OS time/status. COAD 438 (97 events; TCGA-A6-2670 has no RNA-seq), READ 157 (25), with the current `exclude_neoadjuvant = FALSE`.
- **High/Low:** `median_split()` of each score over **all scored patients** (before exclusion), like the .sav groups and 06.
- **Continuous term:** `score_sd` = score standardized within the survival cohort; HR per 1 SD, regardless of `cox_gene_term`.
- **Models per method** (all four reported; the primary method `signature_primary_method` is listed first): KM + log-rank, RMST (High − Low, tau 60), univariable Cox (High vs Low; per 1 SD), multivariable Cox (per 1 SD + `cox_covariates`, `build_multi()` EPV rule, complete case), multivariable + `strata(tss)` (sites < 10 collapsed), cox.zph on every unstratified model, time-split Cox (0–24 / >24 months) only when the score term violates PH. All computed with R/survival.R functions, unchanged.
- **Multiplicity:** never pooled with the gene-level BH in 04. With several signatures and a full run (no SIG argument), `q_across_signatures` = BH across signatures for the primary method, per model (and for the log-rank p); otherwise NA.
