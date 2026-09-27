# Covariates

Every setting below lives in `config.R` (the no-hardcoding / audit-grep rule is in CLAUDE.md). Add new settings to this table.

| Key | Default | Used in | Controls |
|---|---|---|---|
| `exclude_neoadjuvant` | `TRUE` | 01, 02b (→ 03–05) | Drop neoadjuvant-treated patients (applied at the end of 02b; see CLAUDE.md) |
| `analysis_na_values` | `residual_tumor = "RX"` | utils `recode_clinical` (01) → 02b–05 | Values treated as missing in every analysis (RX = residual tumor cannot be assessed); they appear in Table 1's Unknown row and are excluded from tests |
| `schema$neoadjuvant`, `schema$neoadjuvant_yes` | `"history_neoadjuvant_treatment"`, `"Yes"` | utils `recode_clinical` (01) | Which column and which value (case-insensitive) set `neoadjuvant_flag` |
| `group_levels` | `c("Low", "High")` | utils, 01–05 | Exposure levels: the 1st is the **reference**, the 2nd is the comparison (above median). Sets the 0/1 mapping (1 = 2nd level), factor order, Cox/KM reference, RMST and time-split arm, DESeq2 contrast and apeglm coefficient, and plot labels. |
| `reference_levels` | `gender = "Male"`, `pathologic_stage = "I"`, `stage = "I"`, `pathologic_stage_12_34 = "I–II"` | utils `apply_reference_levels` (01) → 03, 04 | Reference level of each categorical covariate |
| `cox_covariates` | `c("age", "gender", "stage")` | 04 | Multivariable Cox adjustment; also the adjustment in the strata(TSS) and time-split models, and the caption text |
| `cox_uni_covariates` | `c("age", "gender", "stage", "pathologic_stage_12_34")` | 04 | Covariates in the univariable Cox table |
| `covariate_scale` | `age = list(by = 10, label = "Age (per 10 years)")` | 04 (via utils `add_scaled_terms`) | Unit of continuous covariates in Cox; the model term becomes `age_per10` |
| `stage_full`, `stage_collapsed` | `"stage"`, `"pathologic_stage_12_34"` | 04; utils (creates the collapsed column if it's missing) | Stage term that is swapped for the collapsed one when EPV < `epv_min` |
| `epv_min` | `10` | 04 | EPV threshold for collapsing stage, then for the exploratory flag |
| `min_events` | `10` | 04 | Genes with fewer events are skipped |
| `strata_var` | `"tss"` | 03 (`<strata_var>_by_group_*.csv`), 04 (`strata()` sensitivity) | Stratification / site variable |
| `tss_barcode_pos` | `c(6, 7)` | utils (01) | Characters of `sample_id` used for `tss` |
| `collapse_small_levels`, `collapse_other_label` | `c(tss = 10)`, `"Other"` | 04 (strata), 05 (design), via utils `collapse_small` | Levels with fewer patients than this are merged into "Other" |
| `gsea_design_covariates` | `c("tss")` | 05 | Main DESeq2 design `~ <covariates> + group`; the first one also colours the 2nd PCA panel |
| `gsea_sensitivity_covariates` | `character()` | 05 | Sensitivity design (default `~ group`), which defines `robust_no_covariate` |
| `table1_vars`, `table1_continuous` | 15 clinical variables; `"age"` | 03 | Table 1 rows; which of them are summarized as median [IQR] with Wilcoxon |
| `value_labels_default`, `label_rules_default` (+ `value_labels$<C>`, `label_rules$<C>`) | stage I–IV, stage I–II/III–IV, CEA ≤5/>5, age ≤65/≥66 | utils (01) | Code → label mapping and consistency cutoffs for derived covariates |
| `drop_columns` | `"rock2g_01"` | utils (01) | Duplicate columns removed on import |
| `geo_dir` | `"data/geo"` | 06 | GEOquery download and eset cache (`<GSE>_<GPL>_eset.rds`) |
| `geo_probe_rule` | `"max_mean"` | 06 | Which probe of the primary gene is used: `max_mean`, `max_iqr` or `max_sd` (computed over tumors) |
| `geo_na_values` | `"", "N/A", "NA", "ND", …` | 06 | Raw pheno strings treated as missing |
| `geo_datasets$<GSE>$fields` | GSE39582: 14 fields | 06 | Canonical name → raw `characteristics_ch1` key. Required: `sample_type`, `age`, `gender`, `stage` and every endpoint's time/event. 06 prints all raw keys first. |
| `geo_datasets$<GSE>$exclude_sample_type` | `"Non Tumoral"` | 06 | Samples dropped before the median split |
| `geo_datasets$<GSE>$stage_na_values` | `"0"` | 06 | Stage codes set to NA (listed with IDs). Other codes must be in `value_labels_default$pathologic_stage` of the `validates` cohort. |
| `geo_datasets$<GSE>$endpoints` | OS (primary), RFS | 06 | time/event fields, `time_unit`, `event_value`, label, `exclude_stage` (RFS: stage IV) |
| `geo_datasets$<GSE>$factor_levels` | `mmr_status = c("pMMR", "dMMR")` | 06 | Allowed values and order (1st = reference). An unlisted value stops the run. |
| `geo_datasets$<GSE>$extra_covariate` | `"mmr_status"` | 06 | Added to `cox_covariates` in the "+ MMR" model; also the group × covariate LRT |
| `geo_datasets$<GSE>$subgroup` | `mmr_status == "pMMR"` | 06 | KM and Cox repeated in this subgroup |
| `geo_datasets$<GSE>$validates` | `"COAD"` | 06 | TCGA cohort for the stage labels and the HR direction comparison |

**Unadjusted by design:** KM and log-rank, RMST, univariable Cox, Table 1 tests, and the `<strata_var>` × group crosstab. 06 uses the same `cox_covariates`, `covariate_scale`, `reference_levels`, `group_levels`, `epv_min`, `min_events`, `km_times`, `ph_split_months` and `rmst_tau` as 04, so TCGA and GEO models match.

**Changing covariates:** edit only `config.R`, then rerun (see CLAUDE.md Conventions; the manifest guards patient counts). Confirm that the forest plots and `cox_multi_*` show the intended terms and reference rows. Changing `group_levels` also changes which group the HRs and log2FC are relative to.
