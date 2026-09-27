# Covariates

Every setting below lives in `config.R` (the no-hardcoding / audit-grep rule is in CLAUDE.md). Add new settings to this table.

| Key | Default | Used in | Controls |
|---|---|---|---|
| `exclude_neoadjuvant` | `TRUE` | 01, 02b (→ 03–05) | Drop neoadjuvant-treated patients (applied at the end of 02b; see CLAUDE.md) |
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

**Unadjusted by design:** KM and log-rank, RMST, univariable Cox, Table 1 tests, and the `<strata_var>` × group crosstab.

**Changing covariates:** edit only `config.R`, then rerun (see CLAUDE.md Conventions; the manifest guards patient counts). Confirm that the forest plots and `cox_multi_*` show the intended terms and reference rows. Changing `group_levels` also changes which group the HRs and log2FC are relative to.
