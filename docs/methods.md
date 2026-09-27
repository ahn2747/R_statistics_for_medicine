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
