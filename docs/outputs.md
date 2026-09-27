# Outputs

## Output helpers (`R/utils.R`)
- All of these take an optional `subdir` (e.g. `"gsea/MS4A1"`) under `output/<type>/<C>/`.
- `save_table()` (CSV, UTF-8 BOM) and `save_df_table()` (.docx + .csv)
- `save_fig()` (cairo PDF + 300-dpi LZW TIFF)
- `fmt_hr()` and `fmt_p()` (2 decimals; `<0.001`)
- `move_stale_outputs()`: gene-specific files that this run didn't recreate are moved to `output/*/<C>/_stale/`, never deleted. It uses the prefixes each script owns:
  - 03: `table1_`, `tss_by_group_`
  - 04: `km_`, `forest_multi_`, `zph_`, `cox_multi_`
- `move_stale_dirs()`: gene subfolders under `gsea/` that aren't in the current run are moved to `gsea/_stale/`.

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
