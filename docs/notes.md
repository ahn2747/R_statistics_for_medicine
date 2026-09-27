# Library and history notes

- **PCA bug (fixed):** with two `intgroup` variables, `DESeq2::plotPCA()` overwrites the `group` column with their interaction ("Low:AF"), so the group panel came out grey. `pca_plot()` passes only `intgroup = "group"` and adds the second colour variable from `colData`. PCA figures made before this fix were wrong.
- **msigdbr 26.x** uses `collection` / `subcollection` (not `category`), and its data download on first use. KEGG uses `CP:KEGG_LEGACY`, falling back to `CP:KEGG_MEDICUS`.
- **survminer 0.5.2** with ggplot2 4.x prints "Ignoring unknown labels". It is harmless, which is why survminer plotting is wrapped in `suppressMessages()`.
- **05 memory and speed:** DESeq2 runs twice, plus apeglm and GO:BP fgsea, taking about 10–15 minutes for COAD. A 4-worker run of both cancers in one process was killed for low memory, hence one cancer per process and `n_cores` = 2. 05 frees `dds`/`vst` and calls `gc()` after each gene and each cancer.
- **`load_gdc_clinical()`:** used when `<C>.sav` doesn't exist, building the canonical columns from GDC. This path is untested.
- **Manifest history:** `db_manifest.json` replaced the old hardcoded `expected_surv`.
