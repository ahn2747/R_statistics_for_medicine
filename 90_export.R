# =============================================================
# 90_export.R [GENE] [--dry-run]
# 이미 실행된 03–06 결과를 학생연구자료 폴더로 내보내기
#   <cfg$export$dest_root>/<NAME>(COAD&READ)/
#     firstline/, figure/            03·04 CSV, TIFF (cfg$export$layout$tcga)
#     gsea/data/, gsea/figure/       05 (gsea/<GENE>/)
#     geo/<GSE>_<validates>/data|figure/   06 (이 유전자의 결과일 때만)
#     <NAME>_tables.docx             양식(cfg$export$template)을 채운 표 + KM 그림
#     <NAME>_key_results.csv         핵심 수치 요약
#     _export_manifest.csv           복사 기록 (md5, git 커밋, 원본 → 대상 이름)
#     _old/<YYYYMMDD_HHMM>/          내용이 바뀐 기존 파일 백업 (삭제하지 않음)
# GENE 생략 = cfg$primary_gene. --dry-run = 계획만 출력 (파일 변경 없음)
# NAME = GENE, 또는 cfg$export$use_custom_filename = TRUE면 cfg$export$custom_filename
#   (복사 파일은 layout 패턴의 {G} 자리만 NAME으로 바뀜. 결과 선택·docx 내용은 실제 GENE 기준)
# 같은 대상 이름 파일: md5가 같으면 건너뜀, 다르면 기존 파일을 _old/로 옮기고 복사
# =============================================================

source("config.R")
source("R/utils.R")
source("R/export.R")

args <- commandArgs(trailingOnly = TRUE)
dry  <- "--dry-run" %in% args
pos  <- args[!startsWith(args, "--")]
gene <- toupper(if (length(pos)) pos[1] else cfg$primary_gene)
ex   <- cfg$export
ename <- export_name(gene)
dest <- export_dest(ename)
ts   <- format(Sys.time(), "%Y%m%d_%H%M")
cat("유전자", if (ename != gene) paste0(gene, " (내보내기 이름 ", ename, ")") else gene, "→", dest,
    if (dry) "(dry-run)" else "", "\n")

plan    <- data.frame(source = character(), rel = character(), stringsAsFactors = FALSE)
skipped <- data.frame(what = character(), reason = character(), stringsAsFactors = FALSE)
add_plan <- function(files, rel_dir, dest_base = basename(files)) {
  if (length(files)) plan <<- rbind(plan, data.frame(source = files, rel = file.path(rel_dir, dest_base),
                                                     stringsAsFactors = FALSE))
}
add_skip <- function(what, reason) {
  skipped <<- rbind(skipped, data.frame(what = what, reason = reason, stringsAsFactors = FALSE))
  warning(what, ": ", reason, call. = FALSE)
}
sub_type <- function(sub) if (sub == "figure") "figures" else "tables"
sub_ext  <- function(sub) if (sub == "figure") ex$fig_ext else ex$table_ext
# TCGA 출력 패밀리 → 만드는 스크립트 (재실행 안내용)
producer <- function(p) if (grepl("^(table1_|correlation_|analysis_info_03)", p)) "03_table1.R" else "04_survival.R"

# ---- 1. 유전자 존재 + 2. 03·04 출력 + 3. 신선도 ----
missing <- character()
for (cancer in cfg$cancers) {
  merged <- processed_path(cancer, "merged.rds")
  if (!file.exists(merged)) stop(merged, " 없음 → 02b_merge_genes.R 먼저 실행")
  if (!gene_key(gene) %in% detect_genes(readRDS(merged))) {
    stop(cancer, ": ", gene, "이(가) merged 데이터에 없음 (", cfg$gene_dir,
         "에 gene CSV를 넣고 02b_merge_genes.R 실행)")
  }
  m_time <- file.mtime(merged)
  for (sub in names(ex$layout$tcga)) {
    tg <- layout_targets(ex$layout$tcga[[sub]], gene, ename)
    mt <- match_targets(src_dir(cancer, sub_type(sub)), cancer, tg, sub_ext(sub))
    for (p in tg$pat[!tg$pat %in% mt$pattern]) {
      missing <- c(missing, sprintf("%s_%s.%s → %s 재실행 필요", cancer, p, sub_ext(sub), producer(p)))
    }
    f <- mt$source
    old <- f[file.mtime(f) < m_time]
    if (length(old)) {
      warning(cancer, ": ", length(old), "개 출력이 ", basename(merged), "보다 오래됨 (재실행하지 않은 결과일 수 있음): ",
              paste(basename(old), collapse = ", "), call. = FALSE)
    }
    add_plan(f, sub, mt$dest_base)
  }
}
if (length(missing)) stop("03·04 출력 없음:\n  ", paste(unique(missing), collapse = "\n  "))

# ---- 4. GSEA (없으면 경고 후 건너뜀) ----
for (cancer in cfg$cancers) {
  gsub_dir <- file.path("gsea", gene)
  if (!dir.exists(src_dir(cancer, "tables", gsub_dir))) {
    stale <- dir.exists(src_dir(cancer, "tables", file.path("gsea", "_stale", gene)))
    add_skip(paste(cancer, "GSEA"), paste0(
      "output/tables/", cancer, "/gsea/", gene, "/ 없음 → Rscript 05_gsea.R ", cancer,
      if (gene != toupper(cfg$primary_gene)) " (primary가 아니면 cfg$gsea_genes 또는 MEDIN_GENE 지정)" else "",
      if (stale) paste0("; gsea/_stale/", gene, "/에 이전 결과가 있으나 자동으로 가져오지 않음") else ""))
    next
  }
  for (sub in names(ex$layout$gsea)) {
    tg <- layout_targets(ex$layout$gsea[[sub]], gene, ename)
    mt <- match_targets(src_dir(cancer, sub_type(sub), gsub_dir), cancer, tg, sub_ext(sub))
    if (!all(tg$pat %in% mt$pattern)) {
      warning(cancer, " GSEA: 없는 파일 패턴 ", paste(tg$pat[!tg$pat %in% mt$pattern], collapse = ", "), call. = FALSE)
    }
    add_plan(mt$source, file.path("gsea", sub), mt$dest_base)
  }
}

# ---- 5. GEO (이 유전자의 결과일 때만) ----
geo_ok <- character()
for (gse in names(cfg$geo_datasets)) {
  if (!dir.exists(src_dir("GEO", "tables", gse))) {
    add_skip(gse, "06 출력 없음 → 06_external_geo.R 실행")
    next
  }
  gr <- geo_result_gene(gse)
  if (is.na(gr$gene)) {
    add_skip(gse, "어느 유전자의 결과인지 확인 불가 (gene 열·probe_QC 없음) → 06_external_geo.R 재실행")
    next
  }
  if (gr$gene != gene) {
    add_skip(gse, paste0("06 결과는 ", gr$gene, "의 것 (", gr$source, ") → MEDIN_GENE=", gene,
                         "로 06_external_geo.R 재실행"))
    next
  }
  geo_ok <- c(geo_ok, gse)
  rel <- file.path("geo", paste0(gse, "_", cfg$geo_datasets[[gse]]$validates))
  for (sub in names(ex$layout$geo)) {
    tg <- layout_targets(ex$layout$geo[[sub]], gene, ename)
    mt <- match_targets(src_dir("GEO", sub_type(sub), gse), gse, tg, sub_ext(sub))
    add_plan(mt$source, file.path(rel, sub), mt$dest_base)
  }
}

# ---- 생성 파일: key_results, tables.docx (임시 폴더에 만든 뒤 같은 규칙으로 복사) ----
tmp <- file.path(tempdir(), "export")
dir.create(tmp, showWarnings = FALSE)
key <- do.call(rbind, c(lapply(cfg$cancers, key_rows_tcga, gene = gene),
                        lapply(geo_ok, key_rows_geo),
                        lapply(cfg$cancers, key_rows_gsea, gene = gene)))
key_f <- file.path(tmp, paste0(ename, "_key_results.csv"))
readr::write_excel_csv(key, key_f, na = "")
generated <- key_f
if (!is.null(ex$template) && file.exists(ex$template)) {
  docx_f <- file.path(tmp, paste0(ename, "_tables.docx"))
  build_tables_docx(gene, docx_f)
  generated <- c(generated, docx_f)
} else {
  add_skip("tables.docx", paste0("양식 없음: ", ex$template, " (cfg$export$template)"))
}
add_plan(generated, ".")
plan$rel <- sub("^\\./", "", plan$rel)
check_plan_collisions(plan)

# ---- 계획 ----
plan <- plan_copy(plan, dest)
cat("\n대상:", dest, "\n")
for (i in seq_len(nrow(plan))) cat(sprintf("  %-8s %s\n", plan$status[i], plan$rel[i]))
for (i in seq_len(nrow(skipped))) cat(sprintf("  %-8s %s — %s\n", "skip", skipped$what[i], skipped$reason[i]))
n_by <- table(factor(plan$status, levels = c("new", "replace", "same")))
cat(sprintf("\nnew %d / replace %d / same %d / skip %d\n", n_by[["new"]], n_by[["replace"]], n_by[["same"]],
            nrow(skipped)))

if (dry) {
  cat("dry-run: 파일 시스템 변경 없음\n")
  quit(save = "no", status = 0)
}

# ---- 실행 ----
old_dir <- file.path(dest, "_old", ts)
fail <- apply_copy(plan, dest, old_dir)
plan <- attr(fail, "plan")

changed <- any(plan$status %in% c("new", "replace"))
man_f <- file.path(dest, "_export_manifest.csv")
if (changed || !file.exists(man_f)) {
  gi <- git_info()
  ct <- tryCatch(read.csv(file.path(src_dir(cfg$cancers[1], "tables"),
                                    out_file(cfg$cancers[1], "survival_summary_raw.csv")))$cox_gene_term[1],
                 error = function(e) NA)
  man <- data.frame(
    source = ifelse(plan$source %in% generated, "generated", normalizePath(plan$source, winslash = "/")),
    dest = plan$rel, status = plan$status, md5 = plan$md5,
    size = file.size(plan$source), source_mtime = format(file.mtime(plan$source), "%Y-%m-%d %H:%M:%S"),
    exported_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    git_head = gi$head, git_dirty = gi$dirty, gene = gene, export_name = ename, cox_gene_term = ct %||% NA)
  if (file.exists(man_f)) {
    dir.create(old_dir, recursive = TRUE, showWarnings = FALSE)
    bak <- file.path(old_dir, basename(man_f))
    if (file.exists(bak)) bak <- paste0(bak, ".", format(Sys.time(), "%H%M%S"))
    file.rename(man_f, bak)
  }
  readr::write_excel_csv(man, man_f, na = "")
  if (gi$dirty %in% TRUE) warning("git 작업 트리에 커밋하지 않은 변경이 있음 (manifest git_dirty = TRUE)", call. = FALSE)
} else {
  cat("변경 없음 → manifest 유지\n")
}

if (length(fail)) {
  cat("\n복사 실패 (", length(fail), "개):\n", paste0("  ", fail, collapse = "\n"), "\n", sep = "")
  quit(save = "no", status = 1)
}
cat("\n완료:", dest, "\n")
