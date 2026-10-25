# =============================================================
# 20_paper_prep.R [C ...] [--refresh-gdc]
# 논문 분석용 데이터 준비 (암종마다 따로) → data/processed/<C>_paper.rds
#   입력: <C>_merged.rds (02b), <C>_sigscores_<SIG>.csv (11), database/TCGA_<C>_RNAseq_Expression.csv (02a)
#   1. signature 점수/그룹 (cfg$paper$scores): .sav 점수는 z 평균 재계산과 대조 (불일치 → 중단;
#      .sav 정의 = 결측 발현은 유전자 평균으로 대치 후 z 평균),
#      11의 점수는 점수 샘플 전체에서 median_split() (12와 같음)
#   2. 흡연 상태: GDC clinical (TCGAbiolinks) → <C>_gdc_clinical.csv 캐시 (--refresh-gdc로 다시 받기)
#      cfg$paper$smoking$map으로 Never/Former/Current. GDC pack-years와 .sav pack-years 대조 (불일치 → 중단)
#   3. tumor purity: ESTIMATE (tidyestimate) on log2(TPM+1) → <C>_purity_estimate.csv
#      tidyestimate는 RNA-seq에서 purity를 주지 않음 → ESTIMATE 원 공식 cos(0.6049872018 + 0.0001467884 × ESTIMATE score)
#      (Affymetrix 기준 공식; Methods에 명시)
# 결과: output/tables/<C>/paper/ <C>_paper_QC.csv, <C>_paper_QC_smoking.csv
# =============================================================

source("config.R")
source("R/utils.R")
source("R/paper.R")

args <- commandArgs(trailingOnly = TRUE)
refresh <- "--refresh-gdc" %in% args
cancers <- if (any(toupper(args) %in% cfg$cancers)) toupper(args[toupper(args) %in% cfg$cancers]) else cfg$cancers
pc <- cfg$paper
sm <- pc$smoking

for (cancer in cancers) {
  cat("\n==========", cancer, "==========\n")
  f_m <- processed_path(cancer, "merged.rds")
  if (!file.exists(f_m)) stop(f_m, " 없음 → 02b_merge_genes.R 먼저 실행")
  d <- readRDS(f_m)
  qc <- list()
  add_qc <- function(item, value) qc[[length(qc) + 1]] <<- data.frame(item = item, value = as.character(value))
  add_qc("merged.rds patients", nrow(d))

  # ---- 1. signature 점수 ----
  for (s in names(pc$scores)) {
    sc <- pc$scores[[s]]
    sc_col <- paste0("score_", s); gr_col <- paste0("group_", s)
    if (sc$source == "sav") {
      miss <- setdiff(c(sc$score, sc$group), names(d))
      if (length(miss)) stop(cancer, ": 점수 열 없음 ", paste(miss, collapse = ", "), paper_where(paste0("scores$", s)))
      cols <- vapply(sc$check_genes, function(g) paper_gene_col(d, g, "_expression", cancer), "")
      # .sav 점수 정의: 결측 발현은 그 유전자 평균으로 채운 뒤 z (LUSC: OncoLnc CSV마다 1명씩 빠짐)
      ex <- as.matrix(d[, cols])
      imputed <- d$sample_id[!stats::complete.cases(ex)]
      ex <- apply(ex, 2, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))
      z <- rowMeans(scale(ex))
      diff <- abs(z - d[[sc$score]])
      bad <- d$sample_id[!is.na(diff) & diff > pc$z_check_tol]
      if (length(bad)) {
        stop(cancer, ": .sav 점수 ", sc$score, "가 ", paste(sc$check_genes, collapse = "/"),
             " 원 발현 z 평균과 다름 (", length(bad), "명) → ", paste(head(bad, 10), collapse = ", "),
             paper_where(paste0("scores$", s)))
      }
      add_qc(paste0(s, ": max |sav - recomputed z mean|"), signif(max(diff, na.rm = TRUE), 3))
      add_qc(paste0(s, ": gene expression missing -> gene mean imputed"), paste(imputed, collapse = " "))
      d[[sc_col]] <- d[[sc$score]]
      d[[gr_col]] <- as_group_factor(d[[sc$group]])
      ms <- check_median_split(d[[sc_col]], d[[gr_col]])
      if (!ms$median_split_pass) stop(cancer, ": ", sc$group, "이 점수의 median split과 ", ms$n_diff_from_median_split, "명 다름")
    } else if (sc$source == "sigscores") {
      sf <- processed_path(cancer, paste0("sigscores_", pc$signature, ".csv"))
      if (!file.exists(sf)) stop(sf, " 없음 → Rscript 11_sig_scores.R ", pc$signature, " ", cancer)
      x <- read.csv(sf, stringsAsFactors = FALSE)
      if (!sc$method %in% names(x)) stop(sf, ": 열 ", sc$method, " 없음", paper_where(paste0("scores$", s)))
      x <- data.frame(sample_id = x$sample_id, score = x[[sc$method]])
      x$group <- median_split(x$score)              # 점수 샘플 전체 (12와 같음)
      names(x) <- c("sample_id", sc_col, gr_col)
      d <- dplyr::left_join(d, x, by = "sample_id")
      add_qc(paste0(s, ": scored patients in cohort"), sum(!is.na(d[[sc_col]])))
    } else {
      stop("알 수 없는 점수 source: ", sc$source, paper_where(paste0("scores$", s, "$source")))
    }
    ms <- table(d[[gr_col]], useNA = "ifany")
    add_qc(paste0(s, ": group sizes"), paste(names(ms), ms, collapse = "; "))
    cat(sprintf("  점수 %s: %s\n", s, paste(names(ms), ms, collapse = " / ")))
  }

  # ---- 2. 흡연 상태 (GDC clinical) ----
  gf <- processed_path(cancer, "gdc_clinical.csv")
  if (!file.exists(gf) || refresh) {
    cat("  GDC clinical 다운로드 (TCGAbiolinks::GDCquery_clinic)\n")
    cl <- TCGAbiolinks::GDCquery_clinic(paste0(cfg$gdc$project_prefix, cancer), "clinical")
    keep <- c("submitter_id", sm$field, sm$pack_field)
    miss <- setdiff(keep, names(cl))
    if (length(miss)) stop(cancer, ": GDC clinical에 열 없음 ", paste(miss, collapse = ", "), paper_where("smoking"))
    readr::write_excel_csv(as.data.frame(cl[, keep]), gf, na = "")
  }
  g <- read.csv(gf, stringsAsFactors = FALSE, check.names = FALSE)
  g$sample_id <- substr(g$submitter_id, 1, 12)
  dup <- g$sample_id[duplicated(g$sample_id)]
  if (length(dup)) stop(gf, ": 중복 환자 ", paste(head(dup, 10), collapse = ", "))
  raw <- trimws(g[[sm$field]]); raw[raw == ""] <- NA
  unknown <- setdiff(unique(na.omit(raw)), c(names(sm$map), sm$na))
  if (length(unknown)) stop(cancer, ": 매핑 없는 흡연 값 ", paste(unknown, collapse = " | "), paper_where("smoking$map"))
  g$smoking <- factor(unname(sm$map[raw]), levels = sm$levels)
  g$smoking_raw <- raw
  g$gdc_pack_years <- suppressWarnings(as.numeric(g[[sm$pack_field]]))
  d <- dplyr::left_join(d, g[, c("sample_id", "smoking", "smoking_raw", "gdc_pack_years")], by = "sample_id")
  attr(d$smoking, "label") <- "Smoking status"

  # pack-years 대조 (.sav vs GDC)
  py <- d$number_pack_years_smoked
  dif <- abs(py - d$gdc_pack_years)
  bad <- setdiff(d$sample_id[!is.na(dif) & dif > sm$pack_tol], sm$pack_known_diff)
  if (length(bad)) {
    stop(cancer, ": .sav pack-years와 GDC pack_years_smoked 불일치 ", length(bad), "명 → ",
         paste(head(bad, 10), collapse = ", "), paper_where("smoking$pack_known_diff"))
  }
  add_qc("smoking: cohort patients not in GDC clinical", sum(!d$sample_id %in% g$sample_id))
  add_qc("pack-years: present in .sav only", sum(!is.na(py) & is.na(d$gdc_pack_years)))
  add_qc("pack-years: present in GDC only", sum(is.na(py) & !is.na(d$gdc_pack_years)))
  add_qc("pack-years: both present, equal", sum(!is.na(dif)))
  never_py <- d$sample_id[d$smoking %in% sm$levels[1] & !is.na(py) & py > 0]
  add_qc(paste0("pack-years > 0 in ", sm$levels[1], " smokers"), paste(never_py, collapse = " "))
  qs <- as.data.frame(table(raw = d$smoking_raw, mapped = d$smoking, useNA = "ifany"))
  qs <- qs[qs$Freq > 0, ]
  save_table(qs, cancer, "paper_QC_smoking.csv", pc$subdir)
  cat("  흡연:", paste(names(table(d$smoking, useNA = "ifany")), table(d$smoking, useNA = "ifany"), collapse = " / "), "\n")

  # ---- 3. tumor purity (ESTIMATE) ----
  ef <- file.path(cfg$clinical_dir, paste0("TCGA_", cancer, "_RNAseq_Expression.csv"))
  if (!file.exists(ef)) stop(ef, " 없음 → 02a_download_tcga.R 먼저 실행")
  m <- readr::read_csv(ef, show_col_types = FALSE, progress = FALSE)
  m <- as.data.frame(m)
  mm <- as.matrix(m[, -1])
  # 같은 기호의 행이 여러 개 → 평균 발현이 가장 큰 행
  ord <- order(m$Gene, -rowMeans(mm, na.rm = TRUE))
  keep <- ord[!duplicated(m$Gene[ord])]
  ex <- as.data.frame(mm[keep, , drop = FALSE]); rownames(ex) <- m$Gene[keep]
  add_qc("purity: expression matrix samples", ncol(ex))
  fg <- suppressMessages(tidyestimate::filter_common_genes(ex, id = "hgnc_symbol", tidy = FALSE,
                                                           tell_missing = FALSE, find_alias = TRUE))
  add_qc("purity: ESTIMATE common genes used", nrow(fg))
  es <- suppressMessages(tidyestimate::estimate_score(fg, is_affymetrix = FALSE))
  es <- data.frame(sample_id = substr(es$sample, 1, 12), stromal_score = es$stromal, immune_score = es$immune,
                   estimate_score = es$estimate, purity = cos(0.6049872018 + 0.0001467884 * es$estimate))
  if (anyDuplicated(es$sample_id)) stop(cancer, ": ESTIMATE 결과 중복 환자")
  write.csv(es, processed_path(cancer, "purity_estimate.csv"), row.names = FALSE)
  d <- dplyr::left_join(d, es, by = "sample_id")
  attr(d$purity, "label") <- "Tumor purity (ESTIMATE)"
  add_qc("purity: cohort patients with purity", sum(!is.na(d$purity)))
  cat(sprintf("  purity: %d명 (중앙값 %.3f, 범위 %.3f–%.3f)\n", sum(!is.na(d$purity)),
              median(d$purity, na.rm = TRUE), min(d$purity, na.rm = TRUE), max(d$purity, na.rm = TRUE)))

  saveRDS(d, processed_path(cancer, "paper.rds"))
  save_table(do.call(rbind, qc), cancer, "paper_QC.csv", pc$subdir)
  cat("  저장:", processed_path(cancer, "paper.rds"), "—", nrow(d), "명\n")
}
