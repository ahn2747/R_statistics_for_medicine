# =============================================================
# 01_import_clinical.R
# <CANCER>.sav 읽기 → 열 이름 통일 → 정리/재코딩 → QC
# 결과:
#   data/processed/<CANCER>_clinical_all.rds  (전체 환자, neoadjuvant_flag 포함)
#   data/processed/<CANCER>_clinical.rds      (분석용: cfg$exclude_neoadjuvant 적용)
#   output/tables/<CANCER>/clinical_QC_median_split.csv
# =============================================================

source("config.R")
source("R/utils.R")

surv_counts <- function(d) {
  ok <- !is.na(d[[cfg$surv_time]]) & !is.na(d[[cfg$surv_status]])
  c(n = sum(ok), ev = sum(d[[cfg$surv_status]][ok] == 1))
}

for (cancer in cfg$cancers) {
  cat("\n==========", cancer, "==========\n")
  df <- load_clinical(cancer) |> recode_clinical()

  genes <- detect_genes(df)
  cat("환자", nrow(df), "명 /", ncol(df), "열\n")
  cat("열:", paste(names(df), collapse = ", "), "\n")
  cat("유전자:", paste(toupper(genes), collapse = ", "), "\n")

  # ---- 생존: Days NA 행 == Status NA 행 ----
  na_time   <- is.na(df[[cfg$surv_time]])
  na_status <- is.na(df[[cfg$surv_status]])
  if (!identical(na_time, na_status)) {
    stop(cancer, ": days/status 결측 위치 불일치 → ",
         paste(df$sample_id[na_time != na_status], collapse = ", "))
  }

  # ---- 신보조요법 환자 ----
  neo_ids <- df$sample_id[df$neoadjuvant_flag]
  cat("신보조요법(neoadjuvant) 환자", length(neo_ids), "명:",
      if (length(neo_ids)) paste(neo_ids, collapse = ", ") else "-", "\n")

  # ---- 생존 n / 사건 수 검증 ----
  counts <- list(incl = surv_counts(df), excl = surv_counts(df[!df$neoadjuvant_flag, ]))
  cat(sprintf("생존분석 n (사건): 제외 전 %d (%d), 신보조요법 제외 후 %d (%d)\n",
              counts$incl["n"], counts$incl["ev"], counts$excl["n"], counts$excl["ev"]))
  expected <- cfg$expected_surv[[cancer]]
  if (!is.null(expected)) {
    for (k in names(counts)) {
      if (any(counts[[k]] != expected[[k]][c("n", "ev")])) {
        stop(cancer, " (", k, "): 생존 n/사건 = ", paste(counts[[k]], collapse = "/"),
             ", 기대값 = ", paste(expected[[k]], collapse = "/"))
      }
    }
    cat("  → cfg$expected_surv 와 일치\n")
  }

  # ---- 기존 유전자 그룹 QC (재계산 값으로 덮어쓰지 않음, 전체 코호트 기준) ----
  qc <- bind_rows(lapply(genes, function(g) {
    data.frame(gene = toupper(g),
               check_median_split(df[[paste0(g, "_expression")]], df[[paste0(g, "_group")]]))
  }))
  print(qc, row.names = FALSE)
  if (any(!qc$median_split_pass)) {
    warning(cancer, ": median split과 ", cfg$median_tie_tolerance, "명 넘게 다른 유전자 → ",
            paste(qc$gene[!qc$median_split_pass], collapse = ", "))
  }
  save_table(qc, cancer, "clinical_QC_median_split.csv")

  # ---- 저장 ----
  saveRDS(df, processed_path(cancer, "clinical_all.rds"))
  analytic <- if (cfg$exclude_neoadjuvant) df[!df$neoadjuvant_flag, ] else df
  cat("분석용 환자:", nrow(df), "→", nrow(analytic),
      if (cfg$exclude_neoadjuvant) "(신보조요법 제외)" else "(제외 없음)", "\n")
  saveRDS(analytic, processed_path(cancer, "clinical.rds"))
}
