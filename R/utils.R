# =============================================================
# R/utils.R
# 모든 스크립트가 공유하는 함수 (암종/유전자 무관)
# 사용: source("config.R"); source("R/utils.R")
# =============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# ---- 상수 --------------------------------------------------------------

# TCGA 결측 표기
na_strings <- c("", "[Not Available]", "[Unknown]", "[Not Applicable]",
                "[Discrepancy]", "[Not Evaluated]")

# 0/1 등 코드 → 라벨 (새 암종/변수는 여기에 추가)
value_labels <- list(
  pathologic_stage       = c(`1` = "I", `2` = "II", `3` = "III", `4` = "IV"),
  pathologic_stage_12_34 = c(`0` = "I–II", `1` = "III–IV"),
  cea_g                  = c(`0` = "≤5 ng/mL", `1` = ">5 ng/mL"),  # 5.0은 0
  age_g                  = c(Low = "≤65", High = "≥66")
)

# Table/그림용 변수 라벨
var_labels <- c(
  age                                 = "Age (years)",
  age_g                               = "Age group",
  gender                              = "Sex",
  stage                               = "Pathologic stage",
  pathologic_stage_12_34              = "Pathologic stage (I–II vs III–IV)",
  pathologic_t                        = "Pathologic T",
  pathologic_n                        = "Pathologic N",
  pathologic_m                        = "Pathologic M",
  cea_g                               = "Pretreatment CEA",
  histologic_diagnosis                = "Histology",
  anatomic_neoplasm_subdivision       = "Tumor site",
  residual_tumor                      = "Residual tumor",
  lymphovascular_invasion_indicator   = "Lymphovascular invasion",
  vascular_invasion_indicator         = "Vascular invasion",
  perineural_invasion                 = "Perineural invasion",
  # 검사 "시행 여부" 변수 — 변이/MMR 결과가 아님
  kras_gene_analysis_indicator        = "KRAS testing performed",
  braf_gene_analysis_indicator        = "BRAF testing performed",
  mismatch_rep_proteins_tested_by_ihc = "MMR IHC performed"
)

# ---- 파일/경로 ---------------------------------------------------------

out_dir <- function(cancer, type = c("tables", "figures")) {
  type <- match.arg(type)
  d <- file.path(cfg$output_dir, type, cancer)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

processed_path <- function(cancer, suffix) {
  dir.create(cfg$processed_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(cfg$processed_dir, paste0(cancer, "_", suffix))
}

save_table <- function(df, cancer, name) {
  path <- file.path(out_dir(cancer, "tables"), name)
  write.csv(df, path, row.names = FALSE, na = "")
  invisible(path)
}

# ---- 열 이름 -----------------------------------------------------------

# 유전자 기호 → 열 이름 키 (MS4A1 → ms4a1, HLA-A → hla_a)
gene_key <- function(gene) tolower(gsub("[^A-Za-z0-9]", "_", gene))

# 유전자 열(<GENE>expression / <GENE>Group 등, 대소문자 무관) → <gene>_expression / <gene>_group
# 나머지 열 → janitor::clean_names (ageG/AgeG → age_g 등)
# clean_names만 쓰면 APCexpression → ap_cexpression 이 되므로 유전자 열을 먼저 처리
standardize_names <- function(df) {
  nm     <- names(df)
  suffix <- str_match(nm, regex("(expression|group)$", ignore_case = TRUE))[, 2]
  stem   <- str_remove(nm, regex("(expression|group)$", ignore_case = TRUE))
  # 유전자 기호: 대문자/숫자로 시작, 밑줄 없음, 대문자 1개 이상
  is_gene <- !is.na(suffix) &
    str_detect(stem, "^[A-Z0-9][A-Za-z0-9.-]*$") & str_detect(stem, "[A-Z]")

  new <- nm
  new[is_gene]  <- paste0(gene_key(stem[is_gene]), "_", tolower(suffix[is_gene]))
  new[!is_gene] <- janitor::make_clean_names(nm[!is_gene])
  if (anyDuplicated(new)) {
    stop("열 이름 표준화 후 중복: ", paste(unique(new[duplicated(new)]), collapse = ", "))
  }
  names(df) <- new
  df
}

# <gene>_expression 과 <gene>_group 이 모두 있는 유전자 키
detect_genes <- function(df) {
  g <- str_remove(grep("_group$", names(df), value = TRUE), "_group$")
  g[paste0(g, "_expression") %in% names(df)]
}

# ---- 값 정리 -----------------------------------------------------------

# 발현 그룹 → factor(Low, High). 0/1 숫자는 1 = High
as_group_factor <- function(x) {
  if (is.numeric(x)) {
    if (!all(x %in% c(0, 1, NA))) stop("그룹 값이 0/1이 아님: ", paste(unique(x), collapse = ", "))
    x <- c("Low", "High")[x + 1]
  }
  x <- str_to_title(str_squish(as.character(x)))
  x[x %in% c("", "Na")] <- NA
  bad <- setdiff(unique(na.omit(x)), c("Low", "High"))
  if (length(bad)) stop("알 수 없는 그룹 값: ", paste(bad, collapse = ", "))
  factor(x, levels = c("Low", "High"))
}

# 공백 제거, TCGA 결측 표기 → NA, 숫자형 문자 → numeric,
# 중복 열 제거, 발현 그룹 → factor, 생존 상태 0/1 확인
clean_tcga <- function(df) {
  df <- df[, setdiff(names(df), cfg$drop_columns), drop = FALSE]

  df[] <- lapply(df, function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) {
      x <- str_squish(x)
      x[x %in% na_strings] <- NA
      num <- suppressWarnings(as.numeric(x))
      if (any(!is.na(x)) && all(is.na(x) | !is.na(num))) x <- num   # "2.00", ".00" → 2, 0
    }
    x
  })

  for (g in detect_genes(df)) {
    df[[paste0(g, "_group")]] <- as_group_factor(df[[paste0(g, "_group")]])
  }

  st <- cfg$surv_status
  if (st %in% names(df)) {
    df[[st]] <- as.numeric(df[[st]])
    if (!all(df[[st]] %in% c(0, 1, NA))) stop("생존 상태가 0/1이 아님")
  }
  df
}

# 코드 → 라벨 factor (매핑에 없는 값이 있으면 중단)
apply_value_labels <- function(df, labels = value_labels) {
  for (v in intersect(names(labels), names(df))) {
    map <- labels[[v]]
    x   <- as.character(df[[v]])
    bad <- setdiff(unique(na.omit(x)), names(map))
    if (length(bad)) stop(v, ": 라벨 매핑에 없는 값 ", paste(bad, collapse = ", "))
    df[[v]] <- factor(unname(map[x]), levels = unname(map))
  }
  df
}

# 숫자 코드 → 접두어 붙인 factor (T 2 → "T2")
code_factor <- function(x, prefix) {
  lv <- sort(unique(na.omit(x)))
  factor(ifelse(is.na(x), NA, paste0(prefix, x)), levels = paste0(prefix, lv))
}

set_var_labels <- function(df, labels = var_labels) {
  for (v in intersect(names(labels), names(df))) attr(df[[v]], "label") <- labels[[v]]
  df
}

# ---- 임상 데이터 -------------------------------------------------------

load_clinical <- function(cancer) {
  f <- file.path(cfg$clinical_dir, paste0(cancer, ".sav"))
  if (file.exists(f)) {
    df <- haven::read_sav(f) |>
      haven::as_factor(only_labelled = TRUE) |>
      haven::zap_formats() |>
      haven::zap_widths() |>
      as.data.frame()
    df <- standardize_names(df)
  } else {
    message(cancer, ": ", f, " 없음 → GDC 임상 데이터 사용")
    df <- load_gdc_clinical(cancer)
  }
  clean_tcga(df)
}

# .sav가 없는 암종: GDC 임상 → 핵심 열(sample_id, age, gender, stage, days, status)
load_gdc_clinical <- function(cancer) {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop(cancer, ".sav가 없고 TCGAbiolinks도 설치되지 않음 → Rscript 00_setup.R 실행 후 재시도")
  }
  cl <- TCGAbiolinks::GDCquery_clinic(paste0("TCGA-", cancer), type = "clinical")
  roman <- c(I = 1, II = 2, III = 3, IV = 4)
  stage_main <- str_match(cl$ajcc_pathologic_stage, "Stage (I{1,3}V?|IV)")[, 2]
  dead <- cl$vital_status == "Dead"
  data.frame(
    sampleID = cl$submitter_id,
    age_at_initial_pathologic_diagnosis = cl$age_at_index,
    gender = toupper(cl$gender),
    history_neoadjuvant_treatment = NA_character_,
    pathologic_stage = unname(roman[stage_main]),
    days = ifelse(dead, cl$days_to_death, cl$days_to_last_follow_up),
    status = ifelse(is.na(cl$vital_status), NA, as.numeric(dead)),
    stringsAsFactors = FALSE
  ) |> standardize_names()
}

# 분석용 변수 생성 (없는 열은 건너뜀 → 다른 암종에도 사용 가능)
recode_clinical <- function(df) {
  df$os_months <- df[[cfg$surv_time]] / 30.44

  if ("age_at_initial_pathologic_diagnosis" %in% names(df)) {
    df$age <- df$age_at_initial_pathologic_diagnosis
  }
  if ("gender" %in% names(df)) {
    df$gender <- factor(str_to_title(df$gender), levels = c("Male", "Female"))
  }
  tnm <- c(pathologic_t = "T", pathologic_n = "N", pathologic_m = "M")
  for (v in intersect(names(tnm), names(df))) df[[v]] <- code_factor(df[[v]], tnm[[v]])

  df <- apply_value_labels(df)
  if ("pathologic_stage" %in% names(df)) df$stage <- df$pathologic_stage

  # 검사 시행 여부 변수: YES/NO → Yes/No (결과 변수로 사용하지 않음)
  for (v in intersect(c("kras_gene_analysis_indicator", "braf_gene_analysis_indicator",
                        "mismatch_rep_proteins_tested_by_ihc"), names(df))) {
    df[[v]] <- factor(str_to_title(df[[v]]), levels = c("No", "Yes"))
  }

  neo <- df[["history_neoadjuvant_treatment"]]
  df$neoadjuvant_flag <- !is.na(neo) & tolower(neo) == "yes"

  set_var_labels(df)
}

# ---- 발현 그룹 ---------------------------------------------------------

# 중앙값 초과 = High (새 유전자에만 사용; 기존 .sav 그룹은 덮어쓰지 않음)
median_split <- function(x) {
  factor(ifelse(is.na(x), NA, ifelse(x > median(x, na.rm = TRUE), "High", "Low")),
         levels = c("Low", "High"))
}

# 기존 그룹이 median split과 일치하는지 QC (홀수 n의 중앙값 동점 1명 차이 허용)
check_median_split <- function(expr, group, tol = cfg$median_tie_tolerance) {
  ok <- !is.na(expr) & !is.na(group)
  n_diff <- sum(as.character(median_split(expr[ok])) != as.character(group[ok]))
  data.frame(n = sum(ok),
             n_low = sum(group[ok] == "Low"), n_high = sum(group[ok] == "High"),
             median = median(expr[ok]),
             n_diff_from_median_split = n_diff,
             median_split_pass = n_diff <= tol)
}

# <CANCER>_<split>_<GENE>.csv → sample_id, expression, group, csv_status(0/1)
read_gene_file <- function(path) {
  parts <- str_match(basename(path), cfg$gene_pattern)
  if (is.na(parts[1])) stop("파일명 규칙 불일치: ", basename(path))

  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("Patient", "Expression", "Group")
  if (!all(need %in% names(d))) {
    stop(basename(path), ": 필요한 열 없음 (", paste(setdiff(need, names(d)), collapse = ", "), ")")
  }

  status <- rep(NA_real_, nrow(d))
  if ("Status" %in% names(d)) {
    s <- str_to_title(str_squish(d$Status))
    bad <- setdiff(unique(na.omit(s)), c("Alive", "Dead", ""))
    if (length(bad)) stop(basename(path), ": 알 수 없는 Status 값 ", paste(bad, collapse = ", "))
    status <- c(Alive = 0, Dead = 1)[s]
  }

  out <- data.frame(
    sample_id  = str_sub(str_squish(d$Patient), 1, 12),
    expression = as.numeric(d$Expression),
    group      = as_group_factor(d$Group),
    csv_status = unname(status),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(out$sample_id)) {
    stop(basename(path), ": 중복 환자 ", paste(unique(out$sample_id[duplicated(out$sample_id)]), collapse = ", "))
  }
  list(cancer = toupper(parts[2]), split = parts[3], gene = parts[4],
       key = gene_key(parts[4]), file = basename(path), data = out)
}

# 02a 발현 행렬(log2(TPM+1))에서 유전자 추출 → TPM으로 되돌리고 median split
# (다른 발현 열과 같은 선형 척도 유지). ids: median 계산에 쓸 환자
extract_gdc_genes <- function(cancer, genes, ids) {
  f <- file.path(cfg$clinical_dir, paste0("TCGA_", cancer, "_RNAseq_Expression.csv"))
  if (!file.exists(f)) stop(f, " 없음 → 02a_download_tcga.R 먼저 실행")
  m <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE)
  miss <- setdiff(genes, m$Gene)
  if (length(miss)) warning(cancer, ": 발현 행렬에 없는 유전자 ", paste(miss, collapse = ", "))
  m <- m[m$Gene %in% genes, ]

  out <- data.frame(sample_id = ids)
  for (i in seq_len(nrow(m))) {
    key  <- gene_key(m$Gene[i])
    expr <- 2^unlist(m[i, -1]) - 1
    x    <- unname(expr[ids])
    out[[paste0(key, "_expression")]] <- x
    out[[paste0(key, "_group")]]      <- median_split(x)
  }
  out
}
