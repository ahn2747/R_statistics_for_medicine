# =============================================================
# R/geo.R
# GEO 데이터셋 공통 함수 (06_external_geo.R, 13_sig_geo.R)
#   다운로드/캐시, log2 확인, probe → 유전자, characteristics_ch1 파싱, 필드 매핑, 재코딩, pheno QC, 종점 데이터
# 사용: source("config.R"); source("R/utils.R"); source("R/geo.R")
# 설정: cfg$geo_dir, cfg$geo_probe_rule, cfg$geo_na_values, cfg$geo_datasets (config.R)
# =============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

options(timeout = max(1200, getOption("timeout")))    # 큰 series matrix 다운로드
Sys.setenv(VROOM_CONNECTION_SIZE = 131072 * 64)      # 긴 series matrix 줄 읽기 오류 방지

geo_where <- function(gse, key) paste0(" → cfg$geo_datasets$", gse, "$", key)
ids_txt <- function(x) paste(head(x, 10), collapse = ", ")

# ---- 1. 다운로드 / 캐시 ------------------------------------------------

load_geo_eset <- function(gse, g) {
  dir.create(cfg$geo_dir, recursive = TRUE, showWarnings = FALSE)
  rds <- file.path(cfg$geo_dir, paste0(gse, "_", g$platform, "_eset.rds"))
  if (file.exists(rds)) {
    cat("[캐시] ", rds, "\n", sep = "")
    eset <- readRDS(rds)
  } else {
    cat("[다운로드] ", gse, " → ", cfg$geo_dir, " (받아 둔 series matrix/GPL 파일이 있으면 재사용)\n", sep = "")
    lst <- getGEO(gse, GSEMatrix = TRUE, destdir = cfg$geo_dir, getGPL = TRUE)
    idx <- if (length(lst) == 1) 1 else grep(g$platform, names(lst), fixed = TRUE)
    if (length(idx) != 1) {
      stop(gse, ": platform '", g$platform, "'의 series matrix를 하나로 특정할 수 없음 (",
           paste(names(lst), collapse = ", "), ")", geo_where(gse, "platform"))
    }
    eset <- lst[[idx]]
    saveRDS(eset, rds)
  }
  if (annotation(eset) != g$platform) {
    stop(gse, ": platform이 ", annotation(eset), " (설정 ", g$platform, ")", geo_where(gse, "platform"))
  }
  eset
}

# log2 척도 확인 (max > 100 → log2 변환)
prepare_expression <- function(eset, gse) {
  ex <- exprs(eset)
  rng <- range(ex, na.rm = TRUE)
  cat(sprintf("발현 행렬: %d probes × %d samples, 범위 %.2f–%.2f, NA %d\n",
              nrow(ex), ncol(ex), rng[1], rng[2], sum(is.na(ex))))
  if (rng[2] > 100) {
    if (rng[1] < 0) stop(gse, ": 선형 척도로 보이는데 음수 값이 있어 log2 변환 불가")
    off <- if (rng[1] < 1) 1 else 0
    ex <- log2(ex + off)
    cat(sprintf("  max > 100 → log2%s 변환, 새 범위 %.2f–%.2f\n",
                if (off) "(x + 1)" else "(x)", min(ex, na.rm = TRUE), max(ex, na.rm = TRUE)))
  } else {
    cat("  이미 log2 척도 (max ≤ 100) → 변환 없음\n")
  }
  ex
}

# ---- 2. probe → 유전자 -------------------------------------------------

# GPL 주석의 유전자 기호 열(g$symbol_col) → probe별 기호 목록 ("A /// B" 다중 매핑 분리)
geo_symbols <- function(fd, g, gse) {
  sc <- g$symbol_col
  if (!sc %in% names(fd)) {
    stop(gse, ": GPL 주석에 '", sc, "' 열 없음 (열: ", paste(names(fd), collapse = ", "), ")",
         geo_where(gse, "symbol_col"))
  }
  lapply(strsplit(as.character(fd[[sc]]), "///", fixed = TRUE), trimws)
}

# cfg$geo_probe_rule → probe 선택에 쓸 통계량 열 이름 (mean / IQR / sd; 가장 큰 probe 선택)
probe_rule_col <- function() {
  rules <- c(max_mean = "mean", max_iqr = "IQR", max_sd = "sd")
  rule  <- cfg$geo_probe_rule
  if (!rule %in% names(rules)) {
    stop("cfg$geo_probe_rule = '", rule, "' 지원 안 함 (", paste(names(rules), collapse = ", "), ")")
  }
  rules[[rule]]
}

# 유전자에 매칭되는 모든 probe의 QC (tumor 샘플 기준) + cfg$geo_probe_rule로 선택
probe_qc <- function(ex, fd, gene, g, gse) {
  sc <- g$symbol_col
  syms <- geo_symbols(fd, g, gse)
  if (!identical(rownames(ex), rownames(fd))) stop(gse, ": 발현 행렬과 GPL 주석의 probe 순서 불일치")
  hit  <- vapply(syms, function(s) gene %in% s, logical(1))
  probes <- rownames(fd)[hit]
  if (!length(probes)) stop(gse, ": '", sc, "'에서 ", gene, " probe를 찾지 못함")

  sub <- ex[probes, , drop = FALSE]
  qc <- data.frame(
    probe_id     = probes,
    symbol       = as.character(fd[probes, sc]),
    multi_mapped = lengths(syms[hit]) > 1,
    n            = rowSums(!is.na(sub)),
    mean         = rowMeans(sub, na.rm = TRUE),
    median       = apply(sub, 1, median, na.rm = TRUE),
    IQR          = apply(sub, 1, IQR, na.rm = TRUE),
    sd           = apply(sub, 1, sd, na.rm = TRUE),
    row.names = NULL, stringsAsFactors = FALSE
  )
  rule  <- cfg$geo_probe_rule
  qc$chosen <- seq_len(nrow(qc)) == which.max(qc[[probe_rule_col()]])
  if (length(probes) > 1) {
    r <- cor(t(sub), use = "pairwise.complete.obs")
    for (p in probes) qc[[paste0("r_", p)]] <- r[probes, p]
  }
  qc$rule <- rule
  qc
}

# ---- 3. phenotype ------------------------------------------------------

# characteristics_ch1.* ("key: value") → 샘플 × key 표 (원본 key 이름 그대로)
parse_characteristics <- function(pd, gse) {
  cols <- grep("^characteristics_ch1", names(pd), value = TRUE)
  if (!length(cols)) stop(gse, ": characteristics_ch1 열 없음")
  long <- do.call(rbind, lapply(cols, function(cl) {
    x  <- as.character(pd[[cl]])
    ok <- !is.na(x) & grepl(":", x, fixed = TRUE)
    data.frame(sample_id = rownames(pd)[ok],
               key   = trimws(sub(":.*$", "", x[ok])),
               value = trimws(sub("^[^:]*:", "", x[ok])),
               stringsAsFactors = FALSE)
  }))
  dup <- duplicated(long[, c("sample_id", "key")])
  if (any(dup)) {
    stop(gse, ": 한 샘플에 같은 key가 여러 번 → ", ids_txt(unique(paste(long$sample_id[dup], long$key[dup]))))
  }
  out <- data.frame(sample_id = rownames(pd), stringsAsFactors = FALSE)
  for (k in unique(long$key)) {
    v <- long[long$key == k, ]
    out[[k]] <- v$value[match(out$sample_id, v$sample_id)]
  }
  out
}

is_geo_missing <- function(x) is.na(x) | trimws(x) %in% cfg$geo_na_values

print_raw_keys <- function(raw) {
  keys <- setdiff(names(raw), "sample_id")
  tab <- data.frame(
    key = keys,
    n_present = vapply(keys, function(k) sum(!is_geo_missing(raw[[k]])), 0L),
    n_unique  = vapply(keys, function(k) length(unique(raw[[k]][!is_geo_missing(raw[[k]])])), 0L),
    examples  = vapply(keys, function(k) {
      u <- unique(raw[[k]][!is_geo_missing(raw[[k]])])
      substr(paste(head(u, 5), collapse = " | "), 1, 60)
    }, ""),
    row.names = NULL
  )
  cat("\n[characteristics_ch1 원본 key] (", nrow(raw), " samples; config fields에 이 key를 지정)\n", sep = "")
  print(tab, right = FALSE, row.names = FALSE)
  invisible(tab)
}

# config fields (표준 이름 → 원본 key) 적용. 결측 표기 → NA
map_fields <- function(raw, g, gse) {
  f <- g$fields
  ep_fields <- unlist(lapply(g$endpoints, function(e) c(e$time, e$event)))
  req <- unique(c("sample_type", "age", "gender", "stage", ep_fields,
                  g$extra_covariate, g$subgroup$var))
  miss <- setdiff(req, names(f))
  if (length(miss)) stop(gse, ": fields에 필수 항목 없음 (", paste(miss, collapse = ", "), ")", geo_where(gse, "fields"))
  for (k in names(f)) {
    if (!f[[k]] %in% names(raw)) {
      near <- agrep(f[[k]], names(raw), max.distance = 0.3, value = TRUE, ignore.case = TRUE)
      stop(gse, ": 원본 key '", f[[k]], "' 없음", geo_where(gse, paste0("fields[\"", k, "\"]")),
           if (length(near)) paste0(". 비슷한 key: ", paste(head(near, 5), collapse = ", ")) else "")
    }
  }
  d <- data.frame(sample_id = raw$sample_id, stringsAsFactors = FALSE)
  for (k in names(f)) {
    x <- trimws(raw[[f[[k]]]])
    x[is_geo_missing(x)] <- NA
    d[[k]] <- x
  }
  d
}

num_field <- function(d, k, g, gse) {
  x <- d[[k]]
  v <- suppressWarnings(as.numeric(x))
  bad <- !is.na(x) & is.na(v)
  if (any(bad)) {
    stop(gse, ": ", k, " (원본 key '", g$fields[[k]], "')에 숫자가 아닌 값 ",
         paste(unique(x[bad]), collapse = ", "), " → ", ids_txt(d$sample_id[bad]),
         " (결측 표기면 cfg$geo_na_values에 추가)")
  }
  v
}

# 비종양 제외 + 공변량 / 종점 변수 생성
recode_geo <- function(d, g, gse) {
  if (any(is.na(d$sample_type))) {
    stop(gse, ": sample_type 결측 → ", ids_txt(d$sample_id[is.na(d$sample_type)]),
         geo_where(gse, "fields[\"sample_type\"]"))
  }
  cat("\n샘플 유형:\n")
  print(table(d$sample_type))
  excl <- d$sample_type %in% g$exclude_sample_type
  cat(sprintf("비종양 제외 (%s): %d → 종양 %d\n",
              paste(g$exclude_sample_type, collapse = ", "), sum(excl), sum(!excl)))
  d <- d[!excl, ]
  if (!nrow(d)) stop(gse, ": 종양 샘플 없음", geo_where(gse, "exclude_sample_type"))

  # 나이
  d$age <- num_field(d, "age", g, gse)
  if (any(d$age < 0 | d$age > 120, na.rm = TRUE)) stop(gse, ": 나이 값 이상 (0–120 아님)")

  # 성별 (기준 수준: cfg$reference_levels)
  d$gender <- factor(str_to_title(d$gender))
  if (nlevels(d$gender) != 2) warning(gse, ": 성별 값이 2개가 아님 (", paste(levels(d$gender), collapse = ", "), ")", call. = FALSE)

  # stage: stage_na_values → NA (명시), 나머지는 1–4 코드 (TCGA와 같은 라벨)
  s <- d$stage
  to_na <- !is.na(s) & s %in% g$stage_na_values
  if (any(to_na)) {
    cat(sprintf("stage %s → NA %d명 (cfg$geo_datasets$%s$stage_na_values): %s\n",
                paste(g$stage_na_values, collapse = "/"), sum(to_na), gse, ids_txt(d$sample_id[to_na])))
    s[to_na] <- NA
  }
  stage_labels <- labels_for(g$validates)
  bad <- !is.na(s) & !s %in% names(stage_labels$pathologic_stage)
  if (any(bad)) {
    stop(gse, ": stage 값 ", paste(unique(s[bad]), collapse = ", "), " 해석 불가 → ", ids_txt(d$sample_id[bad]),
         geo_where(gse, "stage_na_values"))
  }
  d$stage_num <- as.numeric(s)
  d$pathologic_stage <- d$stage_num
  d$pathologic_stage_12_34 <- ifelse(is.na(d$stage_num), NA, as.numeric(d$stage_num >= 3))
  d <- apply_value_labels(d, stage_labels[intersect(c("pathologic_stage", "pathologic_stage_12_34"),
                                                    names(stage_labels))])
  d$stage <- d$pathologic_stage

  # TCGA .sav에만 있는 그룹 변수 (예: age_g): validates 코호트의 label_rules(source > cutoff → high)와
  # value_labels로 생성 → cox_covariates를 TCGA와 같은 정의로 사용
  rules <- rules_for(g$validates)
  for (v in names(rules)) {
    r <- rules[[v]]
    if (v %in% names(d) || !r$source %in% names(d)) next
    codes <- names(stage_labels[[v]])
    if (length(codes) != 2 || !r$high %in% codes) {
      stop(gse, ": ", v, " 생성 불가 — cfg$value_labels_default$", v, "에 high 코드 ", r$high, " 포함 2개 수준 필요")
    }
    x <- d[[r$source]]
    d[[v]] <- ifelse(is.na(x), NA, ifelse(x > r$cutoff, r$high, setdiff(codes, r$high)))
    d <- apply_value_labels(d, stage_labels[v])
    cat(sprintf("%s 생성 (%s > %g → %s, cfg$label_rules_default$%s): %s\n", v, r$source, r$cutoff,
                stage_labels[[v]][[r$high]], v, paste(names(table(d[[v]])), table(d[[v]]), sep = " ", collapse = " / ")))
  }

  # 나머지 범주형 필드 → factor (factor_levels가 있으면 그 순서, 첫 번째 = 기준)
  ep_fields <- unlist(lapply(g$endpoints, function(e) c(e$time, e$event)))
  other <- setdiff(names(g$fields), c("sample_type", "age", "gender", "stage", ep_fields))
  for (k in other) {
    lv <- g$factor_levels[[k]]
    if (!is.null(lv)) {
      bad <- !is.na(d[[k]]) & !d[[k]] %in% lv
      if (any(bad)) {
        stop(gse, ": ", k, " 값 ", paste(unique(d[[k]][bad]), collapse = ", "), "이(가) factor_levels에 없음 → ",
             ids_txt(d$sample_id[bad]), geo_where(gse, paste0("factor_levels$", k)))
      }
      d[[k]] <- factor(d[[k]], levels = lv)
    } else {
      d[[k]] <- factor(d[[k]])
    }
  }

  # 종점: <EP>_months, <EP>_event (0/1)
  for (ep in names(g$endpoints)) {
    e <- g$endpoints[[ep]]
    tm <- num_field(d, e$time, g, gse)
    if (!e$time_unit %in% c("days", "months")) stop(gse, " ", ep, ": time_unit은 days/months", geo_where(gse, paste0("endpoints$", ep)))
    if (e$time_unit == "days") tm <- tm / 30.44
    if (any(tm < 0, na.rm = TRUE)) stop(gse, " ", ep, ": 음수 시간 → ", ids_txt(d$sample_id[which(tm < 0)]))
    ev_raw <- d[[e$event]]
    vals <- unique(na.omit(ev_raw))
    if (length(vals) > 2 || !as.character(e$event_value) %in% vals) {
      stop(gse, " ", ep, ": event 값 (", paste(vals, collapse = ", "), ")이 event_value ",
           e$event_value, "와 맞지 않음", geo_where(gse, paste0("endpoints$", ep, "$event_value")))
    }
    d[[paste0(ep, "_months")]] <- tm
    d[[paste0(ep, "_event")]]  <- code_status(ev_raw, e$event_value)
  }

  d <- apply_reference_levels(d)
  set_var_labels(d)
}

# 필드별 결측 / 값 요약 (종양 샘플)
pheno_qc <- function(d, g) {
  rows <- lapply(names(g$fields), function(k) {
    x <- d[[k]]
    if (k %in% c(unlist(lapply(g$endpoints, `[[`, "time")))) {
      ep <- names(g$endpoints)[vapply(g$endpoints, function(e) e$time == k, TRUE)]
      x  <- d[[paste0(ep, "_months")]]
    }
    if (k %in% c(unlist(lapply(g$endpoints, `[[`, "event")))) {
      ep <- names(g$endpoints)[vapply(g$endpoints, function(e) e$event == k, TRUE)]
      x  <- d[[paste0(ep, "_event")]]
    }
    vals <- if (is.numeric(x) && length(unique(na.omit(x))) > 5) {
      sprintf("median %.1f (range %.1f–%.1f)", median(x, na.rm = TRUE), min(x, na.rm = TRUE), max(x, na.rm = TRUE))
    } else {
      t <- table(x)
      paste(paste0(names(t), " = ", as.integer(t)), collapse = "; ")
    }
    data.frame(field = k, raw_key = g$fields[[k]], n = sum(!is.na(x)), n_missing = sum(is.na(x)),
               values = vals, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---- 데이터셋 불러오기 (06, 13 공통) --------------------------------------
# series matrix + log2 확인 + characteristics_ch1 → fields 매핑 → 비종양 제외/재코딩
# 반환: list(eset, ex (log2 발현, 전체 샘플), d (종양 샘플 phenotype), raw_keys)
geo_load_cohort <- function(gse, g) {
  eset <- load_geo_eset(gse, g)
  ex   <- prepare_expression(eset, gse)
  pd   <- pData(eset)
  cat(sprintf("phenotype: %d samples × %d columns\n", nrow(pd), ncol(pd)))

  # ---- phenotype ----
  raw <- parse_characteristics(pd, gse)
  raw_keys <- print_raw_keys(raw)
  d <- recode_geo(map_fields(raw, g, gse), g, gse)
  list(eset = eset, ex = ex, d = d, raw_keys = raw_keys)
}

# ---- 종점별 분석 데이터 ----------------------------------------------------
# 시간/사건 둘 다 있고 keep인 샘플, 종점의 exclude_stage 제외 (제외 수를 출력)
# keep: 06 = 발현 그룹 있음, 13 = 점수 있음
geo_endpoint_data <- function(d, ep, g, gse, tag, keep) {
  e     <- g$endpoints[[ep]]
  time  <- paste0(ep, "_months")
  event <- paste0(ep, "_event")
  has_ep <- !is.na(d[[time]]) & !is.na(d[[event]])
  one_side <- is.na(d[[time]]) != is.na(d[[event]])
  if (any(one_side)) {
    cat(sprintf("  [%s] 시간/사건 중 하나만 결측 %d명 → 이 종점에서 제외: %s\n",
                tag, sum(one_side), ids_txt(d$sample_id[one_side])))
  }
  de <- d[has_ep & keep, ]
  if (!is.null(e$exclude_stage)) {
    ex_st <- de$stage_num %in% e$exclude_stage
    cat(sprintf("  [%s] stage %s 제외 %d명 (cfg$geo_datasets$%s$endpoints$%s$exclude_stage)\n",
                tag, paste(e$exclude_stage, collapse = "/"), sum(ex_st), gse, ep))
    de <- de[!ex_st, ]
  }
  de
}
