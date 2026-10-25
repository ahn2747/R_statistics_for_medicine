# =============================================================
# R/paper.R
# 논문 분석 트랙 공통 함수 (20_paper_prep.R, 21_paper_analysis.R, 22_paper_geo.R)
# 사용: source("config.R"); source("R/utils.R"); source("R/survival.R"); source("R/paper.R")
# 설정은 cfg$paper (config.R)
# =============================================================

paper_where <- function(key) paste0(" → cfg$paper$", key)

# 논문 스크립트에서만 쓰는 연속형 단위 (pack-years per 10) → cfg$covariate_scale에 추가 (04 설정은 그대로)
paper_use_scales <- function() {
  cfg$covariate_scale <<- modifyList(cfg$covariate_scale, cfg$paper$covariate_scale %||% list())
  invisible(NULL)
}

# 유전자 기호 → 발현 열 (없으면 중단)
paper_gene_col <- function(d, gene, suffix = "_expression_log2", where = "") {
  col <- paste0(gene_key(gene), suffix)
  if (!col %in% names(d)) stop(where, ": 유전자 ", gene, " 열 ", col, " 없음 (merged.rds) → 02b 확인 / cfg$paper")
  col
}

# 종점 열 → <EP>_months / <EP>_event (0/1), 둘 다 있는 행만
paper_endpoint_data <- function(d, ep_name, where = "") {
  ep <- cfg$paper$endpoints[[ep_name]]
  if (is.null(ep)) stop(where, ": 종점 ", ep_name, " 정의 없음", paper_where("endpoints"))
  miss <- setdiff(c(ep$time, ep$event), names(d))
  if (length(miss)) stop(where, ": 종점 ", ep_name, " 열 없음: ", paste(miss, collapse = ", "), paper_where("endpoints"))
  tm <- as.numeric(d[[ep$time]])
  if (identical(ep$time_unit, "days")) tm <- tm / 30.44
  ev <- d[[ep$event]]
  bad <- unique(ev[!is.na(ev) & !ev %in% c(0, 1, ep$event_value)])
  if (length(bad)) stop(where, ": ", ep$event, " 값 ", paste(bad, collapse = ", "), paper_where(paste0("endpoints$", ep_name)))
  d$ep_months <- tm
  d$ep_event  <- as.numeric(ev == ep$event_value)
  neg <- d$sample_id[!is.na(tm) & tm < 0]
  if (length(neg)) stop(where, ": ", ep_name, " 시간 < 0 → ", paste(neg, collapse = ", "))
  d[!is.na(d$ep_months) & !is.na(d$ep_event), , drop = FALSE]
}

# ---- 상관 ---------------------------------------------------------------

# vars = 열 이름 벡터 (이름 = 라벨). 변수쌍마다 결측 제외 → long (var1, var2, r, p, n)
paper_cor <- function(d, vars, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  labs <- names(vars)
  g <- expand.grid(i = seq_along(vars), j = seq_along(vars))
  res <- t(mapply(function(i, j) {
    x <- d[[vars[i]]]; y <- d[[vars[j]]]
    ok <- !is.na(x) & !is.na(y)
    if (i == j) return(c(1, NA, sum(ok)))
    if (sum(ok) < 3) return(c(NA, NA, sum(ok)))
    ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = method, exact = FALSE))
    c(unname(ct$estimate), ct$p.value, sum(ok))
  }, g$i, g$j))
  data.frame(var1 = labs[g$i], var2 = labs[g$j], method = method,
             r = res[, 1], p = res[, 2], n = as.integer(res[, 3]))
}

# long → 표 (행 = var1, 열 = var2). 하삼각 "r (p)", 상삼각 n (pairwise)
paper_cor_wide <- function(long, labs) {
  out <- data.frame(Variable = labs, check.names = FALSE)
  for (j in seq_along(labs)) {
    out[[labs[j]]] <- vapply(seq_along(labs), function(i) {
      r <- long[long$var1 == labs[i] & long$var2 == labs[j], ]
      if (i == j) return("1")
      if (i > j) return(ifelse(is.na(r$r), "", sprintf("%.3f (%s)", r$r, fmt_p(r$p))))
      paste0("n = ", r$n)
    }, "")
  }
  out
}

paper_shapiro <- function(d, vars) {
  do.call(rbind, lapply(seq_along(vars), function(i) {
    x <- d[[vars[i]]]; x <- x[!is.na(x)]
    st <- if (length(x) >= 3 && length(x) <= 5000 && stats::sd(x) > 0) stats::shapiro.test(x) else NULL
    data.frame(variable = names(vars)[i], n = length(x),
               W = if (is.null(st)) NA else unname(st$statistic),
               p = if (is.null(st)) NA else st$p.value,
               normal = if (is.null(st)) NA else st$p.value >= 0.05)
  }))
}

# 부분 상관: x, y를 z(들)로 회귀한 잔차의 상관 (Spearman = 순위로 바꾼 뒤 같은 방식)
paper_partial_cor <- function(d, x, y, z, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  dd <- d[stats::complete.cases(d[, c(x, y, z)]), c(x, y, z), drop = FALSE]
  if (method == "spearman") dd[] <- lapply(dd, rank)
  rx <- stats::resid(stats::lm(dd[[x]] ~ ., data = dd[, z, drop = FALSE]))
  ry <- stats::resid(stats::lm(dd[[y]] ~ ., data = dd[, z, drop = FALSE]))
  n <- nrow(dd); k <- length(z)
  r <- stats::cor(rx, ry)
  tt <- r * sqrt((n - 2 - k) / (1 - r^2))
  data.frame(method = method, r = r, p = 2 * stats::pt(-abs(tt), df = n - 2 - k), n = n, df = n - 2 - k)
}

# ---- Table 1 (범주형, χ² 또는 Fisher) -----------------------------------

paper_table1 <- function(d, group, vars, labels) {
  g <- d[[group]]
  lv_g <- levels(g)
  rows <- list(); tests <- list()
  for (v in vars) {
    x <- d[[v]]
    if (!is.factor(x)) x <- factor(x)
    x <- droplevels(x)
    ok <- !is.na(x) & !is.na(g)
    tb <- table(x[ok], g[ok])
    ex <- suppressWarnings(stats::chisq.test(tb, correct = FALSE)$expected)
    low <- any(ex < 5)
    test <- if (nrow(tb) < 2) NULL else if (low) {
      list(name = "Fisher's exact", p = stats::fisher.test(tb, workspace = 2e8, simulate.p.value = nrow(tb) * ncol(tb) > 20,
                                                          B = 1e5)$p.value)
    } else {
      list(name = "Chi-square", p = stats::chisq.test(tb, correct = FALSE)$p.value)
    }
    miss <- tapply(is.na(x), g, sum)
    lab <- labels[[v]] %||% v
    for (l in rownames(tb)) {
      r <- data.frame(Variable = lab, Level = l, check.names = FALSE)
      for (k in lv_g) r[[k]] <- sprintf("%d (%.1f%%)", tb[l, k], 100 * tb[l, k] / sum(tb[, k]))
      r$Total <- sprintf("%d (%.1f%%)", sum(tb[l, ]), 100 * sum(tb[l, ]) / sum(tb))
      rows[[length(rows) + 1]] <- r
    }
    tests[[v]] <- data.frame(variable = v, label = lab, test = test$name %||% NA, p = test$p %||% NA,
                             n_tested = sum(tb), missing_total = sum(is.na(x) & !is.na(g)),
                             missing_by_group = paste(paste0(lv_g, " ", miss[lv_g]), collapse = "; "),
                             min_expected = min(ex), any_expected_lt5 = low,
                             n_cells_expected_lt5 = sum(ex < 5))
  }
  list(rows = do.call(rbind, rows), tests = do.call(rbind, tests))
}

# ---- 분위 그룹 ----------------------------------------------------------

# k 분위 (Q1 = 가장 낮음 = 기준). 분위 경계 = quantile type 7
paper_quantile_group <- function(x, k, prefix) {
  br <- stats::quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE, type = 7)
  if (anyDuplicated(br)) stop("분위 경계 중복 (값 동점이 많음): ", paste(round(br, 4), collapse = ", "))
  factor(cut(x, br, include.lowest = TRUE, labels = FALSE), levels = seq_len(k), labels = paste0(prefix, seq_len(k)))
}

# k 그룹 KM (팔레트: 2그룹이면 cfg$group_colors, 아니면 순서형 색)
paper_km_plot <- function(d, title, legend_title, ylab, p_lr) {
  fit <- survival::survfit(survival::Surv(ep_months, ep_event) ~ group, data = d)
  lv <- levels(d$group)
  pal <- if (identical(lv, cfg$group_levels)) unname(cfg$group_colors[lv]) else
    grDevices::colorRampPalette(unname(cfg$group_colors[cfg$group_levels]))(length(lv))
  p <- suppressMessages(survminer::ggsurvplot(
    fit, data = d, palette = pal, legend.labs = lv, legend.title = legend_title,
    pval = paste("Log-rank p", ifelse(p_lr < 0.001, "< 0.001", paste("=", sprintf("%.3f", p_lr)))),
    pval.size = 4, censor = TRUE, censor.shape = "|", censor.size = 3,
    risk.table = TRUE, risk.table.title = "Number at risk", risk.table.y.text = FALSE, fontsize = 3.5,
    break.time.by = 12, xlab = "Time (months)", ylab = ylab, title = title,
    ggtheme = ggplot2::theme_classic(base_size = 11), tables.theme = survminer::theme_cleantable()))
  suppressMessages(p$plot / p$table + patchwork::plot_layout(heights = c(3, 1)))
}

paper_logrank <- function(d) {
  sd <- survival::survdiff(survival::Surv(ep_months, ep_event) ~ group, data = d)
  stats::pchisq(sd$chisq, length(sd$n) - 1, lower.tail = FALSE)
}

# ---- GEO (22) -------------------------------------------------------------

# characteristics_ch1.* → 샘플 × key 표. R/geo.R의 parse_characteristics()와 같지만
# 한 샘플에 같은 key가 반복되면 열 순서대로 "key #2", "key #3"으로 구분 (예: GSE31210의 months before relapse/censor)
paper_parse_characteristics <- function(pd, gse) {
  cols <- grep("^characteristics_ch1", names(pd), value = TRUE)
  if (!length(cols)) stop(gse, ": characteristics_ch1 열 없음")
  long <- do.call(rbind, lapply(seq_along(cols), function(i) {
    x  <- as.character(pd[[cols[i]]])
    ok <- !is.na(x) & grepl(":", x, fixed = TRUE)
    data.frame(sample_id = rownames(pd)[ok], col = i,
               key   = trimws(sub(":.*$", "", x[ok])),
               value = trimws(sub("^[^:]*:", "", x[ok])), stringsAsFactors = FALSE)
  }))
  long <- long[order(long$sample_id, long$col), ]
  occ <- stats::ave(seq_len(nrow(long)), long$sample_id, long$key, FUN = seq_along)
  long$key <- ifelse(occ > 1, paste0(long$key, " #", occ), long$key)
  out <- data.frame(sample_id = rownames(pd), stringsAsFactors = FALSE)
  for (k in unique(long$key)) {
    v <- long[long$key == k, ]
    out[[k]] <- v$value[match(out$sample_id, v$sample_id)]
  }
  out
}

# config fields (표준 이름 → 원본 key, 여러 개면 앞에서부터 첫 비결측) → 표. 결측 표기 → NA, value_map 적용
paper_geo_map <- function(raw, g, gse) {
  na <- c(cfg$geo_na_values, g$na_values)
  d <- data.frame(sample_id = raw$sample_id, stringsAsFactors = FALSE)
  for (k in names(g$fields)) {
    keys <- g$fields[[k]]
    miss <- setdiff(keys, names(raw))
    if (length(miss)) stop(gse, ": 원본 key 없음 '", paste(miss, collapse = "', '"), "'", paper_where(paste0("geo$", gse, "$fields$", k)))
    x <- rep(NA_character_, nrow(raw))
    for (kk in keys) {
      v <- trimws(raw[[kk]])
      v[is.na(v) | v %in% na] <- NA
      x <- ifelse(is.na(x), v, x)
    }
    vm <- g$value_map[[k]]
    if (!is.null(vm)) {
      hit <- !is.na(x) & x %in% names(vm)
      x[hit] <- unname(vm[x[hit]])
    }
    d[[k]] <- x
  }
  d
}

# T/N/M 문자열 ("T2", "pN0pT1" 등) → AJCC 7판 병기 코드 "1"–"4" (모형에는 I–II vs III–IV만 사용)
#   M1 → 4; N2–3, T4, T3N1 → 3 (M 미상이어도 III–IV 그룹); M0이고 T1–2N0 → 1, T3N0·T1–2N1 → 2; 그 외 NA
paper_stage_tnm <- function(t, n, m = NULL, assume_m0 = FALSE) {
  dig <- function(x, L) suppressWarnings(as.integer(substr(stringr::str_extract(toupper(x), paste0(L, "[0-4X]")), 2, 2)))
  tt <- dig(t, "T"); nn <- dig(n, "N")
  mm <- if (is.null(m)) rep(if (assume_m0) 0L else NA_integer_, length(tt)) else dig(m, "M")
  st <- rep(NA_character_, length(tt))
  st[mm %in% 1] <- "4"
  iii <- is.na(st) & ((nn %in% 2:3) | (tt %in% 4) | (tt %in% 3 & nn %in% 1))
  st[iii] <- "3"
  m0 <- is.na(st) & mm %in% 0
  st[m0 & tt %in% 1:2 & nn %in% 0] <- "1"
  st[m0 & ((tt %in% 3 & nn %in% 0) | (tt %in% 1:2 & nn %in% 1))] <- "2"
  st
}
