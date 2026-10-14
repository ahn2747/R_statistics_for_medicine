# =============================================================
# 06_external_geo.R
# 주 가설 유전자(cfg$primary_gene)의 GEO 외부 검증 (데이터셋: cfg$geo_datasets)
#   1. GEOquery series matrix + GPL 주석 (cfg$geo_dir 캐시), log2 척도 확인
#   2. probe → 유전자: 매칭 probe 전부 QC, cfg$geo_probe_rule로 하나 선택
#   3. characteristics_ch1 파싱 → config fields 매핑, 비종양 샘플 제외
#   4. 종양 샘플 내 median split (cfg$group_levels, 첫 번째 = 기준)
#   5. 종점별 (첫 번째 = 주 종점): KM + log-rank, 단변량 Cox (유전자 + 임상 변수별), 다변량 Cox (EPV),
#      다변량 + extra_covariate (MMR), 하위군 (pMMR), cox.zph
#   6. TCGA 결과와 방향 비교 (validates 코호트의 survival_summary_raw.csv)
#   1–3단계와 종점 데이터 구성은 R/geo.R (13_sig_geo.R과 공유)
# 사용: Rscript 06_external_geo.R [GSE39582 ...]   (인자 없으면 cfg$geo_datasets 전체)
# 결과 (파일명 앞에 <GSE>_): output/tables/GEO/<GSE>/ (probe_QC, pheno_QC, survival_summary, km_summary, cox_<EP>, ph_tests)
#       output/figures/GEO/<GSE>/ (km_<EP>, km_<EP>_<하위군>, forest, forest_multi_<EP>)
# =============================================================

source("config.R")
source("R/utils.R")
source("R/survival.R")
source("R/geo.R")       # GEO 다운로드/pheno/probe/종점 공통 함수 (13_sig_geo.R과 공유)

# ---- 5. 생존분석 -------------------------------------------------------

model_label_extra <- function(g) var_labels[[g$extra_covariate]] %||% g$extra_covariate

# 한 종점 × 한 집단 (All 또는 하위군)
analyze_endpoint <- function(d, ep, g, gse, gene, pop = "All") {
  e     <- g$endpoints[[ep]]
  time  <- paste0(ep, "_months")
  event <- paste0(ep, "_event")
  out   <- list(summary = NULL, km = NULL, cox = NULL, zph = NULL, files = character(), forest = NULL)
  tag   <- if (pop == "All") ep else paste0(ep, "_", pop)

  de <- add_scaled_terms(geo_endpoint_data(d, ep, g, gse, tag, keep = !is.na(d$group)))
  cnt <- table(factor(de$group, levels = cfg$group_levels))
  n_ev <- sum(de[[event]])
  cat(sprintf("  [%s] n = %d (%s %d / %s %d), 사건 %d\n", tag, nrow(de),
              group_ref(), cnt[[group_ref()]], group_alt(), cnt[[group_alt()]], n_ev))
  # 분석 정보용 (실제 적합한 모형식은 아래 add_model 등에서 out$models에 추가)
  out$n_txt <- sprintf("%s/%s: n=%d (%s %d / %s %d); events=%d", ep, pop, nrow(de), group_ref(),
                       cnt[[group_ref()]], group_alt(), cnt[[group_alt()]], n_ev)
  out$models <- character()
  if (n_ev < cfg$min_events || any(cnt == 0)) {
    warning(gse, " ", tag, ": 사건 ", n_ev, " 또는 빈 그룹 → 건너뜀", call. = FALSE)
    return(out)
  }

  # ---- KM ----
  f  <- as.formula(paste0("Surv(", time, ", ", event, ") ~ group"))
  km <- surv_fit(f, data = de)
  p_lr <- survdiff(f, data = de)$pvalue
  out$km <- data.frame(endpoint = ep, population = pop, km_stats(km, cfg$km_times), logrank_p = p_lr)
  title <- paste0(gene, " expression — ", g$label, if (pop != "All") paste0(", ", pop))
  out$files <- c(out$files, suppressMessages(save_fig(
    km_plot(km, de, gene, NULL, p_lr, title = title, ylab = paste(e$label, "probability")),
    "GEO", paste0("km_", tag), width = 7, height = 6, subdir = gse, prefix = gse)))

  # ---- Cox 모형들 ----
  gt <- gene_term()   # 다변량/시간 분할 모형의 유전자 항 (cfg$cox_gene_term)
  # TCGA cox_covariates 중 이 데이터셋에 없는 변수는 cfg$geo_datasets$<GSE>$cox_covariates_drop에 명시
  cov_terms <- covariate_terms(setdiff(cfg$cox_covariates, g$cox_covariates_drop), de,
                               paste0("geo_datasets$", gse, "$cox_covariates_drop / cox_covariates"), paste(gse, tag))
  # 임상 변수 단변량 (04와 같은 cfg$cox_uni_covariates; 없는 변수는 cox_covariates_drop, 하위군에서는 extra 제외)
  uni_terms <- covariate_terms(setdiff(cfg$cox_uni_covariates, g$cox_covariates_drop), de,
                               paste0("geo_datasets$", gse, "$cox_covariates_drop / cox_uni_covariates"), paste(gse, tag))
  if (!is.null(g$extra_covariate) && pop == "All" && g$extra_covariate %in% names(de)) {
    uni_terms <- c(uni_terms, g$extra_covariate)
  }
  labs <- term_labels(unique(c("group", "expr_log2", cov_terms, uni_terms, g$extra_covariate, cfg$stage_collapsed)), gene)
  rows <- list(); cox_tabs <- list(); zph <- list()

  add_model <- function(name, m, key_term = gt, keep_table = FALSE) {
    fit <- fit_cox(m$d, m$terms, time = time, event = event)
    out$models <<- c(out$models, sprintf("%s: %s (n=%d, events=%d)", name,
                                         paste(deparse(formula(fit), width.cutoff = 500L), collapse = " "),
                                         nrow(m$d), m$events))
    tab <- cox_rows(fit, m$d, m$terms, labs[m$terms], event = event)
    z   <- zph_rows(fit, gene, name)
    k   <- tab[tab$term == key_term & !tab$reference, ]
    rows[[length(rows) + 1]] <<- data.frame(
      endpoint = ep, population = pop, model = name, term = key_term,
      n = nrow(m$d), events = m$events, dropped = m$dropped, params = m$k, epv = m$epv,
      stage_levels = m$stage_levels %||% NA, exploratory = m$exploratory %||% FALSE,
      hr = k$hr, lo = k$lo, hi = k$hi, p = k$p,
      ph_p_term = z$p[z$term == key_term], ph_p_global = z$p[z$term == "GLOBAL"],
      stringsAsFactors = FALSE)
    cox_tabs[[name]] <<- data.frame(Model = name, format_cox(tab), check.names = FALSE)
    zph[[name]] <<- data.frame(endpoint = ep, population = pop, z)
    if (keep_table) out$forest <<- list(tab = tab, name = name, m = m)
    fit
  }
  simple <- function(terms) {
    dm <- complete_cases(de, terms)
    list(terms = terms, d = dm, events = sum(dm[[event]]), k = n_params(dm, terms),
         epv = sum(dm[[event]]) / n_params(dm, terms), dropped = nrow(de) - nrow(dm),
         stage_levels = NA, exploratory = FALSE)
  }

  add_model("Univariable", simple("group"), key_term = "group")
  add_model("Univariable, continuous (per 1 log2 unit)", simple("expr_log2"), key_term = "expr_log2")
  # 임상 변수 단변량: Cox 표에만 (요약/forest/key_results는 유전자 항 기준, 04의 cox_uni와 같은 방식)
  uni_cov <- bind_rows(lapply(uni_terms, function(t) {
    dt <- complete_cases(de, t)
    cat(sprintf("  [%s] 단변량 %-24s n = %d, 사건 %d\n", tag, t, nrow(dt), sum(dt[[event]])))
    cox_rows(fit_cox(dt, t, time = time, event = event), dt, t, labs[t], event = event)
  }))
  if (nrow(uni_cov)) {
    cox_tabs[["Univariable, clinical"]] <- data.frame(Model = "Univariable", format_cox(uni_cov), check.names = FALSE)
    out$models <- c(out$models, sprintf("Univariable (clinical, complete case per term): Surv(%s, %s) ~ <term>; terms: %s",
                                        time, event, paste(uni_terms, collapse = ", ")))
  }

  m <- build_multi(de, c(gt, cov_terms), event = event)
  add_model("Multivariable", m, keep_table = is.null(g$extra_covariate) || pop != "All")
  cat(sprintf("  [%s] 다변량: n = %d (결측 제외 %d), 사건 %d, 모수 %d, EPV %.1f, stage = %s%s\n",
              tag, nrow(m$d), m$dropped, m$events, m$k, m$epv, m$stage_levels,
              if (m$exploratory) " → 탐색적(EPV 부족)" else ""))

  # ---- + extra covariate (MMR): 같은 환자에서 보정 전/후 비교 + 상호작용 ----
  x <- g$extra_covariate
  if (!is.null(x) && pop == "All" && x %in% names(de)) {
    xl  <- model_label_extra(g)
    mk  <- build_multi(de[!is.na(de[[x]]), ], c(gt, cov_terms), event = event)
    add_model(paste0("Multivariable (", xl, " known)"), mk)
    mx  <- build_multi(de, c(gt, cov_terms, x), event = event)
    fx  <- add_model(paste0("Multivariable + ", xl), mx, keep_table = TRUE)
    fi  <- fit_cox(mx$d, c(mx$terms, paste0(gt, ":", x)), time = time, event = event)
    p_int <- anova(fx, fi)[["Pr(>|Chi|)"]][2]
    rows[[length(rows)]]$p_interaction <- p_int
    out$models <- c(out$models, paste0("Interaction LRT: anova(+ ", xl, " model, ",
                                       paste(deparse(formula(fi), width.cutoff = 500L), collapse = " "), ")"))
    cat(sprintf("  [%s] + %s: n = %d, 사건 %d, EPV %.1f; %s × %s 상호작용 p = %s\n",
                tag, xl, nrow(mx$d), mx$events, mx$epv, gt, x, fmt_p(p_int)))
  }

  # ---- 시간 분할 Cox (04와 같은 규칙: 유전자 항이 PH 위반일 때만, 다변량 모형 기준) ----
  zz <- bind_rows(zph)
  if (any(zz$flag & zz$term == gt)) {
    cut <- cfg$ph_split_months
    tsc <- time_split_cox(m$d, m$terms, cut, time = time, event = event)
    cat("  [", tag, "] ", tsc$msg, "\n", sep = "")
    out$models <- c(out$models, sprintf("Time-split Cox (gene-term PH violation): Multivariable split at %d months", cut))
    for (w in 1:2) {
      sfx <- c("early", "late")[w]
      rows[[length(rows) + 1]] <- data.frame(
        endpoint = ep, population = pop,
        model = if (w == 1) sprintf("Multivariable, 0–%d months", cut) else sprintf("Multivariable, >%d months", cut),
        term = gt,
        n = if (w == 1) nrow(m$d) else sum(m$d[[time]] > cut), events = sum(tsc$events[w, ]),
        dropped = m$dropped, params = m$k, epv = m$epv, stage_levels = m$stage_levels,
        exploratory = m$exploratory,
        hr = tsc$est[[paste0("hr_", sfx)]], lo = tsc$est[[paste0("lo_", sfx)]],
        hi = tsc$est[[paste0("hi_", sfx)]], p = tsc$est[[paste0("p_", sfx)]],
        ph_p_term = NA, ph_p_global = NA, stringsAsFactors = FALSE)
    }
  }

  # ---- RMST 차이 (High − Low, 비례위험 가정 불필요) ----
  rr <- rmst_diff(de, time = time, event = event)
  if (nzchar(rr$note)) warning(gse, " ", tag, ": ", rr$note, call. = FALSE)
  cat(sprintf("  [%s] RMST 차이 (%s − %s, tau %g개월): %.2f (%.2f, %.2f), p %s\n", tag, group_alt(), group_ref(),
              rr$est$rmst_tau, rr$est$rmst_diff, rr$est$rmst_lo, rr$est$rmst_hi, fmt_p(rr$est$rmst_p)))

  out$summary <- cbind(bind_rows(rows), rr$est)
  out$summary$logrank_p <- p_lr
  out$cox <- bind_rows(cox_tabs)
  out$zph <- bind_rows(zph)
  out
}

# TCGA 결과 (04의 survival_summary_raw.csv) — 방향 비교용
tcga_reference <- function(cancer, gene) {
  f <- file.path(cfg$output_dir, "tables", cancer, out_file(cancer, "survival_summary_raw.csv"))
  if (!file.exists(f)) {
    warning(f, " 없음 → TCGA 방향 비교 생략 (04_survival.R 먼저 실행)", call. = FALSE)
    return(NULL)
  }
  s <- read.csv(f, stringsAsFactors = FALSE)
  r <- s[s$gene == gene, ]
  if (!nrow(r) || is.na(r$multi_hr)) return(NULL)
  # 04가 다른 cfg$cox_gene_term으로 실행됐으면 HR 척도가 달라 비교 불가 (열 없는 옛 파일 = group)
  tg <- r$cox_gene_term %||% "group"
  if (tg != (cfg$cox_gene_term %||% "group")) {
    warning(f, ": cox_gene_term = ", tg, " (현재 cfg$cox_gene_term = ", cfg$cox_gene_term,
            ") → TCGA 방향 비교 생략 (04_survival.R 다시 실행)", call. = FALSE)
    return(NULL)
  }
  r
}

# ---- 데이터셋 1개 ------------------------------------------------------

run_geo_validation <- function(gse, cfg_geo) {
  g    <- cfg_geo
  gene <- toupper(cfg$primary_gene)
  created <- character()
  cat("\n==========", gse, "(", g$platform, ", 검증 대상 TCGA-", g$validates, ") ==========\n", sep = "")

  coh  <- geo_load_cohort(gse, g)   # series matrix, log2 확인, phenotype 매핑/재코딩 (R/geo.R)
  eset <- coh$eset
  ex   <- coh$ex
  d    <- coh$d

  # ---- probe 선택 (종양 샘플 기준) ----
  pq <- probe_qc(ex[, d$sample_id, drop = FALSE], fData(eset), gene, g, gse)
  cat("\n[", gene, " probe] ", nrow(pq), "개, 선택 규칙 cfg$geo_probe_rule = ", cfg$geo_probe_rule, "\n", sep = "")
  print(transform(pq, mean = round(mean, 3), median = round(median, 3), IQR = round(IQR, 3), sd = round(sd, 3)),
        row.names = FALSE, digits = 3)
  f_pq <- save_table(pq, "GEO", "probe_QC.csv", subdir = gse, prefix = gse)
  created <- c(created, f_pq)
  probe <- pq$probe_id[pq$chosen]

  # ---- median split (종양 샘플 전체) ----
  d$expr_log2 <- unname(ex[probe, d$sample_id])
  d$group <- median_split(d$expr_log2)
  cat(sprintf("\nmedian split (%s, 종양 %d명): median %.3f → %s %d / %s %d, 발현 결측 %d\n",
              probe, sum(!is.na(d$expr_log2)), median(d$expr_log2, na.rm = TRUE),
              group_ref(), sum(d$group == group_ref(), na.rm = TRUE),
              group_alt(), sum(d$group == group_alt(), na.rm = TRUE), sum(is.na(d$expr_log2))))

  # ---- pheno QC ----
  pq_tab <- pheno_qc(d, g)
  for (ep in names(g$endpoints)) {
    ok <- !is.na(d[[paste0(ep, "_months")]]) & !is.na(d[[paste0(ep, "_event")]])
    pq_tab <- rbind(pq_tab, data.frame(
      field = paste0(ep, " (usable)"), raw_key = "", n = sum(ok), n_missing = sum(!ok),
      values = paste0("events = ", sum(d[[paste0(ep, "_event")]][ok])), stringsAsFactors = FALSE))
  }
  cat("\n[phenotype QC] 종양", nrow(d), "명\n")
  print(pq_tab, right = FALSE, row.names = FALSE)
  f_ph <- save_table(pq_tab, "GEO", "pheno_QC.csv", subdir = gse, prefix = gse)
  created <- c(created, f_ph)

  # ---- 생존분석: 종점 × (전체, 하위군) ----
  pops <- list(All = d)
  sg <- g$subgroup
  if (!is.null(sg)) pops[[sg$level]] <- d[!is.na(d[[sg$var]]) & d[[sg$var]] == sg$level, ]
  res <- list()
  for (ep in names(g$endpoints)) {
    for (pop in names(pops)) {
      cat("\n--", ep, "/", pop, "--\n")
      r <- analyze_endpoint(pops[[pop]], ep, g, gse, gene, pop)
      res[[paste(ep, pop)]] <- r
      created <- c(created, r$files)
    }
  }

  # ---- TCGA 방향 비교 ----
  tcga <- tcga_reference(g$validates, gene)

  # ---- 분석 정보용 공통 값 (출력 CSV별 → <GSE>_analysis_info_06.csv) ----
  ds_set <- function(k) {   # cfg$geo_datasets$<GSE>$<k> (이름 있으면 이름:값)
    v <- unlist(g[[k]])
    if (!is.null(names(v))) v <- paste0(names(v), ":", v)
    paste0("geo_datasets$", gse, "$", k, "=", paste(v, collapse = ","))
  }
  geo_settings <- c("primary_gene", "group_levels", "cox_gene_term", "cox_covariates", "cox_uni_covariates",
                    "epv_min", "stage_full", "stage_collapsed", "min_events", "ph_split_months", "rmst_tau",
                    "km_times", ds_set("cox_covariates_drop"), ds_set("extra_covariate"), ds_set("subgroup"),
                    paste0("exclude_stage=", paste(sprintf("%s:%s", names(g$endpoints), vapply(g$endpoints, function(e)
                      paste(e$exclude_stage %||% "none", collapse = ","), "")), collapse = ",")))
  tcga_f <- out_file(g$validates, "survival_summary_raw.csv")
  geo_data <- c(sprintf("%s / %s (series matrix as deposited, GEOquery; %d tumors)", gse, annotation(eset), nrow(d)),
                sprintf("probe %s (%s)", probe, cfg$geo_probe_rule),
                if (is.null(tcga)) "TCGA direction reference: none"
                else paste0("TCGA direction reference: ", tcga_f, " (multi_hr ", signif(tcga$multi_hr, 4), ")"))
  n_all <- vapply(res, function(r) r$n_txt %||% NA_character_, "")
  n_all <- paste(n_all[!is.na(n_all)], collapse = " | ")
  all_models <- unlist(lapply(names(res), function(k) if (length(res[[k]]$models)) paste0("[", k, "] ", res[[k]]$models)))
  s <- bind_rows(lapply(res, `[[`, "summary"))
  if (!"p_interaction" %in% names(s)) s$p_interaction <- NA
  s$tcga_multi_hr <- if (is.null(tcga)) NA else tcga$multi_hr
  s$direction_match <- if (is.null(tcga)) NA else sign(log(s$hr)) == sign(log(tcga$multi_hr))
  s$gene          <- gene                              # 90_export.R가 이 유전자의 결과인지 확인
  s$cox_gene_term <- cfg$cox_gene_term %||% "group"
  f_raw <- save_table(s, "GEO", "survival_summary_raw.csv", subdir = gse, prefix = gse)
  created <- c(created, f_raw)

  ep_lab <- vapply(g$endpoints, `[[`, "", "label")
  tab <- data.frame(
    Endpoint = unname(ep_lab[s$endpoint]), Population = s$population, Model = s$model,
    N = s$n, Events = s$events,
    `HR (95% CI)` = fmt_hr(s$hr, s$lo, s$hi), p = fmt_p(s$p),
    EPV = sprintf("%.1f", s$epv), `Stage variable` = s$stage_levels,
    `PH p (gene term)` = fmt_p(s$ph_p_term), `Interaction p` = fmt_p(s$p_interaction),
    `Log-rank p` = ifelse(duplicated(paste(s$endpoint, s$population)), "", fmt_p(s$logrank_p)),
    check.names = FALSE)
  first <- !duplicated(paste(s$endpoint, s$population))
  tab[[paste0("RMST difference at ", cfg$rmst_tau, " mo, months (95% CI)")]] <-
    ifelse(first, sprintf("%.2f (%.2f, %.2f)", s$rmst_diff, s$rmst_lo, s$rmst_hi), "")
  tab[["RMST p"]] <- ifelse(first, fmt_p(s$rmst_p), "")
  tab[[paste0("Direction as TCGA-", g$validates)]] <- ifelse(is.na(s$direction_match), "",
                                                             ifelse(s$direction_match, "Yes", "No"))
  tab[[paste0("Exploratory (EPV < ", cfg$epv_min, ")")]] <- s$exploratory
  created <- c(created, save_df_table(
    tab, "GEO", "survival_summary", landscape = TRUE, font_size = 7, subdir = gse, prefix = gse,
    caption = paste0(
      "External validation of ", gene, " expression (", group_contrast_label(), ", median split of probe ", probe,
      " within tumors) in ", g$label, ". HR from Cox regression",
      if (gene_term() != "group") "; multivariable and time-split gene HRs are per 1 log2 unit of expression" else "",
      "; multivariable models adjusted for ", adjust_text(setdiff(cfg$cox_covariates, g$cox_covariates_drop)),
      if (!is.null(g$extra_covariate)) paste0("; '+ ", model_label_extra(g), "' additionally adjusts for ",
                                              tolower(model_label_extra(g)), " and '(", model_label_extra(g),
                                              " known)' is the model without it on the same patients") else "",
      if (!is.null(g$extra_covariate)) paste0(". Interaction p = likelihood-ratio test for ", gene,
                                              if (gene_term() == "group") " group × " else " expression × ",
                                              tolower(model_label_extra(g))) else "",
      ". Time-split HRs (0–", cfg$ph_split_months, " / >", cfg$ph_split_months,
      " months, adjusted) are shown only when the gene term violated proportional hazards. RMST difference = ",
      group_alt(), " − ", group_ref(), " restricted mean survival time up to ", cfg$rmst_tau, " months. ",
      paste(vapply(names(g$endpoints), function(ep) {
        e <- g$endpoints[[ep]]
        if (is.null(e$exclude_stage)) "" else paste0(e$label, " excludes stage ", paste(e$exclude_stage, collapse = "/"), ". ")
      }, ""), collapse = ""),
      if (!is.null(tcga)) sprintf("Direction compared with the TCGA-%s multivariable HR %.2f.", g$validates, tcga$multi_hr) else "")))

  # ---- KM 요약 ----
  km <- bind_rows(lapply(res, `[[`, "km"))
  km_out <- data.frame(Endpoint = unname(ep_lab[km$endpoint]), Population = km$population,
                       Group = km$group, N = km$n, Events = km$events,
                       `Median, months (95% CI)` = fmt_median(km$median, km$median_lo, km$median_hi),
                       check.names = FALSE)
  for (tm in cfg$km_times) {
    col <- paste0("surv_", tm / 12, "y")
    km_out[[paste0(tm / 12, "-year (95% CI)")]] <- fmt_surv(km[[col]], km[[paste0(col, "_lo")]], km[[paste0(col, "_hi")]])
  }
  km_out$`Log-rank p` <- ifelse(duplicated(paste(km$endpoint, km$population)), "", fmt_p(km$logrank_p))
  created <- c(created, save_df_table(km_out, "GEO", "km_summary", subdir = gse, prefix = gse,
                                      caption = paste0("Kaplan–Meier estimates by ", gene, " expression, ",
                                                       g$label, ". NR = not reached.")))

  # ---- Cox 표 (종점 × 집단) ----
  cox_info <- list()
  for (k in names(res)) {
    r <- res[[k]]
    if (is.null(r$cox)) next
    ep  <- sub(" .*$", "", k); pop <- sub("^\\S+ ", "", k)
    stem <- paste0("cox_", ep, if (pop != "All") paste0("_", pop))
    f <- save_df_table(
      r$cox, "GEO", stem, subdir = gse, prefix = gse,
      caption = paste0("Cox regression for ", tolower(ep_lab[[ep]]), ": ", gene, " expression, ", g$label,
                       if (pop != "All") paste0(", ", pop, " only") else ""))
    created <- c(created, f)
    cox_info[[length(cox_info) + 1]] <- list(
      Output = f[grepl("\\.csv$", f)],
      Analysis = paste0("Cox regression, ", ep_lab[[ep]], ", ", pop, " (all models, all covariate rows)"),
      Model = r$models, Settings = geo_settings, n = r$n_txt, Data = geo_data, Packages = "survival")
  }
  zph <- bind_rows(lapply(res, `[[`, "zph"))
  f_zph <- save_table(transform(zph, p = signif(p, 4), chisq = round(chisq, 3)),
                      "GEO", "ph_tests.csv", subdir = gse, prefix = gse)
  created <- c(created, f_zph)
  if (any(zph$flag)) {
    fl <- zph[zph$flag, ]
    cat("\ncox.zph p < 0.05:", paste(sprintf("%s/%s %s: %s", fl$endpoint, fl$population, fl$model, fl$term),
                                    collapse = "; "), "\n")
  }

  # ---- forest: 모형별 유전자 HR ----
  fs <- s[s$term == gene_term(), ]
  ft <- data.frame(
    label = sprintf("%s, %s: %s (n = %d, events = %d)", fs$endpoint, fs$population, fs$model, fs$n, fs$events),
    hr = fs$hr, lo = fs$lo, hi = fs$hi, p = fs$p, reference = FALSE,
    shape = ifelse(fs$exploratory, "exploratory", "est"))
  created <- c(created, save_fig(
    forest_plot(ft, paste0(gene, " ", gene_term_label(), " — ", g$label),
                paste0("Probe ", probe, "; multivariable models adjusted for ", adjust_text(setdiff(cfg$cox_covariates, g$cox_covariates_drop)))),
    "GEO", "forest", width = 13, height = 1.6 + 0.35 * nrow(ft), subdir = gse, prefix = gse))

  # 전체 공변량 forest (종점별, 가장 많이 보정한 모형)
  for (ep in names(g$endpoints)) {
    fo <- res[[paste(ep, "All")]]$forest
    if (is.null(fo)) next
    created <- c(created, save_fig(
      forest_plot(fo$tab, paste0(fo$name, " Cox: ", ep_lab[[ep]], " — ", g$label),
                  sprintf("n = %d, events = %d, EPV = %.1f", nrow(fo$m$d), fo$m$events, fo$m$epv)),
      "GEO", paste0("forest_multi_", ep), width = 9, height = 1.6 + 0.35 * nrow(fo$tab), subdir = gse, prefix = gse))
  }

  # ---- 분석 정보 ----
  tab_path <- function(name) file.path(out_dir("GEO", "tables", gse), out_file(gse, name))
  surv_model <- c("KM + log-rank: Surv(<EP>_months, <EP>_event) ~ group", all_models,
                  paste0("RMST: survRM2::rmst2 (", group_contrast_label(), ")"), "cox.zph for every model",
                  paste0("direction vs TCGA-", g$validates, " multivariable HR"))
  info <- c(
    list(
      list(Output = f_pq, Analysis = paste0("Probe QC and selection for ", gene),
           Model = c(paste0("all probes with ", g$symbol_col, " = ", gene, " (split on ///)"),
                     paste0("chosen by geo_probe_rule over tumors; median split of probe ", probe, " within tumors")),
           Settings = c("geo_probe_rule", ds_set("symbol_col"), ds_set("platform")),
           n = sprintf("%d probes; %d tumors", nrow(pq), nrow(d)), Data = geo_data, Packages = c("GEOquery", "Biobase")),
      list(Output = f_ph, Analysis = "Phenotype mapping QC (characteristics_ch1 → canonical fields)",
           Model = "n, missing and values per mapped field; usable n/events per endpoint",
           Settings = c("geo_na_values", ds_set("fields"), ds_set("exclude_sample_type"), ds_set("stage_na_values")),
           n = sprintf("%d tumors", nrow(d)), Data = geo_data, Packages = c("GEOquery", "Biobase")),
      list(Output = f_raw, Analysis = "External validation summary per endpoint × population × model",
           Model = surv_model, Settings = geo_settings, n = n_all, Data = geo_data,
           Packages = c("survival", "survRM2", "GEOquery", "Biobase")),
      list(Output = tab_path("survival_summary.csv"),
           Analysis = "Formatted external validation summary (from survival_summary_raw)",
           Model = surv_model, Settings = geo_settings, n = n_all, Data = geo_data,
           Packages = c("survival", "survRM2", "GEOquery", "Biobase")),
      list(Output = tab_path("km_summary.csv"), Analysis = "Kaplan–Meier estimates by expression group",
           Model = c("survfit(Surv(<EP>_months, <EP>_event) ~ group)", "survdiff (log-rank)"),
           Settings = c("group_levels", "km_times", ds_set("subgroup")), n = n_all, Data = geo_data,
           Packages = "survival")),
    cox_info,
    list(list(Output = f_zph, Analysis = "Proportional hazards test (Schoenfeld residuals)",
              Model = paste0("cox.zph of: ", paste(all_models[!grepl("^\\[[^]]+\\] (Univariable \\(clinical|Interaction|Time-split)", all_models)],
                                                   collapse = "; ")),
              Settings = c("cox_gene_term", "cox_covariates"), n = n_all, Data = geo_data, Packages = "survival")))
  created <- c(created, write_analysis_info(info, "GEO", "06", subdir = gse, prefix = gse))

  move_stale_outputs("GEO", c("km_", "forest", "cox_", "zph_"), created, subdir = gse, prefix = gse)
  list(files = created, summary = s, probe = pq[pq$chosen, ], tcga = tcga, g = g, d = d)
}

# ---- 실행 -------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
datasets <- if (length(args)) args else names(cfg$geo_datasets)
unknown <- setdiff(datasets, names(cfg$geo_datasets))
if (length(unknown)) stop("cfg$geo_datasets에 없는 GSE: ", paste(unknown, collapse = ", "))

created  <- character()
failures <- character()
console  <- list()
for (gse in datasets) {
  r <- tryCatch(run_geo_validation(gse, cfg$geo_datasets[[gse]]), error = function(e) {
    failures <<- c(failures, paste0(gse, ": ", conditionMessage(e)))
    cat("  !!", gse, "실패:", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(r)) {
    created <- c(created, r$files)
    console[[gse]] <- r
  }
}

cat("\n생성된 파일 (", length(created), "개):\n", paste0("  ", created, collapse = "\n"), "\n", sep = "")

cat("\n===== 요약 =====\n")
for (gse in names(console)) {
  r <- console[[gse]]; s <- r$summary; g <- r$g
  cat(sprintf("%s (검증 대상 TCGA-%s), %s probe %s (%s; mean %.2f, IQR %.2f)\n", gse, g$validates,
              toupper(cfg$primary_gene), r$probe$probe_id, r$probe$rule, r$probe$mean, r$probe$IQR))
  if (!is.null(r$tcga)) {
    cat(sprintf("  TCGA-%s 다변량 HR %s, p %s\n", g$validates,
                fmt_hr(r$tcga$multi_hr, r$tcga$multi_lo, r$tcga$multi_hi), fmt_p(r$tcga$multi_p)))
  }
  for (i in which(s$term == gene_term())) {
    cat(sprintf("  %-4s %-5s %-38s n = %3d, 사건 %3d: HR %s, p %s%s%s\n",
                s$endpoint[i], s$population[i], s$model[i], s$n[i], s$events[i],
                fmt_hr(s$hr[i], s$lo[i], s$hi[i]), fmt_p(s$p[i]),
                if (is.na(s$direction_match[i])) "" else if (s$direction_match[i]) "  [TCGA와 같은 방향]" else "  [TCGA와 반대 방향]",
                if (s$exploratory[i]) " (탐색적)" else ""))
  }
  lr <- unique(s[, c("endpoint", "population", "logrank_p")])
  cat("  log-rank:", paste(sprintf("%s/%s p %s", lr$endpoint, lr$population, fmt_p(lr$logrank_p)), collapse = ", "), "\n")
}
if (length(failures)) {
  cat("\n실패한 분석 (", length(failures), "개):\n", paste0("  ", failures, collapse = "\n"), "\n", sep = "")
} else {
  cat("\n실패한 분석 없음\n")
}
