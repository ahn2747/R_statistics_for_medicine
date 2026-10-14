# =============================================================
# 13_sig_geo.R [GSE ...] [SIG ...]
# signature 점수의 GEO 외부 검증 (데이터셋: cfg$geo_datasets, signature: cfg$signatures)
#   1. R/geo.R: series matrix + phenotype (06과 같은 정리), 종양 샘플
#   2. 전 유전자 probe 축약 (cfg$geo_probe_rule, 06의 probe 선택 규칙과 같음) → 유전자 × 종양 행렬
#   3. 저발현 필터 (cfg$signature_expr_filter$geo) → 점수 4종 + QC (11과 같은 함수; z/PC1은 이 코호트 안에서 계산)
#   4. 방법별 median split (전체 종양, 06과 같음) → 종점 × (전체, 하위군):
#      KM/log-rank/RMST, Cox per 1 SD (단변량, 다변량 EPV), + extra_covariate 모형과 상호작용 (06과 같음), cox.zph, 시간 분할
#      1 SD = 각 종점 × 집단 분석 코호트 안의 SD
#   5. TCGA signature 결과(12)와 HR 방향만 비교 (점수 값 자체는 비교하지 않음)
# 결과 (파일명 앞에 <GSE>_): data/processed/<GSE>_sigscores_<SIG>.csv
#   output/tables/GEO/<GSE>/signature/<SIG>/  sig_genes_, sig_qc_, sig_score_cor_, sig_survival_summary(_raw)_, sig_cox_<SIG>_<EP>[_<sub>],
#                                             sig_km_summary_, sig_ph_tests_, analysis_info_13
#   output/figures/GEO/<GSE>/signature/<SIG>/ sig_score_pairs_, sig_random_pc1_, sig_km_<SIG>_<method>_<EP>[_<sub>], sig_forest_, sig_forest_multi_<SIG>_<EP>
# =============================================================

source("config.R")
source("R/utils.R")
source("R/survival.R")
source("R/geo.R")
source("R/signature.R")

args <- commandArgs(trailingOnly = TRUE)
datasets <- if (any(args %in% names(cfg$geo_datasets))) args[args %in% names(cfg$geo_datasets)] else names(cfg$geo_datasets)
sigs <- args[!args %in% names(cfg$geo_datasets)]
all_sigs <- !length(sigs)
if (all_sigs) sigs <- names(cfg$signatures)
if (!length(sigs)) stop("cfg$signatures가 비어 있음 → config.R에 signature 정의 추가")
unknown <- setdiff(sigs, names(cfg$signatures))
if (length(unknown)) stop("cfg$signatures 또는 cfg$geo_datasets에 없는 이름: ", paste(unknown, collapse = ", "))
pm <- cfg$signature_primary_method
if (length(pm) != 1 || is.na(pm) || !pm %in% sig_methods) {
  stop("cfg$signature_primary_method는 ", paste(sig_methods, collapse = ", "), " 중 하나 (결과를 보기 전에 지정)")
}
defs <- lapply(setNames(sigs, sigs), sig_def)
cat("GSE:", paste(datasets, collapse = ", "), "/ signature:", paste(sigs, collapse = ", "), "/ 메인 방법:", pm, "\n")

fmt_rmst <- function(x, lo, hi) ifelse(is.na(x), NA, sprintf("%.2f (%.2f, %.2f)", x, lo, hi))

# TCGA signature 결과 (12의 sig_survival_summary_raw) → 방법별 다변량 HR
tcga_sig_reference <- function(cancer, sig) {
  f <- file.path(cfg$output_dir, "tables", cancer, sig_subdir(sig), out_file(cancer, paste0("sig_survival_summary_raw_", sig, ".csv")))
  if (!file.exists(f)) {
    warning(f, " 없음 → TCGA 방향 비교 생략 (12_sig_survival.R ", sig, " ", cancer, " 먼저 실행)", call. = FALSE)
    return(NULL)
  }
  s <- read.csv(f, stringsAsFactors = FALSE)
  s <- s[!is.na(s$model) & s$model == "Multivariable (per 1 SD)", c("method", "hr")]
  list(file = f, hr = setNames(s$hr, s$method))
}

# + extra_covariate 모형 (06과 같은 구성): 같은 환자에서 보정 전/후 + score × extra 상호작용 LRT
extra_models <- function(de, cov_terms, x, xl, labs, time, event, label) {
  rows <- list(); cox <- list(); zph <- list(); models <- character()
  add <- function(name, m) {
    fit <- fit_cox(m$d, m$terms, time = time, event = event)
    models <<- c(models, sprintf("%s: %s (n=%d, events=%d)", name,
                                 paste(deparse(formula(fit), width.cutoff = 500L), collapse = " "), nrow(m$d), m$events))
    tab <- cox_rows(fit, m$d, m$terms, labs[m$terms], event = event)
    z <- zph_rows(fit, label, name)
    k <- tab[tab$term == "score_sd", ]
    rows[[name]] <<- data.frame(model = name, term = "score_sd", n = nrow(m$d), events = m$events, dropped = m$dropped,
                                params = m$k, epv = m$epv, stage_levels = m$stage_levels %||% NA, exploratory = m$exploratory,
                                hr = k$hr, lo = k$lo, hi = k$hi, p = k$p,
                                ph_p_term = z$p[z$term == "score_sd"], ph_p_global = z$p[z$term == "GLOBAL"],
                                stringsAsFactors = FALSE)
    cox[[name]] <<- data.frame(Model = name, format_cox(tab), check.names = FALSE)
    zph[[name]] <<- z
    list(fit = fit, tab = tab, m = m)
  }
  mk <- build_multi(de[!is.na(de[[x]]), ], c("score_sd", cov_terms), event = event)
  add(paste0("Multivariable (", xl, " known) (per 1 SD)"), mk)
  mx <- build_multi(de, c("score_sd", cov_terms, x), event = event)
  fx <- add(paste0("Multivariable + ", xl, " (per 1 SD)"), mx)
  fi <- fit_cox(mx$d, c(mx$terms, paste0("score_sd:", x)), time = time, event = event)
  rows[[length(rows)]]$p_interaction <- anova(fx$fit, fi)[["Pr(>|Chi|)"]][2]
  models <- c(models, paste0("Interaction LRT: anova(+ ", xl, " model, ", paste(deparse(formula(fi), width.cutoff = 500L), collapse = " "), ")"))
  list(rows = dplyr::bind_rows(rows), cox = dplyr::bind_rows(cox), zph = dplyr::bind_rows(zph), models = models, most = fx)
}

created  <- character()
failures <- character()
console  <- list()

for (gse in datasets) {
  g <- cfg$geo_datasets[[gse]]
  cat("\n==========", gse, "(", g$platform, ", 비교 대상 TCGA-", g$validates, ") ==========\n", sep = "")
  coh <- tryCatch(geo_load_cohort(gse, g), error = function(e) {
    failures <<- c(failures, paste0(gse, ": ", conditionMessage(e))); NULL
  })
  if (is.null(coh)) next
  d  <- coh$d
  gm <- geo_gene_matrix(coh$ex[, d$sample_id, drop = FALSE], fData(coh$eset), g, gse)
  rule <- cfg$signature_expr_filter$geo
  flt <- sig_filter_genes(gm$mat, rule)
  cat(sprintf("종양 %d명; 저발현 필터 (%s): %d → %d유전자\n", nrow(d), sig_filter_text(rule), nrow(gm$mat), nrow(flt$expr)))
  ranks <- singscore::rankGenes(flt$expr)
  ep_lab <- vapply(g$endpoints, `[[`, "", "label")
  xcov <- g$extra_covariate
  xl <- if (is.null(xcov)) NULL else (var_labels[[xcov]] %||% xcov)
  geo_data <- c(sprintf("%s / %s (series matrix as deposited, GEOquery; %d tumors)", gse, annotation(coh$eset), nrow(d)),
                sprintf("all-gene probe collapse by geo_probe_rule = %s (%d symbols)", cfg$geo_probe_rule, nrow(gm$mat)))

  res <- list()
  for (sig in sigs) {
    def <- defs[[sig]]
    sub <- file.path(gse, sig_subdir(sig))
    cat("\n--", sig, "--\n")
    r <- tryCatch({
      # ---- 점수 + QC ----
      tab  <- sig_match(def, rownames(gm$mat), rownames(flt$expr), gse, missing_label = "not_on_platform")
      used <- sig_used_def(def, tab)
      lost <- tab[tab$status != "matched", ]
      cat(sprintf("  사용 유전자: up %d/%d, down %d/%d%s\n", length(used$up), length(def$up), length(used$down), length(def$down),
                  if (nrow(lost)) paste0(" (제외: ", paste(sprintf("%s %s", lost$gene, lost$status), collapse = ", "), ")") else ""))
      scores <- sig_scores(flt$expr, used, ranks)
      qc <- sig_qc(flt$expr, used, tab, flt$stats, scores)
      pmap <- gm$probe_map[match(qc$genes$gene, gm$probe_map$symbol), c("probe_id", "n_probes", "multi_mapped")]
      qc$genes <- cbind(qc$genes, pmap, row.names = NULL)
      qc$summary <- cbind(qc$summary, random_n = cfg$signature_random_n, seed = cfg$seed, n_samples = nrow(d),
                          n_genes_background = nrow(flt$expr), filter_rule = sig_filter_text(rule))

      # ---- 생존: 방법 × 종점 × 집단 ----
      for (m in sig_methods) d[[paste0(m, "_group")]] <- median_split(scores[[m]][match(d$sample_id, scores$sample_id)])
      for (m in sig_methods) d[[m]] <- scores[[m]][match(d$sample_id, scores$sample_id)]
      pops <- list(All = d)
      sg <- g$subgroup
      if (!is.null(sg)) pops[[sg$level]] <- d[!is.na(d[[sg$var]]) & d[[sg$var]] == sg$level, ]
      tcga <- tcga_sig_reference(g$validates, sig)
      out <- list()
      for (ep in names(g$endpoints)) for (pop in names(pops)) for (m in sig_methods) {
        time <- paste0(ep, "_months"); event <- paste0(ep, "_event")
        tag <- paste0(if (pop == "All") ep else paste0(ep, "_", pop), " ", m)
        dp <- pops[[pop]]
        dp$group <- dp[[paste0(m, "_group")]]
        de <- add_scaled_terms(geo_endpoint_data(dp, ep, g, gse, tag, keep = !is.na(dp[[m]])))
        sd_raw <- stats::sd(de[[m]])
        de$score_sd <- (de[[m]] - mean(de[[m]])) / sd_raw
        cov_terms <- covariate_terms(setdiff(cfg$cox_covariates, g$cox_covariates_drop), de,
                                     paste0("geo_datasets$", gse, "$cox_covariates_drop / cox_covariates"), paste(gse, tag))
        label <- paste(def$label, sig_method_labels[[m]])
        ss <- sig_survival_set(de, label, cov_terms, time = time, event = event, where = paste(gse, sig, tag))
        if (!is.null(ss$rows) && !is.null(xcov) && pop == "All" && xcov %in% names(de)) {
          xm <- extra_models(de, cov_terms, xcov, xl, c(ss$labels, term_labels(xcov)), time, event, label)
          ss$rows <- dplyr::bind_rows(ss$rows, xm$rows)
          ss$cox <- dplyr::bind_rows(ss$cox, xm$cox)
          ss$zph <- dplyr::bind_rows(ss$zph, xm$zph)
          ss$models <- c(ss$models, xm$models)
          ss$multi <- list(m = xm$most$m, tab = xm$most$tab)   # forest: 가장 많이 보정한 모형
        }
        ss$score_sd_raw <- sd_raw
        ss$data <- de
        out[[paste(ep, pop, m)]] <- c(ss, list(ep = ep, pop = pop, method = m))
      }
      list(def = def, sub = sub, tab = tab, used = used, scores = scores, qc = qc, out = out, tcga = tcga)
    }, error = function(e) {
      msg <- conditionMessage(e)
      failures <<- c(failures, if (startsWith(msg, paste(gse, sig))) msg else paste0(gse, " ", sig, ": ", msg))
      cat("  !!", sig, "실패:", msg, "\n")
      NULL
    })
    if (!is.null(r)) res[[sig]] <- r
  }

  # ---- raw 요약 + signature 간 BH (메인 방법, 전체 실행, signature 2개 이상) ----
  raw_all <- lapply(res, function(r) dplyr::bind_rows(lapply(r$out, function(o) {
    if (is.null(o$rows)) return(data.frame(endpoint = o$ep, population = o$pop, method = o$method, note = paste("skipped:", o$skipped)))
    x <- data.frame(cohort = gse, sig = r$def$sig, method = o$method, primary = o$method == pm, endpoint = o$ep,
                    population = o$pop, o$rows, logrank_p = o$p_lr, o$rmst$est, n_low = o$n_low, n_high = o$n_high,
                    score_sd_raw = o$score_sd_raw, note = o$rmst$note, stringsAsFactors = FALSE)
    x$tcga_multi_hr <- if (is.null(r$tcga)) NA else unname(r$tcga$hr[o$method])
    x$direction_match <- ifelse(is.na(x$tcga_multi_hr), NA, sign(log(x$hr)) == sign(log(x$tcga_multi_hr)))
    x
  })))
  do_bh <- all_sigs && length(raw_all) > 1
  for (sig in names(raw_all)) {
    raw_all[[sig]]$q_across_signatures <- NA_real_
    if (!"p_interaction" %in% names(raw_all[[sig]])) raw_all[[sig]]$p_interaction <- NA_real_
  }
  if (do_bh) {
    prim <- dplyr::bind_rows(lapply(raw_all, function(x) x[x$primary %in% TRUE, c("sig", "endpoint", "population", "model", "p")]))
    key <- paste(prim$endpoint, prim$population, prim$model)
    prim$q <- ave(prim$p, key, FUN = function(p) p.adjust(p, "BH"))
    for (sig in names(raw_all)) {
      x <- raw_all[[sig]]; i <- which(x$primary %in% TRUE)
      x$q_across_signatures[i] <- prim$q[match(paste(sig, x$endpoint[i], x$population[i], x$model[i]),
                                               paste(prim$sig, key))]
      raw_all[[sig]] <- x
    }
  }

  # ---- 출력 ----
  for (sig in names(res)) {
    r <- res[[sig]]; def <- r$def; sub <- r$sub; raw <- raw_all[[sig]]
    sv <- function(df, name) save_table(df, "GEO", name, subdir = sub, prefix = gse)
    files <- character()
    tryCatch({
      cov_txt <- adjust_text(setdiff(cfg$cox_covariates, g$cox_covariates_drop))
      sf <- sig_scores_path(gse, sig)
      readr::write_csv(r$scores, sf)
      f_genes <- sv(r$qc$genes, paste0("sig_genes_", sig, ".csv"))
      f_qc    <- sv(r$qc$summary, paste0("sig_qc_", sig, ".csv"))
      f_cor   <- sv(r$qc$cor, paste0("sig_score_cor_", sig, ".csv"))
      f_raw   <- sv(raw, paste0("sig_survival_summary_raw_", sig, ".csv"))
      files <- c(sf, f_genes, f_qc, f_cor, f_raw)
      files <- c(files, save_fig(sig_pairs_plot(r$scores, paste0(def$label, " signature scores — ", g$label, " (n = ", nrow(d), ")")),
                                 "GEO", paste0("sig_score_pairs_", sig), width = 9, height = 6, subdir = sub, prefix = gse))
      rp <- sig_random_plot(r$qc$null, r$qc$summary, paste0(def$label, " PC1 vs random gene sets — ", g$label))
      if (!is.null(rp)) files <- c(files, save_fig(rp, "GEO", paste0("sig_random_pc1_", sig), width = 8, height = 4,
                                                   subdir = sub, prefix = gse))

      # 요약 표 (메인 방법 먼저)
      sm <- raw[!is.na(raw$model), ]
      sm <- sm[order(match(sm$endpoint, names(g$endpoints)), sm$population != "All", sm$method != pm,
                     match(sm$method, sig_methods)), ]
      grp <- paste(sm$endpoint, sm$population, sm$method)
      first <- !duplicated(grp)
      tab <- data.frame(Endpoint = unname(ep_lab[sm$endpoint]), Population = sm$population,
                        Method = ifelse(first, paste0(sig_method_labels[sm$method], ifelse(sm$primary, " (primary)", "")), ""),
                        Model = sm$model, N = sm$n, Events = sm$events, `HR (95% CI)` = fmt_hr(sm$hr, sm$lo, sm$hi),
                        p = fmt_p(sm$p), check.names = FALSE)
      if (do_bh) tab$`q (across signatures, primary method)` <- ifelse(sm$primary, fmt_p(sm$q_across_signatures), "")
      tab$EPV <- sprintf("%.1f", sm$epv)
      tab$`PH p (score term)` <- fmt_p(sm$ph_p_term)
      tab$`Interaction p` <- fmt_p(sm$p_interaction)
      tab$`Log-rank p` <- ifelse(first, fmt_p(sm$logrank_p), "")
      tab[[paste0("RMST difference at ", cfg$rmst_tau, " mo, months (95% CI)")]] <- ifelse(first, fmt_rmst(sm$rmst_diff, sm$rmst_lo, sm$rmst_hi), "")
      tab$`RMST p` <- ifelse(first, fmt_p(sm$rmst_p), "")
      tab[[paste0("Direction as TCGA-", g$validates)]] <- ifelse(is.na(sm$direction_match), "", ifelse(sm$direction_match, "Yes", "No"))
      tab[[paste0("Exploratory (EPV < ", cfg$epv_min, ")")]] <- sm$exploratory
      files <- c(files, save_df_table(
        tab, "GEO", paste0("sig_survival_summary_", sig), landscape = TRUE, font_size = 6, subdir = sub, prefix = gse,
        caption = paste0(
          "External validation of the ", def$label, " signature score (", def$source, ") in ", g$label,
          ". Primary method: ", sig_method_labels[[pm]], " (prespecified). Scores computed within this cohort; continuous HRs per 1 SD ",
          "within each endpoint/population; High vs Low = median split over all tumors. Multivariable models adjusted for ", cov_txt,
          if (!is.null(xcov)) paste0("; '+ ", xl, "' additionally adjusts for ", tolower(xl), "; interaction p = likelihood-ratio test for score × ", tolower(xl)) else "",
          ". ", paste(vapply(names(g$endpoints), function(ep) {
            e <- g$endpoints[[ep]]
            if (is.null(e$exclude_stage)) "" else paste0(e$label, " excludes stage ", paste(e$exclude_stage, collapse = "/"), ". ")
          }, ""), collapse = ""),
          if (!is.null(r$tcga)) paste0("Direction compared with the TCGA-", g$validates, " multivariable HR of the same method.") else "")))

      # Cox 표 (종점 × 집단)
      for (ep in names(g$endpoints)) for (pop in unique(sm$population)) {
        cx <- dplyr::bind_rows(lapply(sig_methods, function(m) {
          o <- r$out[[paste(ep, pop, m)]]
          if (is.null(o$cox)) NULL else data.frame(Method = sig_method_labels[[m]], o$cox, check.names = FALSE)
        }))
        if (!nrow(cx)) next
        files <- c(files, save_df_table(cx, "GEO", paste0("sig_cox_", sig, "_", ep, if (pop != "All") paste0("_", pop)),
                                        landscape = TRUE, font_size = 7, subdir = sub, prefix = gse,
                                        caption = paste0("Cox regression for ", tolower(ep_lab[[ep]]), ": ", def$label,
                                                         " signature score, ", g$label, if (pop != "All") paste0(", ", pop, " only") else "")))
      }
      # KM 요약 + 그림
      km <- dplyr::bind_rows(lapply(r$out, function(o) if (is.null(o$km)) NULL else
        data.frame(endpoint = o$ep, population = o$pop, method = o$method, o$km, logrank_p = o$p_lr)))
      km_out <- data.frame(Endpoint = unname(ep_lab[km$endpoint]), Population = km$population,
                           Method = unname(sig_method_labels[km$method]), Group = km$group, N = km$n, Events = km$events,
                           `Median, months (95% CI)` = fmt_median(km$median, km$median_lo, km$median_hi), check.names = FALSE)
      for (tm in cfg$km_times) {
        col <- paste0("surv_", tm / 12, "y")
        km_out[[paste0(tm / 12, "-year (95% CI)")]] <- fmt_surv(km[[col]], km[[paste0(col, "_lo")]], km[[paste0(col, "_hi")]])
      }
      km_out$`Log-rank p` <- ifelse(duplicated(paste(km$endpoint, km$population, km$method)), "", fmt_p(km$logrank_p))
      files <- c(files, save_df_table(km_out, "GEO", paste0("sig_km_summary_", sig), subdir = sub, prefix = gse,
                                      caption = paste0("Kaplan–Meier estimates by ", def$label, " signature score (",
                                                       group_contrast_label(), ", median split), ", g$label, ". NR = not reached.")))
      for (o in r$out) {
        if (is.null(o$fit_km)) next
        e <- g$endpoints[[o$ep]]
        p <- km_plot(o$fit_km, o$data, paste(def$label, sig_method_labels[[o$method]]), NULL, o$p_lr,
                     title = paste0(def$label, " score (", sig_method_labels[[o$method]], ") — ", g$label,
                                    if (o$pop != "All") paste0(", ", o$pop)),
                     ylab = paste(e$label, "probability"))
        files <- c(files, suppressMessages(save_fig(p, "GEO", paste0("sig_km_", sig, "_", o$method, "_", o$ep,
                                                                     if (o$pop != "All") paste0("_", o$pop)),
                                                    width = 7, height = 6, subdir = sub, prefix = gse)))
      }
      zph <- dplyr::bind_rows(lapply(r$out, function(o) if (is.null(o$zph)) NULL else
        data.frame(endpoint = o$ep, population = o$pop, method = o$method, o$zph)))
      f_zph <- sv(transform(zph, p = signif(p, 4), chisq = round(chisq, 3)), paste0("sig_ph_tests_", sig, ".csv"))
      files <- c(files, f_zph)

      # forest: 메인 방법, 모든 종점 × 집단 × 모형의 score HR
      fr <- sm[sm$primary & sm$term == "score_sd", ]
      if (nrow(fr)) {
        ft <- data.frame(label = sprintf("%s, %s: %s (n = %d, events = %d)", fr$endpoint, fr$population, fr$model, fr$n, fr$events),
                         hr = fr$hr, lo = fr$lo, hi = fr$hi, p = fr$p, reference = FALSE,
                         shape = ifelse(fr$exploratory, "exploratory", "est"))
        files <- c(files, save_fig(
          forest_plot(ft, paste0(def$label, " score (", sig_method_labels[[pm]], ", per 1 SD) — ", g$label),
                      paste0("Multivariable models adjusted for ", cov_txt)),
          "GEO", paste0("sig_forest_", sig), width = 13, height = 1.6 + 0.35 * nrow(ft), subdir = sub, prefix = gse))
      }
      for (ep in names(g$endpoints)) {
        mp <- r$out[[paste(ep, "All", pm)]]$multi
        if (is.null(mp)) next
        files <- c(files, save_fig(
          forest_plot(mp$tab, paste0("Multivariable Cox: ", def$label, " (", sig_method_labels[[pm]], "), ", ep_lab[[ep]], " — ", g$label),
                      sprintf("n = %d, events = %d, EPV = %.1f", nrow(mp$m$d), mp$m$events, mp$m$epv)),
          "GEO", paste0("sig_forest_multi_", sig, "_", ep), width = 9, height = 1.6 + 0.35 * nrow(mp$tab), subdir = sub, prefix = gse))
      }

      # 분석 정보
      def_txt <- paste0("signatures$", sig, "=up:(", paste(def$up, collapse = ","), ");down:(", paste(def$down, collapse = ","), ")")
      ds_set <- function(k) {
        v <- unlist(g[[k]]); if (!is.null(names(v))) v <- paste0(names(v), ":", v)
        paste0("geo_datasets$", gse, "$", k, "=", paste(v, collapse = ","))
      }
      settings <- c(def_txt, "signature_primary_method", "signature_expr_filter", "signature_min_match", "signature_min_genes",
                    "geo_probe_rule", "group_levels", "cox_covariates", "epv_min", "stage_full", "stage_collapsed", "min_events",
                    "ph_split_months", "rmst_tau", "km_times", ds_set("cox_covariates_drop"), ds_set("extra_covariate"),
                    ds_set("subgroup"), paste0("bh_across_signatures=", do_bh))
      dat <- c(geo_data, paste0("signature source: ", def$source),
               if (is.null(r$tcga)) "TCGA direction reference: none" else paste0("TCGA direction reference: ", basename(r$tcga$file)))
      n_txt <- sprintf("%d tumors; genes used up %d/%d, down %d/%d; background %d genes", nrow(d), length(r$used$up), length(def$up),
                       length(r$used$down), length(def$down), nrow(flt$expr))
      models <- unlist(lapply(names(r$out), function(k) if (length(r$out[[k]]$models)) paste0("[", k, "] ", r$out[[k]]$models)))
      surv_model <- c("KM + log-rank: Surv(<EP>_months, <EP>_event) ~ group (median split over all tumors)",
                      "score_sd = (score - mean) / SD within each endpoint x population", models,
                      paste0("RMST: survRM2::rmst2 (", group_contrast_label(), ")"), "cox.zph for every model")
      tp <- function(name) file.path(out_dir("GEO", "tables", sub), out_file(gse, name))
      info <- list(
        list(Output = f_genes, Analysis = "Signature gene matching (platform) and per-gene QC (no outcome used)",
             Model = c("all-gene probe collapse by geo_probe_rule; status matched / not_on_platform / low_expression",
                       "psych::alpha(check.keys = FALSE) on z (down genes x -1)"),
             Settings = settings, n = n_txt, Data = dat, Packages = c("psych", "GEOquery", "Biobase")),
        list(Output = f_qc, Analysis = "Signature QC per set: matching, Cronbach's alpha, PC1 variance vs random gene sets",
             Model = "random gene sets of the same size from the filtered background", Settings = c(settings, "signature_random_n", "seed"),
             n = n_txt, Data = dat, Packages = c("psych", "stats")),
        list(Output = f_cor, Analysis = "Spearman correlation between the four signature scores (computed within this cohort)",
             Model = "cor(method = \"spearman\")", Settings = c(settings, "signature_ssgsea_normalize"), n = n_txt, Data = dat,
             Packages = c("GSVA", "singscore")),
        list(Output = f_raw, Analysis = "Signature external validation summary per method x endpoint x population x model",
             Model = surv_model, Settings = settings, n = n_txt, Data = dat, Packages = c("survival", "survRM2")),
        list(Output = tp(paste0("sig_survival_summary_", sig, ".csv")), Analysis = "Formatted signature external validation summary",
             Model = surv_model, Settings = settings, n = n_txt, Data = dat, Packages = c("survival", "survRM2")),
        list(Output = tp(paste0("sig_km_summary_", sig, ".csv")), Analysis = "Kaplan-Meier estimates by score group",
             Model = "survfit(Surv(<EP>_months, <EP>_event) ~ group)", Settings = c("group_levels", "km_times", ds_set("subgroup")),
             n = n_txt, Data = dat, Packages = "survival"),
        list(Output = f_zph, Analysis = "Proportional hazards test (Schoenfeld residuals)", Model = "cox.zph",
             Settings = "cox_covariates", n = n_txt, Data = dat, Packages = "survival"))
      cox_f <- files[grepl(paste0("_sig_cox_", sig, ".*\\.csv$"), files)]
      for (f in cox_f) info[[length(info) + 1]] <- list(Output = f, Analysis = "Cox regression, all models, all covariate rows",
                                                         Model = models, Settings = settings, n = n_txt, Data = dat, Packages = "survival")
      files <- c(files, write_analysis_info(info, "GEO", "13", subdir = sub, prefix = gse))
      move_stale_outputs("GEO", c("sig_"), files, subdir = sub, prefix = gse)
      created <- c(created, files)
      console[[paste(gse, sig)]] <- list(sm = sm, qc = r$qc$summary)
    }, error = function(e) {
      failures <<- c(failures, paste0(gse, " ", sig, " (출력): ", conditionMessage(e)))
      cat("  !!", sig, "출력 실패:", conditionMessage(e), "\n")
    })
  }
}

cat("\n생성된 파일 (", length(created), "개):\n", paste0("  ", created, collapse = "\n"), "\n", sep = "")
cat("\n===== 요약 (메인 방법 ", pm, ", score 항) =====\n", sep = "")
for (k in names(console)) {
  s <- console[[k]]$sm
  q <- console[[k]]$qc
  cat(k, " — 사용 유전자 ", paste(sprintf("%s %d/%d", q$set, q$n_used, q$n_defined), collapse = ", "), "\n", sep = "")
  for (i in which(s$primary & s$term == "score_sd")) {
    cat(sprintf("  %-4s %-5s %-44s n = %3d, 사건 %3d: HR %s, p %s%s%s\n", s$endpoint[i], s$population[i], s$model[i], s$n[i],
                s$events[i], fmt_hr(s$hr[i], s$lo[i], s$hi[i]), fmt_p(s$p[i]),
                if (is.na(s$direction_match[i])) "" else if (s$direction_match[i]) "  [TCGA와 같은 방향]" else "  [TCGA와 반대 방향]",
                if (s$exploratory[i]) " (탐색적)" else ""))
  }
}
if (length(failures)) {
  cat("\n실패한 분석 (", length(failures), "개):\n", paste0("  ", failures, collapse = "\n"), "\n", sep = "")
  quit(save = "no", status = 1)
}
cat("\n실패한 분석 없음\n")
