# =============================================================
# 12_sig_survival.R [SIG ...] [COAD READ ...]
# signature 점수의 전체 생존(OS) 분석 — 암종마다 따로, 점수 방법 4개 모두 (메인 = cfg$signature_primary_method)
#   입력: data/processed/<C>_sigscores_<SIG>.csv (11), <C>_clinical_all.rds, <C>_clinical.rds (01)
#   High/Low = 방법별 median_split(), 점수 샘플 전체 (신보조요법 제외 전; .sav 그룹·06과 같은 방식)
#   생존 코호트 = _clinical.rds (cfg$exclude_neoadjuvant 적용) ∩ 점수, os_months/status 있음
#   1 SD = 생존 코호트 안의 SD (score_sd); 연속형 Cox는 1 SD당 HR (cfg$cox_gene_term과 무관)
#   모형: KM/log-rank/RMST (High vs Low), 단변량 Cox (High vs Low, per 1 SD), 다변량 (per 1 SD + cfg$cox_covariates,
#         EPV 규칙), 다변량 + strata(cfg$strata_var), cox.zph, score 항 PH 위반 시 시간 분할 Cox
#   다중비교: 기존 유전자 BH와 섞지 않음. signature가 여러 개이고 전체 실행일 때만 메인 방법 행에서 signature 간 BH
# 결과 (파일명 앞에 <C>_): output/tables/<C>/signature/<SIG>/ sig_survival_summary_raw_<SIG>, sig_survival_summary_<SIG>,
#       sig_cox_<SIG>, sig_km_summary_<SIG>, sig_ph_tests_<SIG>, analysis_info_12
#       output/figures/<C>/signature/<SIG>/ sig_km_<SIG>_<method>, sig_forest_<SIG>, sig_forest_multi_<SIG>
# =============================================================

source("config.R")
source("R/utils.R")
source("R/survival.R")
source("R/signature.R")

args <- commandArgs(trailingOnly = TRUE)
cancers <- if (any(toupper(args) %in% cfg$cancers)) toupper(args[toupper(args) %in% cfg$cancers]) else cfg$cancers
sigs <- args[!toupper(args) %in% cfg$cancers]
all_sigs <- !length(sigs)
if (all_sigs) sigs <- names(cfg$signatures)
if (!length(sigs)) stop("cfg$signatures가 비어 있음 → config.R에 signature 정의 추가")
unknown <- setdiff(sigs, names(cfg$signatures))
if (length(unknown)) stop("cfg$signatures에 없는 signature: ", paste(unknown, collapse = ", "))
pm <- cfg$signature_primary_method
if (length(pm) != 1 || is.na(pm) || !pm %in% sig_methods) {
  stop("cfg$signature_primary_method는 ", paste(sig_methods, collapse = ", "), " 중 하나 (결과를 보기 전에 지정)")
}
defs <- lapply(setNames(sigs, sigs), sig_def)
cat("signature:", paste(sigs, collapse = ", "), "/ 암종:", paste(cancers, collapse = ", "), "/ 메인 방법:", pm, "\n")

md5 <- function(f) unname(tools::md5sum(f))
fmt_rmst <- function(x, lo, hi) ifelse(is.na(x), NA, sprintf("%.2f (%.2f, %.2f)", x, lo, hi))

created  <- character()
failures <- character()
console  <- list()

for (cancer in cancers) {
  cat("\n==========", cancer, "==========\n")
  f_all <- processed_path(cancer, "clinical_all.rds")
  f_cl  <- processed_path(cancer, "clinical.rds")
  for (f in c(f_all, f_cl)) if (!file.exists(f)) stop(f, " 없음 → 01_import_clinical.R 먼저 실행")
  clin_all <- readRDS(f_all)
  clin     <- readRDS(f_cl)
  # _clinical.rds가 현재 cfg$exclude_neoadjuvant로 만들어졌는지 (01 이후 설정이 바뀌면 중단)
  if (isTRUE(cfg$exclude_neoadjuvant) && any(clin$neoadjuvant_flag)) {
    stop(cancer, ": cfg$exclude_neoadjuvant = TRUE인데 ", basename(f_cl), "에 신보조요법 환자 있음 → 01_import_clinical.R 재실행")
  }
  if (!isTRUE(cfg$exclude_neoadjuvant) && nrow(clin) != nrow(clin_all)) {
    stop(cancer, ": cfg$exclude_neoadjuvant = FALSE인데 ", basename(f_cl), "(", nrow(clin), "명) ≠ ", basename(f_all),
         "(", nrow(clin_all), "명) → 01_import_clinical.R 재실행")
  }
  sav_txt <- sprintf("%s (md5 %s)", basename(sav_path(cancer)), md5(sav_path(cancer)))

  res <- list()
  for (sig in sigs) {
    def <- defs[[sig]]
    sub <- sig_subdir(sig)
    cat("\n--", sig, "--\n")
    r <- tryCatch({
      sf <- sig_scores_path(cancer, sig)
      if (!file.exists(sf)) stop(sf, " 없음 → Rscript 11_sig_scores.R ", sig, " ", cancer)
      sc <- read.csv(sf, stringsAsFactors = FALSE)
      if (!identical(names(sc), c("sample_id", sig_methods))) {
        stop(sf, ": 열이 sample_id, ", paste(sig_methods, collapse = ", "), "이 아님 → 11 재실행")
      }
      stray <- setdiff(sc$sample_id, clin_all$sample_id)
      if (length(stray)) stop(sf, ": 임상(전체)에 없는 환자 ", paste(head(stray, 10), collapse = ", "), " → 11 재실행")
      if (file.mtime(sf) < file.mtime(f_all)) {
        warning(cancer, " ", sig, ": 점수 파일이 ", basename(f_all), "보다 오래됨 → 11 재실행 권장", call. = FALSE)
      }
      # 방법별 High/Low: 점수 샘플 전체에서 median split (제외 전)
      for (m in sig_methods) sc[[paste0(m, "_group")]] <- median_split(sc[[m]])
      d0 <- dplyr::inner_join(clin, sc, by = "sample_id")
      base <- add_scaled_terms(d0[!is.na(d0$os_months) & !is.na(d0$status), ])
      cat(sprintf("  점수 %d명 → 분석용 임상과 결합 %d명 → 생존 정보 있음 %d명 (사건 %d)\n",
                  nrow(sc), nrow(d0), nrow(base), sum(base$status)))
      cov_terms <- covariate_terms(cfg$cox_covariates, base, "cox_covariates", paste(cancer, sig))

      per <- list()
      for (m in sig_methods) {
        d <- base
        sd_raw <- stats::sd(d[[m]])
        d$score_sd <- (d[[m]] - mean(d[[m]])) / sd_raw
        d$group <- d[[paste0(m, "_group")]]
        cat(sprintf("  [%s] %s %d / %s %d\n", m, group_ref(), sum(d$group == group_ref()), group_alt(), sum(d$group == group_alt())))
        ss <- sig_survival_set(d, paste(def$label, sig_method_labels[[m]]), cov_terms, strata_var = cfg$strata_var,
                               where = paste(cancer, sig, m))
        ss$score_sd_raw <- sd_raw
        ss$data <- d
        per[[m]] <- ss
      }
      list(def = def, per = per, n_scores = nrow(sc), n_join = nrow(d0), base_n = nrow(base), base_ev = sum(base$status),
           cov_terms = cov_terms, sf = sf)
    }, error = function(e) {
      failures <<- c(failures, paste0(cancer, " ", sig, ": ", conditionMessage(e)))
      cat("  !!", sig, "실패:", conditionMessage(e), "\n")
      NULL
    })
    if (!is.null(r)) res[[sig]] <- r
  }

  # ---- signature 간 BH (메인 방법만, 전체 실행이고 signature가 2개 이상일 때) ----
  raw_all <- list()
  for (sig in names(res)) {
    r <- res[[sig]]
    raw_all[[sig]] <- dplyr::bind_rows(lapply(sig_methods, function(m) {
      ss <- r$per[[m]]
      if (is.null(ss$rows)) return(data.frame(method = m, note = paste("skipped:", ss$skipped)))
      data.frame(cohort = cancer, sig = sig, method = m, primary = m == pm, endpoint = "OS", population = "All",
                 ss$rows, logrank_p = ss$p_lr, ss$rmst$est, n_low = ss$n_low, n_high = ss$n_high,
                 score_sd_raw = ss$score_sd_raw, note = ss$rmst$note, stringsAsFactors = FALSE)
    }))
  }
  do_bh <- all_sigs && length(raw_all) > 1
  if (do_bh) {
    prim <- dplyr::bind_rows(lapply(raw_all, function(x) x[x$primary %in% TRUE, c("sig", "model", "p", "logrank_p")]))
    for (mdl in unique(prim$model)) {
      i <- prim$model == mdl
      prim$q[i] <- p.adjust(prim$p[i], "BH")
    }
    lr <- unique(prim[, c("sig", "logrank_p")]); lr$logrank_q <- p.adjust(lr$logrank_p, "BH")
  }
  for (sig in names(raw_all)) {
    x <- raw_all[[sig]]
    x$q_across_signatures <- NA_real_; x$logrank_q_across_signatures <- NA_real_
    if (do_bh) {
      i <- which(x$primary %in% TRUE)
      x$q_across_signatures[i] <- prim$q[match(paste(sig, x$model[i]), paste(prim$sig, prim$model))]
      x$logrank_q_across_signatures[i] <- lr$logrank_q[match(sig, lr$sig)]
    }
    raw_all[[sig]] <- x
  }

  # ---- 출력 (signature마다) ----
  for (sig in names(res)) {
    r <- res[[sig]]; def <- r$def; sub <- sig_subdir(sig); raw <- raw_all[[sig]]
    files <- character()
    tryCatch({
      cov_txt <- adjust_text()
      f_raw <- save_table(raw, cancer, paste0("sig_survival_summary_raw_", sig, ".csv"), sub)
      files <- c(files, f_raw)

      ok <- !is.na(raw$model %||% NA)
      sm <- raw[ok, ]
      first <- !duplicated(sm$method)
      tab <- data.frame(
        Method = ifelse(first, unname(sig_method_labels[sm$method]), ""),
        Primary = ifelse(first, ifelse(sm$primary, "Yes", "No"), ""),
        Model = sm$model, N = sm$n, Events = sm$events,
        `HR (95% CI)` = fmt_hr(sm$hr, sm$lo, sm$hi), p = fmt_p(sm$p),
        check.names = FALSE)
      if (do_bh) tab$`q (across signatures, primary method)` <- ifelse(sm$primary, fmt_p(sm$q_across_signatures), "")
      tab$EPV <- sprintf("%.1f", sm$epv)
      tab$`Stage variable` <- sm$stage_levels
      tab$`PH p (score term)` <- fmt_p(sm$ph_p_term)
      tab$`Log-rank p` <- ifelse(first, fmt_p(sm$logrank_p), "")
      tab[[paste0("RMST difference at ", cfg$rmst_tau, " mo, months (95% CI)")]] <-
        ifelse(first, fmt_rmst(sm$rmst_diff, sm$rmst_lo, sm$rmst_hi), "")
      tab$`RMST p` <- ifelse(first, fmt_p(sm$rmst_p), "")
      tab[[paste0("Exploratory (EPV < ", cfg$epv_min, ")")]] <- sm$exploratory
      ord <- order(sm$method != pm, match(sm$method, sig_methods))   # 메인 방법 먼저
      tab <- tab[ord, ]
      tab$Method <- ifelse(!duplicated(sm$method[ord]), unname(sig_method_labels[sm$method[ord]]), "")
      tab$Primary <- ifelse(!duplicated(sm$method[ord]), ifelse(sm$primary[ord], "Yes", "No"), "")
      files <- c(files, save_df_table(
        tab, cancer, paste0("sig_survival_summary_", sig), landscape = TRUE, font_size = 7, subdir = sub,
        caption = paste0(
          "Overall survival by ", def$label, " signature score, TCGA-", cancer, " (", def$source, "). ",
          "Primary method: ", sig_method_labels[[pm]], " (prespecified). Continuous HRs are per 1 SD of the score within the survival cohort; ",
          "High vs Low = median split of the score over all scored patients. Multivariable models adjusted for ", cov_txt,
          "; sensitivity stratified by ", cfg$strata_var, ". Time-split HRs (0–", cfg$ph_split_months, " / >",
          cfg$ph_split_months, " months) only when the score term violated proportional hazards. RMST difference = ",
          group_alt(), " − ", group_ref(), " up to ", cfg$rmst_tau, " months. Not adjusted together with the gene-level analyses",
          if (do_bh) "; q = Benjamini–Hochberg across signatures for the primary method." else ".")))

      # Cox 표 (모든 공변량 행)
      cx <- dplyr::bind_rows(lapply(sig_methods, function(m) {
        x <- r$per[[m]]$cox
        if (is.null(x)) NULL else data.frame(Method = sig_method_labels[[m]], x, check.names = FALSE)
      }))
      files <- c(files, save_df_table(cx, cancer, paste0("sig_cox_", sig), landscape = TRUE, font_size = 7, subdir = sub,
                                      caption = paste0("Cox regression for overall survival: ", def$label,
                                                       " signature score, TCGA-", cancer, " (all models, all covariate rows)")))

      # KM 요약
      km <- dplyr::bind_rows(lapply(sig_methods, function(m) {
        ss <- r$per[[m]]
        if (is.null(ss$km)) NULL else data.frame(method = m, ss$km, logrank_p = ss$p_lr)
      }))
      km_out <- data.frame(Method = ifelse(duplicated(km$method), "", unname(sig_method_labels[km$method])),
                           Group = km$group, N = km$n, Events = km$events,
                           `Median OS, months (95% CI)` = fmt_median(km$median, km$median_lo, km$median_hi),
                           check.names = FALSE)
      for (tm in cfg$km_times) {
        col <- paste0("surv_", tm / 12, "y")
        km_out[[paste0(tm / 12, "-year OS (95% CI)")]] <- fmt_surv(km[[col]], km[[paste0(col, "_lo")]], km[[paste0(col, "_hi")]])
      }
      km_out$`Log-rank p` <- ifelse(duplicated(km$method), "", fmt_p(km$logrank_p))
      files <- c(files, save_df_table(km_out, cancer, paste0("sig_km_summary_", sig), subdir = sub,
                                      caption = paste0("Kaplan–Meier estimates by ", def$label, " signature score (",
                                                       group_contrast_label(), ", median split), TCGA-", cancer, ". NR = not reached.")))
      zph <- dplyr::bind_rows(lapply(sig_methods, function(m) {
        z <- r$per[[m]]$zph
        if (is.null(z)) NULL else data.frame(method = m, z)
      }))
      f_zph <- save_table(transform(zph, p = signif(p, 4), chisq = round(chisq, 3)), cancer,
                          paste0("sig_ph_tests_", sig, ".csv"), sub)
      files <- c(files, f_zph)

      # 그림: 방법별 KM, 방법 × 모형 forest, 메인 방법 다변량 forest
      for (m in sig_methods) {
        ss <- r$per[[m]]
        if (is.null(ss$fit_km)) next
        p <- km_plot(ss$fit_km, ss$data, paste(def$label, sig_method_labels[[m]]), cancer, ss$p_lr,
                     title = paste0(def$label, " score (", sig_method_labels[[m]], ") — TCGA-", cancer))
        files <- c(files, suppressMessages(save_fig(p, cancer, paste0("sig_km_", sig, "_", m), width = 7, height = 6, subdir = sub)))
      }
      fr <- sm[sm$term == "score_sd", ]
      fr <- fr[order(fr$method != pm, match(fr$method, sig_methods)), ]
      if (nrow(fr)) {
        ft <- data.frame(label = sprintf("%s%s: %s (n = %d, events = %d)", sig_method_labels[fr$method],
                                         ifelse(fr$primary, " [primary]", ""), fr$model, fr$n, fr$events),
                         hr = fr$hr, lo = fr$lo, hi = fr$hi, p = fr$p, reference = FALSE,
                         shape = ifelse(fr$exploratory, "exploratory", "est"))
        files <- c(files, save_fig(
          forest_plot(ft, paste0(def$label, " signature score (per 1 SD) — TCGA-", cancer),
                      paste0("Multivariable models adjusted for ", cov_txt)),
          cancer, paste0("sig_forest_", sig), width = 12, height = 1.6 + 0.35 * nrow(ft), subdir = sub))
      }
      mp <- r$per[[pm]]$multi
      if (!is.null(mp)) {
        files <- c(files, save_fig(
          forest_plot(mp$tab, paste0("Multivariable Cox: ", def$label, " (", sig_method_labels[[pm]], ") — TCGA-", cancer),
                      sprintf("n = %d, events = %d, EPV = %.1f%s", nrow(mp$m$d), mp$m$events, mp$m$epv,
                              if (mp$m$exploratory) " (exploratory)" else "")),
          cancer, paste0("sig_forest_multi_", sig), width = 9, height = 1.6 + 0.35 * nrow(mp$tab), subdir = sub))
      }

      # 분석 정보
      models <- unlist(lapply(sig_methods, function(m) if (length(r$per[[m]]$models)) paste0("[", m, "] ", r$per[[m]]$models)))
      settings <- c(paste0("signatures$", sig, "=up:(", paste(def$up, collapse = ","), ");down:(", paste(def$down, collapse = ","), ")"),
                    "signature_primary_method", "group_levels", "cox_covariates", "epv_min", "stage_full", "stage_collapsed",
                    "min_events", "strata_var", "collapse_small_levels", "ph_split_months", "rmst_tau", "km_times",
                    "exclude_neoadjuvant", paste0("bh_across_signatures=", do_bh))
      dat <- c(sav_txt, sprintf("%s (md5 %s)", basename(r$sf), md5(r$sf)), sprintf("%s (md5 %s)", basename(f_cl), md5(f_cl)),
               paste0("signature source: ", def$source))
      n_txt <- sprintf("scores n=%d; joined with analytic clinical n=%d; survival cohort n=%d, events=%d",
                       r$n_scores, r$n_join, r$base_n, r$base_ev)
      surv_model <- c("KM + log-rank: Surv(os_months, status) ~ group (median split over all scored patients)",
                      "score_sd = (score - mean) / SD within the survival cohort", models,
                      paste0("RMST: survRM2::rmst2 (", group_contrast_label(), ")"), "cox.zph for every unstratified model")
      tp <- function(name) file.path(out_dir(cancer, "tables", sub), out_file(cancer, name))
      info <- list(
        list(Output = f_raw, Analysis = "Signature score survival summary per method x model", Model = surv_model,
             Settings = settings, n = n_txt, Data = dat, Packages = c("survival", "survRM2")),
        list(Output = tp(paste0("sig_survival_summary_", sig, ".csv")), Analysis = "Formatted signature survival summary",
             Model = surv_model, Settings = settings, n = n_txt, Data = dat, Packages = c("survival", "survRM2")),
        list(Output = tp(paste0("sig_cox_", sig, ".csv")), Analysis = "Cox regression, all models, all covariate rows",
             Model = models, Settings = settings, n = n_txt, Data = dat, Packages = "survival"),
        list(Output = tp(paste0("sig_km_summary_", sig, ".csv")), Analysis = "Kaplan-Meier estimates by score group",
             Model = "survfit(Surv(os_months, status) ~ group)", Settings = c("group_levels", "km_times"),
             n = n_txt, Data = dat, Packages = "survival"),
        list(Output = f_zph, Analysis = "Proportional hazards test (Schoenfeld residuals)", Model = "cox.zph",
             Settings = c("cox_covariates"), n = n_txt, Data = dat, Packages = "survival"))
      files <- c(files, write_analysis_info(info, cancer, "12", subdir = sub))
      move_stale_outputs(cancer, c("sig_survival_summary", "sig_cox_", "sig_km_", "sig_ph_tests_", "sig_forest",
                                   "analysis_info_12"), files, subdir = sub)
      created <<- c(created, files)
      console[[paste(cancer, sig)]] <- sm
    }, error = function(e) {
      failures <<- c(failures, paste0(cancer, " ", sig, " (출력): ", conditionMessage(e)))
      cat("  !!", sig, "출력 실패:", conditionMessage(e), "\n")
    })
  }
}

cat("\n생성된 파일 (", length(created), "개):\n", paste0("  ", created, collapse = "\n"), "\n", sep = "")
cat("\n===== 요약 (score 항) =====\n")
for (k in names(console)) {
  s <- console[[k]]
  cat(k, "\n")
  for (i in which(s$term == "score_sd" | s$model == "Univariable (High vs Low)")) {
    cat(sprintf("  %-9s%s %-44s n = %3d, 사건 %3d: HR %s, p %s%s\n", s$method[i], if (s$primary[i]) "*" else " ",
                s$model[i], s$n[i], s$events[i], fmt_hr(s$hr[i], s$lo[i], s$hi[i]), fmt_p(s$p[i]),
                if (s$exploratory[i]) " (탐색적)" else ""))
  }
  lr <- unique(s[, c("method", "logrank_p")])
  cat("  log-rank:", paste(sprintf("%s %s", lr$method, fmt_p(lr$logrank_p)), collapse = ", "), "\n")
}
if (length(failures)) {
  cat("\n실패한 분석 (", length(failures), "개):\n", paste0("  ", failures, collapse = "\n"), "\n", sep = "")
  quit(save = "no", status = 1)
}
cat("\n실패한 분석 없음\n")
