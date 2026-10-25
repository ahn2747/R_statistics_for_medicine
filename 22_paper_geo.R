# =============================================================
# 22_paper_geo.R [GSE ...]
# 9번 GEO 외부 검증 (cfg$paper$geo; 06/13의 cfg$geo_datasets와 별개)
#   1. series matrix (R/geo.R 캐시 공유: data/geo/<GSE>_<GPL>_eset.rds) + characteristics 파싱 (중복 key → "key #2")
#   2. fields 매핑 → keep_sample_type / exclude_values로 분석 코호트 → 시간 유도 → stage (stage_map 또는 T/N/M)
#      → recode_geo() (나이 그룹·성별·stage 그룹·종점 = TCGA와 같은 정의)
#   3. 점수: 분석 코호트 안에서 signature 유전자 log2 발현의 z 평균 (probe 축약 = geo_gene_matrix, cfg$geo_probe_rule),
#      코호트 median 분할 (High > median)
#   4. 종점별: KM log-rank + 그림, Cox 단변량 (그룹, 각 공변량), 다변량 (그룹 + cfg$paper$cox_vars − drop + add, 완전 사례)
# 결과 (파일명 앞에 <GSE>_): output/tables/GEO/paper/ paper_cox_<EP>, paper_pheno_QC, paper_scores;
#       output/figures/GEO/paper/ paper_km_<EP>, paper_forest_<EP>; GEO_paper_summary (전체), analysis_info_22
# =============================================================

source("config.R")
source("R/utils.R")
source("R/survival.R")
source("R/geo.R")
source("R/paper.R")

args <- commandArgs(trailingOnly = TRUE)
pc <- cfg$paper
paper_use_scales()
sets <- pc$geo[vapply(pc$geo, function(g) g$validates %in% cfg$cancers, TRUE)]
gses <- if (length(args)) args else names(sets)
unknown <- setdiff(gses, names(pc$geo))
if (length(unknown)) stop("cfg$paper$geo에 없는 데이터셋: ", paste(unknown, collapse = ", "))
genes <- cfg$signatures[[pc$signature]]$up
if (!length(genes)) stop("cfg$signatures$", pc$signature, "$up 없음", paper_where("signature"))
slab <- pc$scores[[1]]$label
sub <- pc$subdir
md5 <- function(f) unname(tools::md5sum(f))
save_fig_q <- function(...) suppressMessages(save_fig(...))

created <- character(); failures <- character(); summary_rows <- list(); info <- list()

for (gse in gses) {
  g <- pc$geo[[gse]]
  cat("\n==========", gse, "(", g$platform, ", TCGA-", g$validates, ") ==========\n", sep = "")
  r <- tryCatch({
    eset <- load_geo_eset(gse, g)
    ex <- prepare_expression(eset, gse)
    raw <- paper_parse_characteristics(Biobase::pData(eset), gse)
    print_raw_keys(raw)
    d <- paper_geo_map(raw, g, gse)
    qc <- list()
    add_qc <- function(item, value) qc[[length(qc) + 1]] <<- data.frame(item = item, value = as.character(value))
    st <- table(d$sample_type, useNA = "ifany")
    add_qc("samples (series matrix)", nrow(d))
    add_qc("sample types", paste(names(st), st, collapse = "; "))

    # ---- 분석 코호트 ----
    d <- d[d$sample_type %in% g$keep_sample_type, ]
    add_qc(paste0("kept sample type (", paste(g$keep_sample_type, collapse = ", "), ")"), nrow(d))
    for (k in names(g$exclude_values)) {
      hit <- d[[k]] %in% g$exclude_values[[k]]
      add_qc(paste0("excluded: ", k, " = ", paste(g$exclude_values[[k]], collapse = "/")), sum(hit))
      cat(sprintf("제외 (%s = %s): %d명\n", k, paste(g$exclude_values[[k]], collapse = "/"), sum(hit)))
      d <- d[!hit, ]
    }
    if (!nrow(d)) stop(gse, ": 분석 코호트 0명", paper_where(paste0("geo$", gse, "$keep_sample_type")))
    cat("분석 코호트:", nrow(d), "명\n")

    # ---- 시간 유도 (사건이면 if_event, 아니면 if_censored) ----
    for (k in names(g$derive_time)) {
      dv <- g$derive_time[[k]]
      ev <- d[[dv$event]]
      d[[k]] <- ifelse(is.na(ev), NA, ifelse(ev == dv$event_value, d[[dv$if_event]], d[[dv$if_censored]]))
      add_qc(paste0(k, " derived (event: ", dv$if_event, ", censored: ", dv$if_censored, ")"), sum(!is.na(d[[k]])))
    }

    # ---- stage ----
    if (!is.null(g$stage_map)) {
      bad <- unique(d$stage[!is.na(d$stage) & !d$stage %in% names(g$stage_map)])
      if (length(bad)) stop(gse, ": stage_map에 없는 값 ", paste(bad, collapse = ", "), paper_where(paste0("geo$", gse, "$stage_map")))
      d$stage <- unname(g$stage_map[d$stage])
    } else if (!is.null(g$stage_tnm)) {
      tn <- g$stage_tnm
      d$stage <- paper_stage_tnm(d[[tn$t]], d[[tn$n]], if (!is.null(tn$m)) d[[tn$m]], isTRUE(tn$assume_m0))
      add_qc("stage derived from T/N/M (AJCC 7th grouping)", paste(names(table(d$stage, useNA = "ifany")),
                                                                  table(d$stage, useNA = "ifany"), collapse = "; "))
    } else {
      stop(gse, ": stage_map 또는 stage_tnm 필요", paper_where(paste0("geo$", gse)))
    }
    g2 <- g
    g2$exclude_sample_type <- character()
    g2$fields <- g$fields[setdiff(names(g$fields), c("stage"))]
    d <- recode_geo(d, g2, gse)

    # ---- 점수 ----
    gm <- geo_gene_matrix(ex[, d$sample_id, drop = FALSE], Biobase::fData(eset), g, gse)
    miss <- setdiff(genes, rownames(gm$mat))
    if (length(miss)) stop(gse, ": 플랫폼에 없는 signature 유전자 ", paste(miss, collapse = ", "))
    z <- t(scale(t(gm$mat[genes, d$sample_id, drop = FALSE])))
    if (anyNA(z)) stop(gse, ": signature 유전자 발현 결측/분산 0")
    d$score <- colMeans(z)
    d$group <- median_split(d$score)
    pm <- gm$probe_map[gm$probe_map$symbol %in% genes, ]
    add_qc("signature probes (symbol:probe)", paste(pm$symbol, pm$probe_id, sep = ":", collapse = "; "))
    add_qc("score groups", paste(names(table(d$group)), table(d$group), collapse = "; "))
    for (k in c("age_g", "gender", "pathologic_stage", "pathologic_stage_12_34")) if (k %in% names(d)) {
      add_qc(k, paste(names(table(d[[k]], useNA = "ifany")), table(d[[k]], useNA = "ifany"), collapse = "; "))
    }
    p_sc <- save_table(d[, c("sample_id", "score", "group")], "GEO", "paper_scores.csv", sub, prefix = gse)

    cov_vars <- c(setdiff(pc$cox_vars, g$cox_covariates_drop), g$cox_covariates_add)
    cov_terms <- covariate_terms(cov_vars, d, paste0("paper$geo$", gse, "$cox_covariates_drop/add"), gse)
    labels <- c(term_labels(cov_terms), group = paste0(slab, " group"))
    files <- p_sc
    for (ep in names(g$endpoints)) {
      de <- geo_endpoint_data(d, ep, g2, gse, ep, keep = !is.na(d$group))
      de$ep_months <- de[[paste0(ep, "_months")]]; de$ep_event <- de[[paste0(ep, "_event")]]
      p_lr <- paper_logrank(de)
      one <- function(terms, model) {
        dm <- complete_cases(de, terms)
        fit <- fit_cox(dm, terms, time = "ep_months", event = "ep_event")
        data.frame(model = model, model_n = nrow(dm), model_events = sum(dm$ep_event),
                   epv = sum(dm$ep_event) / n_params(dm, terms), cox_rows(fit, dm, terms, labels, event = "ep_event"))
      }
      uni <- do.call(rbind, lapply(c("group", cov_terms), function(t) one(t, "Univariable")))
      mul <- one(c("group", cov_terms), "Multivariable")
      raw_c <- rbind(uni, mul)
      ep_lab <- g$endpoints[[ep]]$label %||% ep
      f1 <- save_table(raw_c, "GEO", paste0("paper_cox_", ep, "_raw.csv"), sub, prefix = gse)
      tab <- cbind(Model = ifelse(duplicated(raw_c$model), "", raw_c$model), format_cox(raw_c))
      f2 <- save_df_table(tab, "GEO", paste0("paper_cox_", ep), subdir = sub, prefix = gse, caption = paste0(
        ep_lab, ", ", g$label, ". ", slab, " = mean z-score of ", paste(genes, collapse = ", "),
        " (log2 expression, z within the cohort); High vs Low = cohort median split. Multivariable: ",
        adjust_text(cov_vars), sprintf(" (complete cases n = %d, events = %d). ", mul$model_n[1], mul$model_events[1]),
        if (length(g$cox_covariates_drop)) paste0("Not available in this dataset: ", paste(g$cox_covariates_drop, collapse = ", "), ".") else ""))
      f3 <- save_fig_q(paper_km_plot(de, paste0(g$label), paste(slab, "group"), paste(ep_lab, "probability"), p_lr),
                       "GEO", paste0("paper_km_", ep), width = 7, height = 7, subdir = sub, prefix = gse)
      f4 <- save_fig_q(forest_plot(mul, paste0(ep_lab, " — ", g$label), sprintf("Multivariable Cox; n = %d, events = %d",
                                                                                  mul$model_n[1], mul$model_events[1])),
                       "GEO", paste0("paper_forest_", ep), width = 10, height = 1.2 + 0.35 * nrow(mul), subdir = sub, prefix = gse)
      files <- c(files, f1, f2, f3, f4)
      gu <- uni[uni$term == "group" & !uni$reference, ]; gm_ <- mul[mul$term == "group" & !mul$reference, ]
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        dataset = gse, validates = g$validates, endpoint = ep, n = nrow(de), events = sum(de$ep_event),
        n_low = sum(de$group == group_ref()), n_high = sum(de$group == group_alt()), logrank_p = p_lr,
        uni_hr = gu$hr, uni_lo = gu$lo, uni_hi = gu$hi, uni_p = gu$p,
        multi_n = mul$model_n[1], multi_events = mul$model_events[1], multi_hr = gm_$hr, multi_lo = gm_$lo,
        multi_hi = gm_$hi, multi_p = gm_$p, covariates = paste(cov_terms, collapse = ","),
        covariates_dropped = paste(g$cox_covariates_drop, collapse = ","))
      cat(sprintf("  %s: n %d (사건 %d), log-rank p %s | 단변량 HR %s p %s | 다변량 (n %d) HR %s p %s\n", ep, nrow(de),
                  sum(de$ep_event), fmt_p(p_lr), fmt_hr(gu$hr, gu$lo, gu$hi), fmt_p(gu$p), mul$model_n[1],
                  fmt_hr(gm_$hr, gm_$lo, gm_$hi), fmt_p(gm_$p)))
    }
    files <- c(files, save_table(do.call(rbind, qc), "GEO", "paper_pheno_QC.csv", sub, prefix = gse))
    info[[length(info) + 1]] <- list(
      Output = files[length(files)], Analysis = paste("GEO validation,", gse),
      Model = list("survdiff(Surv(<EP>_months, <EP>_event) ~ group)", paste0("coxph(... ~ group + ", paste(cov_terms, collapse = " + "), ")")),
      Settings = c("geo_probe_rule", "geo_na_values", paste0("genes=", paste(genes, collapse = ",")),
                   paste0("cox_drop=", paste(g$cox_covariates_drop, collapse = ",")), paste0("cox_add=", paste(g$cox_covariates_add, collapse = ","))),
      n = nrow(d), Data = sprintf("%s / %s series matrix (GEOquery)", gse, g$platform), Packages = c("GEOquery", "survival"))
    created <- c(created, files)
    TRUE
  }, error = function(e) {
    failures <<- c(failures, paste0(gse, ": ", conditionMessage(e)))
    cat("  !!", gse, "실패:", conditionMessage(e), "\n")
    NULL
  })
}

if (length(summary_rows)) {
  sm <- do.call(rbind, summary_rows)
  created <- c(created, save_table(sm, "GEO", "paper_summary_raw.csv", sub, prefix = "GEO"))
  tab <- data.frame(Dataset = sm$dataset, `TCGA cohort` = sm$validates, Endpoint = sm$endpoint, N = sm$n, Events = sm$events,
                    `Log-rank p` = fmt_p(sm$logrank_p), `Univariable HR (95% CI)` = fmt_hr(sm$uni_hr, sm$uni_lo, sm$uni_hi),
                    `Univariable p` = fmt_p(sm$uni_p), `Multivariable n (events)` = paste0(sm$multi_n, " (", sm$multi_events, ")"),
                    `Multivariable HR (95% CI)` = fmt_hr(sm$multi_hr, sm$multi_lo, sm$multi_hi), `Multivariable p` = fmt_p(sm$multi_p),
                    check.names = FALSE)
  created <- c(created, save_df_table(tab, "GEO", "paper_summary", subdir = sub, prefix = "GEO", landscape = TRUE, font_size = 8,
    caption = paste0("External validation of the ", slab, " in GEO cohorts (High vs Low, cohort median split). ",
                     "Multivariable models use the TCGA covariates available in each dataset.")))
}
if (length(info)) created <- c(created, write_analysis_info(info, "GEO", "22", sub, prefix = "GEO"))

cat("\n생성된 파일", length(created), "개\n")
if (length(failures)) {
  cat("\n실패:\n", paste0("  ", failures, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("실패한 분석 없음\n")
