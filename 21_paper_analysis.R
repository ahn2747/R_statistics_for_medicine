# =============================================================
# 21_paper_analysis.R [C ...]
# 논문 분석 (암종마다 따로; 입력 = 20의 data/processed/<C>_paper.rds, 설정 = cfg$paper)
#   점수별 (cfg$paper$scores, 첫 번째 = 메인) → output/{tables,figures}/<C>/paper/<점수>/
#     1. Pearson 상관 (점수, pack-years, 나이, 유전자 log2 발현)          cor_pearson
#     2. Shapiro–Wilk + Spearman 상관 (같은 변수, pairwise n)             shapiro, cor_spearman
#     4. Table 1: 점수 High vs Low × 범주형 (χ² / Fisher, 결측 n, 기대빈도<5)  table1_signature, table1_tests
#     5. Cox (종점별): 단변량, 다변량 (완전 사례), 다변량 − sens_drop (민감도) + forest, KM
#     8. purity: 점수 ~ purity 상관, 부분 상관 (purity 보정), 다변량 + purity
#    10. cut-off 민감도: 연속형 (1 SD당), 삼분위/사분위 KM + Cox
#   점수와 무관 → output/{tables,figures}/<C>/paper/
#     3. GSEA 원 결과 재출력 (cfg$paper$gsea_file, 재계산 없음)              gsea_table6
#     6. 후보 유전자 KM log-rank (merged.rds의 median 그룹), BH q는 보고용     km_candidates_<EP>
#     7. 개별 유전자 KM                                                     km_<GENE>_<EP>, km_genes_summary
# 다중비교 보정 없음 (6번 BH q 제외). 점수 그룹 = 20에서 만든 group_<점수>
# =============================================================

source("config.R")
source("R/utils.R")
source("R/survival.R")
source("R/paper.R")

args <- commandArgs(trailingOnly = TRUE)
cancers <- if (any(toupper(args) %in% cfg$cancers)) toupper(args[toupper(args) %in% cfg$cancers]) else cfg$cancers
pc <- cfg$paper
paper_use_scales()
set.seed(cfg$seed)

created  <- character()
failures <- character()
run <- function(what, expr) {
  tryCatch(expr, error = function(e) {
    failures <<- c(failures, paste0(what, ": ", conditionMessage(e)))
    cat("  !!", what, "실패:", conditionMessage(e), "\n")
    NULL
  })
}
md5 <- function(f) unname(tools::md5sum(f))
save_fig_q <- function(...) suppressMessages(save_fig(...))   # survminer "Ignoring unknown labels" 메시지 억제
ep_label <- function(ep) pc$endpoints[[ep]]$label
ep_ylab  <- function(ep) paste(ep_label(ep), "probability")

# 모형 결과 → raw 행 (모형 이름, n, 사건, EPV 포함)
cox_fit_rows <- function(de, terms, labels, model) {
  de <- complete_cases(de, terms)
  fit <- fit_cox(de, terms, time = "ep_months", event = "ep_event")
  r <- cox_rows(fit, de, terms, labels, event = "ep_event")
  data.frame(model = model, model_n = nrow(de), model_events = sum(de$ep_event),
             epv = sum(de$ep_event) / n_params(de, terms), r, check.names = FALSE)
}
fmt_model_tab <- function(raw) {
  out <- format_cox(raw)
  cbind(Model = ifelse(duplicated(raw$model), "", raw$model), out)
}

info <- list()
add_info <- function(path, analysis, model, settings, n, data, pkgs = c("survival")) {
  info[[length(info) + 1]] <<- list(Output = path[length(path)], Analysis = analysis, Model = model,
                                    Settings = settings, n = n, Data = data, Packages = pkgs)
}

for (cancer in cancers) {
  cat("\n==========", cancer, "==========\n")
  f_in <- processed_path(cancer, "paper.rds")
  if (!file.exists(f_in)) stop(f_in, " 없음 → 20_paper_prep.R 먼저 실행")
  d0 <- readRDS(f_in)
  d0$purity_pct <- 100 * d0$purity
  d0 <- apply_reference_levels(add_scaled_terms(d0))
  info <- list()
  data_txt <- sprintf("%s (md5 %s)", basename(f_in), md5(f_in))
  top <- pc$subdir
  labels_all <- c(var_labels, smoking = "Smoking status")

  cov_terms <- covariate_terms(pc$cox_vars, d0, "paper$cox_vars", cancer)
  sens_terms <- setdiff(cov_terms, unname(sapply(pc$sens_drop, model_term)))

  # ---- 3. GSEA 원 결과 재출력 ----
  gf <- pc$gsea_file[[cancer]]
  if (!is.null(gf)) run(paste(cancer, "GSEA"), {
    if (!file.exists(gf)) stop(gf, " 없음", paper_where("gsea_file"))
    g <- read.csv(gf, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
    need <- c("Term", "ES", "NES", "NOM p-val", "FDR q-val")
    if (length(setdiff(need, names(g)))) stop(gf, ": 열 없음 ", paste(setdiff(need, names(g)), collapse = ", "))
    g <- g[order(g[["FDR q-val"]], -abs(g$NES)), ]
    tab <- data.frame(Pathway = g$Term, ES = sprintf("%.3f", g$ES), NES = sprintf("%.3f", g$NES),
                      `Nominal p` = fmt_p(g[["NOM p-val"]]), `FDR q` = fmt_p(g[["FDR q-val"]]),
                      check.names = FALSE)
    if ("Tag %" %in% names(g)) tab$`Leading edge (tag %)` <- g[["Tag %"]]
    p <- save_df_table(tab, cancer, "gsea_table6", subdir = top, font_size = 8, caption = paste0(
      "GSEA of MSigDB Hallmark gene sets, High vs Low signature group, TCGA-", cancer,
      ". Values re-printed from the original result file (", basename(gf), "); not recomputed."))
    created <<- c(created, p)
    add_info(p, "GSEA re-print (Table 6)", "re-print of gseapy result file", paste0("gsea_file=", basename(gf)),
             nrow(g), sprintf("%s (md5 %s)", basename(gf), md5(gf)), character())
    cat("  3. GSEA:", nrow(g), "경로 (FDR q < 0.05:", sum(g[["FDR q-val"]] < 0.05), ")\n")
  })

  # ---- 6. 후보 유전자 KM log-rank ----
  for (ep in pc$candidate_endpoints) run(paste(cancer, "candidates", ep), {
    de0 <- paper_endpoint_data(d0, ep, cancer)
    miss <- pc$candidate_genes[!paste0(gene_key(pc$candidate_genes), "_group") %in% names(de0)]
    if (length(miss)) stop("후보 유전자 그룹 열 없음: ", paste(miss, collapse = ", "),
                           " → csv_download.py + 02b, 또는", paper_where("candidate_genes"))
    res <- do.call(rbind, lapply(pc$candidate_genes, function(gn) {
      de <- de0; de$group <- de[[paste0(gene_key(gn), "_group")]]
      de <- de[!is.na(de$group), ]
      data.frame(gene = gn, n = nrow(de), events = sum(de$ep_event),
                 n_low = sum(de$group == group_ref()), n_high = sum(de$group == group_alt()),
                 logrank_p = paper_logrank(de))
    }))
    res$bh_q_report_only <- p.adjust(res$logrank_p, "BH")
    res$pass_p_lt_alpha <- res$logrank_p < pc$candidate_alpha
    p1 <- save_table(res, cancer, paste0("km_candidates_", ep, "_raw.csv"), top)
    tab <- data.frame(Gene = res$gene, N = res$n, Events = res$events,
                      `Low / High` = paste(res$n_low, res$n_high, sep = " / "),
                      `Log-rank p` = fmt_p(res$logrank_p), `BH q (reporting only)` = fmt_p(res$bh_q_report_only),
                      `p < 0.05` = ifelse(res$pass_p_lt_alpha, "Yes", "No"), check.names = FALSE)
    p2 <- save_df_table(tab, cancer, paste0("km_candidates_", ep), subdir = top, caption = paste0(
      ep_label(ep), " by candidate gene expression (High vs Low, median split), TCGA-", cancer,
      ". Log-rank test; selection criterion p < ", pc$candidate_alpha,
      ". Benjamini–Hochberg q across the ", nrow(res), " genes is shown for reporting only."))
    created <<- c(created, p1, p2)
    add_info(c(p1, p2), "Candidate gene screening (KM log-rank)", "survdiff(Surv(time, event) ~ gene group)",
             c(paste0("candidate_genes=", paste(pc$candidate_genes, collapse = ",")), paste0("alpha=", pc$candidate_alpha),
               "BH=report only"), nrow(de0), data_txt)
    cat(sprintf("  6. 후보 %s: p < %.2f %d/%d (%s)\n", ep, pc$candidate_alpha, sum(res$pass_p_lt_alpha), nrow(res),
                paste(res$gene[res$pass_p_lt_alpha], collapse = ", ")))
  })

  # ---- 7. 개별 유전자 KM ----
  run(paste(cancer, "gene KM"), {
    rows <- list()
    for (gn in pc$km_genes) for (ep in pc$km_endpoints) {
      de <- paper_endpoint_data(d0, ep, cancer)
      gcol <- paper_gene_col(de, gn, "_group", cancer)
      de$group <- de[[gcol]]; de <- de[!is.na(de$group), ]
      p_lr <- paper_logrank(de)
      fit <- survfit(Surv(ep_months, ep_event) ~ group, data = de)
      ks <- km_stats(fit, cfg$km_times)
      rows[[length(rows) + 1]] <- data.frame(gene = gn, endpoint = ep, ks, logrank_p = p_lr)
      created <<- c(created, save_fig_q(
        paper_km_plot(de, paste0(gn, " expression — TCGA-", cancer), paste(gn, "expression"), ep_ylab(ep), p_lr),
        cancer, paste0("km_", gn, "_", ep), width = 7, height = 7, subdir = top))
    }
    sm <- do.call(rbind, rows)
    p <- save_table(sm, cancer, "km_genes_summary.csv", top)
    created <<- c(created, p)
    add_info(p, "Gene-level KM (High vs Low, median split)", "survfit/survdiff(Surv(time, event) ~ gene group)",
             c(paste0("km_genes=", paste(pc$km_genes, collapse = ",")), "km_times"), nrow(d0), data_txt)
    cat("  7. 유전자 KM:", nrow(sm) / 2, "곡선 (log-rank p < 0.05:",
        paste(unique(paste(sm$gene, sm$endpoint)[sm$logrank_p < 0.05]), collapse = ", "), ")\n")
  })

  # ---- 점수별 ----
  for (s in names(pc$scores)) {
    sub <- file.path(top, s)
    slab <- pc$scores[[s]]$label
    cat("\n  --", s, "(", slab, ")", if (s == names(pc$scores)[1]) "[메인]" else "[민감도]", "--\n")
    d <- d0
    d$score <- d[[paste0("score_", s)]]
    d$group <- d[[paste0("group_", s)]]
    if (is.null(d$score) || is.null(d$group)) { failures <- c(failures, paste(cancer, s, ": 점수 열 없음 → 20 재실행")); next }
    labels <- c(term_labels(unique(c(cov_terms, "purity_pct_per10"))),
                group = paste0(slab, " group"), score_sd = paste0(slab, " (per 1 SD)"))
    score_txt <- paste0(slab, if (pc$scores[[s]]$source == "sav") " (mean z-score of the 5 genes, .sav)" else
                        paste0(" (", pc$scores[[s]]$method, ", 11)"))

    # 1·2. 상관
    run(paste(cancer, s, "correlation"), {
      genes <- c(pc$cor_genes[[cancer]], pc$cor_sig_genes)
      vars <- c(setNames("score", slab), setNames(names(pc$cor_clinical), pc$cor_clinical),
                setNames(vapply(genes, function(g) paper_gene_col(d, g, where = cancer), ""), genes))
      for (m in c("pearson", "spearman")) {
        long <- paper_cor(d, vars, m)
        p1 <- save_table(long, cancer, paste0("cor_", m, "_long.csv"), sub)
        p2 <- save_df_table(paper_cor_wide(long, names(vars)), cancer, paste0("cor_", m), subdir = sub,
                            landscape = TRUE, font_size = 6, caption = paste0(
          if (m == "pearson") "Pearson" else "Spearman", " correlation, TCGA-", cancer,
          ". Lower triangle: r (p); upper triangle: pairwise n. Gene expression as log2(expression + 1). Score: ", score_txt,
          ". Not adjusted for multiple comparisons."))
        created <<- c(created, p1, p2)
        add_info(c(p1, p2), paste(m, "correlation"), paste0("cor.test(method = '", m, "'), pairwise complete"),
                 c(paste0("cor_clinical=", paste(names(pc$cor_clinical), collapse = ",")), paste0("cor_genes=", paste(genes, collapse = ","))), nrow(d), data_txt, "stats")
      }
      sw <- paper_shapiro(d, vars)
      p3 <- save_df_table(transform(sw, W = sprintf("%.3f", W), p = fmt_p(p)), cancer, "shapiro", subdir = sub,
                          caption = paste0("Shapiro–Wilk normality test, TCGA-", cancer, " (p < 0.05 = non-normal)."))
      created <<- c(created, p3)
      add_info(p3, "Shapiro–Wilk", "shapiro.test", character(), nrow(d), data_txt, "stats")
      cat("  1·2. 상관 변수", length(vars), "개; 정규성 기각 (p < 0.05):", sum(!sw$normal, na.rm = TRUE), "\n")
    })

    # 4. Table 1
    run(paste(cancer, s, "table1"), {
      vars <- intersect(pc$table1_vars, names(d))
      miss <- setdiff(pc$table1_vars, names(d))
      if (length(miss)) stop("Table 1 변수 없음: ", paste(miss, collapse = ", "), paper_where("table1_vars"))
      t1 <- paper_table1(d[!is.na(d$group), ], "group", vars, labels_all)
      tests <- t1$tests
      rows <- t1$rows
      rows$p <- ""; rows$Test <- ""
      first <- !duplicated(rows$Variable)
      rows$p[first] <- fmt_p(tests$p[match(rows$Variable[first], tests$label)])
      rows$Test[first] <- tests$test[match(rows$Variable[first], tests$label)]
      rows$`Missing (n)` <- ""
      rows$`Missing (n)`[first] <- tests$missing_total[match(rows$Variable[first], tests$label)]
      rows$`Expected < 5` <- ""
      rows$`Expected < 5`[first] <- ifelse(tests$any_expected_lt5[match(rows$Variable[first], tests$label)], "Yes", "No")
      rows$Variable[!first] <- ""
      ng <- table(d$group)
      names(rows)[names(rows) %in% names(ng)] <- paste0(names(rows)[names(rows) %in% names(ng)], " (n = ", ng[names(rows)[names(rows) %in% names(ng)]], ")")
      p1 <- save_df_table(rows, cancer, "table1_signature", subdir = sub, font_size = 8, caption = paste0(
        "Clinical characteristics by ", slab, " group (High vs Low, median split), TCGA-", cancer,
        ". n (% of non-missing). Chi-square test without continuity correction; Fisher's exact test when any expected count < 5. ",
        "Missing values excluded from percentages and tests."))
      p2 <- save_table(tests, cancer, "table1_tests.csv", sub)
      created <<- c(created, p1, p2)
      add_info(c(p1, p2), "Table 1 (High vs Low signature group)", "chisq.test(correct = FALSE) / fisher.test",
               paste0("table1_vars=", paste(vars, collapse = ",")), sum(!is.na(d$group)), data_txt, "stats")
      cat("  4. Table 1: p < 0.05 →", paste(tests$label[tests$p < 0.05], collapse = ", "), "\n")
    })

    # 5. Cox (종점별)
    for (ep in pc$cox_endpoints) run(paste(cancer, s, "Cox", ep), {
      de <- paper_endpoint_data(d, ep, cancer)
      de <- de[!is.na(de$group), ]
      uni <- do.call(rbind, lapply(c("group", cov_terms), function(t) cox_fit_rows(de, t, labels, "Univariable")))
      mul <- cox_fit_rows(de, c("group", cov_terms), labels, "Multivariable")
      sen <- cox_fit_rows(de, c("group", sens_terms), labels, paste0("Multivariable without ", paste(labels[setdiff(cov_terms, sens_terms)], collapse = ", ")))
      raw <- rbind(uni, mul, sen)
      raw$endpoint <- ep
      p0 <- save_table(raw, cancer, paste0("cox_", ep, "_raw.csv"), sub)
      cap <- function(what, x) paste0(ep_label(ep), ", ", what, " Cox regression, TCGA-", cancer, ". ", slab,
                                       " group: ", group_contrast_label(), " (median split). ", x)
      p1 <- save_df_table(fmt_model_tab(uni), cancer, paste0("cox_uni_", ep), subdir = sub, caption = cap(
        "univariable", "Each variable in its own model; n and events are per model (complete cases)."))
      p2 <- save_df_table(fmt_model_tab(mul), cancer, paste0("cox_multi_", ep), subdir = sub, caption = cap(
        "multivariable", sprintf("All variables in one model, complete cases: n = %d, events = %d, EPV = %.1f.",
                                 mul$model_n[1], mul$model_events[1], mul$epv[1])))
      p3 <- save_df_table(fmt_model_tab(sen), cancer, paste0("cox_multi_nopy_", ep), subdir = sub, caption = cap(
        "sensitivity multivariable", sprintf("Same model without %s (larger complete-case sample): n = %d, events = %d, EPV = %.1f.",
                                             paste(tolower(labels[setdiff(cov_terms, sens_terms)]), collapse = ", "),
                                             sen$model_n[1], sen$model_events[1], sen$epv[1])))
      fp <- save_fig_q(forest_plot(mul, paste0(ep_label(ep), " — TCGA-", cancer, ", multivariable Cox"),
                                 sprintf("%s group %s; n = %d, events = %d", slab, group_contrast_label(),
                                         mul$model_n[1], mul$model_events[1])),
                     cancer, paste0("forest_", ep), width = 10, height = 1.2 + 0.35 * nrow(mul), subdir = sub)
      p_lr <- paper_logrank(de)
      km <- save_fig_q(paper_km_plot(de, paste0(slab, " — TCGA-", cancer), paste(slab, "group"), ep_ylab(ep), p_lr),
                     cancer, paste0("km_signature_", ep), width = 7, height = 7, subdir = sub)
      created <<- c(created, p0, p1, p2, p3, fp, km)
      add_info(c(p0, p1, p2, p3), paste("Cox regression,", ep),
               list(paste0("coxph(Surv(", ep, "_months, ", ep, "_event) ~ group + ", paste(cov_terms, collapse = " + "), ")"),
                    "univariable: one term per model", paste0("sensitivity: without ", paste(pc$sens_drop, collapse = ","))),
               c(paste0("cox_vars=", paste(pc$cox_vars, collapse = ",")), paste0("sens_drop=", paste(pc$sens_drop, collapse = ",")), "group_levels", "reference_levels"),
               nrow(de), data_txt)
      g <- mul[mul$term == "group" & !mul$reference, ]
      gu <- uni[uni$term == "group" & !uni$reference, ]
      cat(sprintf("  5. %s: n %d (사건 %d), log-rank p %s | 단변량 HR %s p %s | 다변량 (n %d) HR %s p %s\n", ep, nrow(de),
                  sum(de$ep_event), fmt_p(p_lr), fmt_hr(gu$hr, gu$lo, gu$hi), fmt_p(gu$p), g$model_n,
                  fmt_hr(g$hr, g$lo, g$hi), fmt_p(g$p)))
    })

    # 8. purity
    run(paste(cancer, s, "purity"), {
      pv <- c(setNames("score", slab), `Tumor purity (ESTIMATE)` = "purity", `ESTIMATE score` = "estimate_score",
              `Stromal score` = "stromal_score", `Immune score` = "immune_score")
      pcor <- rbind(paper_cor(d, pv, "pearson"), paper_cor(d, pv, "spearman"))
      pcor <- pcor[pcor$var1 == slab & pcor$var2 != slab, ]
      p1 <- save_table(pcor, cancer, "purity_cor.csv", sub)
      pg <- pc$purity$partial_genes[[cancer]]
      part <- do.call(rbind, lapply(pg, function(gn) {
        gcol <- paper_gene_col(d, gn, where = cancer)
        do.call(rbind, lapply(c("pearson", "spearman"), function(m) {
          un <- paper_cor(d, setNames(c("score", gcol), c("a", "b")), m)
          un <- un[un$var1 == "a" & un$var2 == "b", ]
          pp <- paper_partial_cor(d, "score", gcol, "purity", m)
          data.frame(gene = gn, method = m, r = un$r, p = un$p, n = un$n,
                     partial_r = pp$r, partial_p = pp$p, partial_n = pp$n)
        }))
      }))
      p2 <- save_table(part, cancer, "purity_partial_cor.csv", sub)
      tab <- data.frame(Gene = part$gene, Method = tools::toTitleCase(part$method),
                        `r (unadjusted)` = sprintf("%.3f", part$r), p = fmt_p(part$p), n = part$n,
                        `Partial r (purity-adjusted)` = sprintf("%.3f", part$partial_r), `Partial p` = fmt_p(part$partial_p),
                        `Partial n` = part$partial_n, check.names = FALSE)
      p3 <- save_df_table(tab, cancer, "purity_partial_cor", subdir = sub, caption = paste0(
        slab, " vs gene expression (log2), unadjusted and partial correlation adjusted for ESTIMATE tumor purity, TCGA-", cancer,
        ". Partial correlation = correlation of residuals after regressing both variables on purity (Spearman: on ranks)."))
      created <<- c(created, p1, p2, p3)
      add_info(c(p1, p2, p3), "Tumor purity correlation / partial correlation",
               "cor.test; partial = cor(resid(x ~ purity), resid(y ~ purity))",
               paste0("partial_genes=", paste(pg, collapse = ",")), nrow(d), data_txt, c("stats", "tidyestimate"))
      for (ep in pc$purity$cox_endpoints) {
        de <- paper_endpoint_data(d, ep, cancer); de <- de[!is.na(de$group), ]
        mp <- cox_fit_rows(de, c("group", cov_terms, "purity_pct_per10"), labels, "Multivariable + tumor purity")
        p4 <- save_df_table(fmt_model_tab(mp), cancer, paste0("cox_multi_purity_", ep), subdir = sub, caption = paste0(
          ep_label(ep), ", multivariable Cox with ESTIMATE tumor purity added, TCGA-", cancer,
          sprintf(". Complete cases: n = %d, events = %d, EPV = %.1f.", mp$model_n[1], mp$model_events[1], mp$epv[1])))
        created <<- c(created, p4)
        add_info(p4, paste("Cox + purity,", ep), paste0("coxph(... ~ group + ", paste(cov_terms, collapse = " + "), " + purity_pct_per10)"),
                 "covariate_scale", nrow(de), data_txt)
        g <- mp[mp$term == "group" & !mp$reference, ]
        cat(sprintf("  8. %s 다변량 + purity: HR %s p %s (n %d)\n", ep, fmt_hr(g$hr, g$lo, g$hi), fmt_p(g$p), g$model_n))
      }
      r <- pcor[pcor$var2 == "Tumor purity (ESTIMATE)", ]
      cat(sprintf("  8. 점수 ~ purity: Pearson r %.3f (p %s), Spearman r %.3f (p %s)\n", r$r[1], fmt_p(r$p[1]), r$r[2], fmt_p(r$p[2])))
    })

    # 10. cut-off 민감도
    for (ep in pc$cutoff_endpoints) run(paste(cancer, s, "cutoff", ep), {
      de <- paper_endpoint_data(d, ep, cancer)
      de <- de[!is.na(de$score), ]
      de$score_sd <- (de$score - mean(de$score)) / stats::sd(de$score)
      cu <- rbind(cox_fit_rows(de, "score_sd", labels, "Univariable"),
                  cox_fit_rows(de, c("score_sd", cov_terms), labels, "Multivariable"))
      p1 <- save_df_table(fmt_model_tab(cu), cancer, paste0("cox_continuous_", ep), subdir = sub, caption = paste0(
        ep_label(ep), ", ", slab, " as a continuous variable (HR per 1 SD within the survival cohort), TCGA-", cancer, "."))
      p0 <- save_table(cu, cancer, paste0("cox_continuous_", ep, "_raw.csv"), sub)
      created <<- c(created, p0, p1)
      rows <- list()
      for (cn in names(pc$cutpoints)) {
        k <- pc$cutpoints[[cn]]
        # 분위 경계 = 점수 샘플 전체 (median split과 같은 기준)
        qd <- d[!is.na(d$score), c("sample_id", "score")]
        qd$q <- paper_quantile_group(qd$score, k, substr(cn, 1, 1))
        dq <- de; dq$group <- qd$q[match(dq$sample_id, qd$sample_id)]
        labq <- c(labels, group = paste(slab, tolower(cn)))
        p_lr <- paper_logrank(dq)
        r <- rbind(cox_fit_rows(dq, "group", labq, "Univariable"),
                   cox_fit_rows(dq, c("group", cov_terms), labq, "Multivariable"))
        r <- r[r$term == "group", ]
        rows[[cn]] <- data.frame(cut = cn, logrank_p = p_lr, r)
        created <<- c(created, save_fig_q(paper_km_plot(dq, paste0(slab, " ", tolower(cn), "s — TCGA-", cancer),
                                                      paste(slab, tolower(cn)), ep_ylab(ep), p_lr),
                                        cancer, paste0("km_", tolower(cn), "_", ep), width = 7, height = 7.5, subdir = sub))
      }
      cs <- do.call(rbind, rows)
      p2 <- save_table(cs, cancer, paste0("cutoff_summary_", ep, "_raw.csv"), sub)
      tab <- data.frame(`Cut-off` = ifelse(duplicated(paste(cs$cut, cs$model)), "", paste(cs$cut, "-", cs$model)),
                        `Log-rank p` = ifelse(duplicated(cs$cut), "", fmt_p(cs$logrank_p)),
                        format_cox(cs)[, c("Level", "N", "Events", "HR (95% CI)", "p", "Overall p")], check.names = FALSE)
      p3 <- save_df_table(tab, cancer, paste0("cutoff_summary_", ep), subdir = sub, caption = paste0(
        ep_label(ep), " by ", slab, " tertiles and quartiles (lowest group = reference), TCGA-", cancer,
        ". Multivariable models adjusted for ", adjust_text(pc$cox_vars), " (complete cases). Overall p = likelihood-ratio test."))
      created <<- c(created, p2, p3)
      add_info(c(p0, p1, p2, p3), paste("Cut-off sensitivity,", ep), "coxph(... ~ score_sd [+ covariates]); quantile groups",
               paste0("cutpoints=", paste(names(pc$cutpoints), pc$cutpoints, sep = ":", collapse = ",")), nrow(de), data_txt)
      g <- cu[cu$model == "Multivariable" & cu$term == "score_sd", ]
      cat(sprintf("  10. %s 연속형 다변량 HR %s p %s; log-rank 삼분위 %s, 사분위 %s\n", ep, fmt_hr(g$hr, g$lo, g$hi), fmt_p(g$p),
                  fmt_p(cs$logrank_p[cs$cut == names(pc$cutpoints)[1]][1]), fmt_p(cs$logrank_p[cs$cut == names(pc$cutpoints)[2]][1])))
    })
  }

  if (length(info)) created <- c(created, write_analysis_info(info, cancer, "21", top))
}

cat("\n생성된 파일", length(created), "개\n")
if (length(failures)) {
  cat("\n실패:\n", paste0("  ", failures, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("실패한 분석 없음\n")
