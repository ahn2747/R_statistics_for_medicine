# =============================================================
# 11_sig_scores.R [SIG ...] [COAD READ ...]
# signature 점수 (z, PC1, ssGSEA, singscore) + QC — 암종마다 따로 (COAD/READ 절대 합치지 않음)
#   입력: database/TCGA_<C>_RNAseq_Expression.csv (02a, log2(TPM+1)), data/processed/<C>_clinical_all.rds (01)
#   점수 샘플: clinical_all ∩ 발현 행렬 (신보조요법 제외 전; 제외는 12에서 기존 규칙대로)
#   1. 저발현 필터 (cfg$signature_expr_filter$tcga) → 통과 행렬 하나를 네 점수와 무작위 배경에 공통 사용
#   2. 유전자 매칭 (정확한 기호만; 사용 비율 < cfg$signature_min_match → 중단)
#   3. 점수: 집합별 계산 후 up − down (PC1은 z-score와 양의 상관이 되도록 부호 맞춤)
#   4. QC (outcome 사용 안 함): Cronbach's alpha, r.drop, alpha-if-dropped, PC1 설명력 vs 무작위 gene set, 점수 간 Spearman
# 인자: signature 이름 (생략 = cfg$signatures 전체), 암종 코드 (생략 = cfg$cancers)
# 결과:
#   data/processed/<C>_sigscores_<SIG>.csv     sample_id, z, pc1, ssgsea, singscore (생존 정보 없음)
#   output/tables/<C>/signature/<SIG>/          <C>_sig_genes_<SIG>, sig_qc_<SIG>, sig_score_cor_<SIG>, analysis_info_11
#   output/figures/<C>/signature/<SIG>/         <C>_sig_score_pairs_<SIG>, sig_random_pc1_<SIG>
# =============================================================

source("config.R")
source("R/utils.R")
source("R/signature.R")

args <- commandArgs(trailingOnly = TRUE)
cancers <- if (any(toupper(args) %in% cfg$cancers)) toupper(args[toupper(args) %in% cfg$cancers]) else cfg$cancers
sigs <- args[!toupper(args) %in% cfg$cancers]
all_sigs <- !length(sigs)
if (all_sigs) sigs <- names(cfg$signatures)
if (!length(sigs)) stop("cfg$signatures가 비어 있음 → config.R에 signature 정의 추가")
unknown <- setdiff(sigs, names(cfg$signatures))
if (length(unknown)) stop("cfg$signatures에 없는 signature: ", paste(unknown, collapse = ", "))
defs <- lapply(setNames(sigs, sigs), sig_def)   # 정의 오류는 계산 전에 모두 확인
cat("signature:", paste(sigs, collapse = ", "), "/ 암종:", paste(cancers, collapse = ", "), "\n")

# 파일 출처 (analysis_info Data)
md5 <- function(f) unname(tools::md5sum(f))
gdc_release <- function(cancer) {
  f <- processed_path(cancer, "gdc_provenance.json")
  if (!file.exists(f)) return("GDC release: unknown (no gdc_provenance.json)")
  paste0("GDC Data Release ", jsonlite::read_json(f)$gdc_data_release$release %||% "unknown")
}

created  <- character()
failures <- character()
console  <- list()

for (cancer in cancers) {
  cat("\n==========", cancer, "==========\n")
  expr_all <- sig_read_tcga_expr(cancer)
  clin_f <- processed_path(cancer, "clinical_all.rds")
  if (!file.exists(clin_f)) stop(clin_f, " 없음 → 01_import_clinical.R 먼저 실행")
  clin <- readRDS(clin_f)
  ids <- intersect(colnames(expr_all), clin$sample_id)
  no_expr <- setdiff(clin$sample_id, colnames(expr_all))
  cat(sprintf("발현 행렬 %d유전자 × %d명, 임상(전체) %d명 → 점수 샘플 %d명 (임상에 있으나 발현 없음 %d명, 발현만 있음 %d명)\n",
              nrow(expr_all), ncol(expr_all), nrow(clin), length(ids), length(no_expr),
              ncol(expr_all) - length(ids)))
  if (length(no_expr)) cat("  발현 없음 (점수 없음):", paste(no_expr, collapse = ", "), "\n")
  if (!length(ids)) stop(cancer, ": 임상과 발현 행렬의 환자 ID가 겹치지 않음")
  rule <- cfg$signature_expr_filter$tcga
  flt <- sig_filter_genes(expr_all[, ids, drop = FALSE], rule)
  rm(expr_all)
  cat(sprintf("저발현 필터 (%s): %d → %d유전자\n", sig_filter_text(rule), nrow(flt$stats), nrow(flt$expr)))
  ranks <- NULL   # singscore 순위 (signature 간 공유)

  data_txt <- c(sprintf("%s (md5 %s; %s)", basename(sig_tcga_expr_path(cancer)), md5(sig_tcga_expr_path(cancer)),
                        gdc_release(cancer)),
                sprintf("%s (md5 %s)", basename(clin_f), md5(clin_f)))

  for (sig in sigs) {
    def <- defs[[sig]]
    sub <- sig_subdir(sig)
    cat("\n--", sig, "(up", length(def$up), "/ down", length(def$down), ") --\n")
    r <- tryCatch({
      tab  <- sig_match(def, flt$stats$gene, rownames(flt$expr), cancer)
      used <- sig_used_def(def, tab)
      cat(sprintf("  사용 유전자: up %d/%d, down %d/%d\n", length(used$up), length(def$up), length(used$down), length(def$down)))
      lost <- tab[tab$status != "matched", ]
      if (nrow(lost)) cat("  제외:", paste(sprintf("%s (%s)", lost$gene, lost$status), collapse = ", "), "\n")

      if (is.null(ranks)) ranks <- singscore::rankGenes(flt$expr)
      scores <- sig_scores(flt$expr, used, ranks)
      qc <- sig_qc(flt$expr, used, tab, flt$stats, scores)
      qc$summary <- cbind(qc$summary, random_n = cfg$signature_random_n, seed = cfg$seed, n_samples = length(ids),
                          n_genes_background = nrow(flt$expr), filter_rule = sig_filter_text(rule),
                          n_clinical = nrow(clin), n_clinical_no_expression = length(no_expr))

      # ---- 저장 ----
      files <- character()
      sf <- sig_scores_path(cancer, sig)
      readr::write_csv(scores, sf)
      f_genes <- save_table(qc$genes, cancer, paste0("sig_genes_", sig, ".csv"), sub)
      f_qc    <- save_table(qc$summary, cancer, paste0("sig_qc_", sig, ".csv"), sub)
      f_cor   <- save_table(qc$cor, cancer, paste0("sig_score_cor_", sig, ".csv"), sub)
      files <- c(sf, f_genes, f_qc, f_cor)
      ttl <- paste0(def$label, " signature scores — TCGA-", cancer, " (n = ", length(ids), ")")
      files <- c(files, save_fig(sig_pairs_plot(scores, ttl), cancer, paste0("sig_score_pairs_", sig),
                                 width = 9, height = 6, subdir = sub))
      rp <- sig_random_plot(qc$null, qc$summary, paste0(def$label, " PC1 vs random gene sets — TCGA-", cancer))
      if (!is.null(rp)) files <- c(files, save_fig(rp, cancer, paste0("sig_random_pc1_", sig), width = 8, height = 4,
                                                   subdir = sub))

      # ---- 분석 정보 ----
      def_txt <- paste0("signatures$", sig, "=up:(", paste(def$up, collapse = ","), ");down:(",
                        paste(def$down, collapse = ","), ")")
      settings <- c(def_txt, "signature_expr_filter", "signature_min_match", "signature_min_genes")
      dat <- c(data_txt, paste0("signature source: ", def$source))
      n_txt <- sprintf("n=%d samples (clinical_all %d, no expression %d); genes used up %d/%d, down %d/%d; background %d genes",
                       length(ids), nrow(clin), length(no_expr), length(used$up), length(def$up),
                       length(used$down), length(def$down), nrow(flt$expr))
      score_model <- c("z: per-gene z over score samples, mean per set",
                       "pc1: prcomp(scale.=TRUE) per set, sign aligned to positive correlation with z",
                       paste0("ssgsea: GSVA::ssgseaParam(normalize = ", isTRUE(cfg$signature_ssgsea_normalize),
                              ", alpha = 0.25), up and down in one call"),
                       "singscore: rankGenes + simpleScore(upSet = set, centerScore = TRUE) per set",
                       "score = score(up) - score(down)")
      info <- list(
        list(Output = f_genes, Analysis = "Signature gene matching and per-gene QC (no outcome used)",
             Model = c("exact symbol match; status matched / not_in_matrix / low_expression",
                       "psych::alpha(check.keys = FALSE) on z (down genes x -1): r.drop, alpha if dropped",
                       "PC1 loading per set (sign aligned)"),
             Settings = settings, n = n_txt, Data = dat, Packages = c("psych", "stats")),
        list(Output = f_qc, Analysis = "Signature QC per set: matching, Cronbach's alpha, PC1 variance vs random gene sets",
             Model = c("PC1 variance explained = largest eigenvalue of gene correlation / k",
                       "random gene sets of the same size from the filtered background; percentile = 100 x mean(null < observed); p = (1 + #null >= observed) / (n + 1)"),
             Settings = c(settings, "signature_random_n", "seed"), n = n_txt, Data = dat, Packages = c("psych", "stats")),
        list(Output = f_cor, Analysis = "Spearman correlation between the four signature scores",
             Model = c(score_model, "cor(method = \"spearman\")"),
             Settings = c(settings, "signature_ssgsea_normalize"), n = n_txt, Data = dat,
             Packages = c("GSVA", "singscore", "stats")))
      files <- c(files, write_analysis_info(info, cancer, "11", subdir = sub))
      move_stale_outputs(cancer, c("sig_genes_", "sig_qc_", "sig_score_cor_", "sig_score_pairs_", "sig_random_pc1_",
                                   "analysis_info_11"), files, subdir = sub)
      list(files = files, qc = qc, scores = scores)
    }, error = function(e) {
      msg <- conditionMessage(e)
      failures <<- c(failures, if (startsWith(msg, paste(cancer, sig))) msg else paste0(cancer, " ", sig, ": ", msg))
      cat("  !!", sig, "실패:", conditionMessage(e), "\n")
      NULL
    })
    if (!is.null(r)) {
      created <- c(created, r$files)
      console[[paste(cancer, sig)]] <- r$qc
    }
  }
  # cfg$signatures에서 빠진 signature의 이전 결과 폴더 → signature/_stale/ (전체 실행일 때만)
  if (all_sigs) move_stale_dirs(cancer, "signature", names(cfg$signatures))
}

cat("\n생성된 파일 (", length(created), "개):\n", paste0("  ", created, collapse = "\n"), "\n", sep = "")
cat("\n===== 요약 =====\n")
for (k in names(console)) {
  s <- console[[k]]$summary
  cat(k, "\n")
  for (i in seq_len(nrow(s))) {
    cat(sprintf("  %-4s 사용 %d/%d, alpha %s, PC1 설명력 %s (무작위 대비 백분위 %s)\n", s$set[i], s$n_used[i], s$n_defined[i],
                sprintf("%.2f", s$alpha[i]), sprintf("%.3f", s$pc1_var_explained[i]), sprintf("%.1f", s$random_pc1_pct[i])))
  }
  cr <- console[[k]]$cor
  cr <- cr[match(cr$method1, sig_methods) < match(cr$method2, sig_methods), ]
  cat("  Spearman:", paste(sprintf("%s-%s %.2f", cr$method1, cr$method2, cr$spearman_rho), collapse = ", "), "\n")
}
if (length(failures)) {
  cat("\n실패한 분석 (", length(failures), "개):\n", paste0("  ", failures, collapse = "\n"), "\n", sep = "")
  quit(save = "no", status = 1)
}
cat("\n실패한 분석 없음\n")
