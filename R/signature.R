# =============================================================
# R/signature.R
# signature score 트랙 (11_sig_scores.R, 12_sig_survival.R, 13_sig_geo.R) 공통 함수
# 사용: source("config.R"); source("R/utils.R"); source("R/signature.R")
# 설정: cfg$signatures, cfg$signature_* (config.R)
# 기존 유전자 파이프라인(merged.rds, detect_genes, 04 BH, 05, 90)과 분리: 이 파일의 함수는 그쪽에서 쓰지 않음
# 점수 파일 = 공통 인터페이스: sample_id + sig_methods 열 (생존 정보 없음)
# =============================================================

# 점수 방법 (점수 파일 열 이름) / 표·그림 라벨
sig_methods <- c("z", "pc1", "ssgsea", "singscore")
sig_method_labels <- c(z = "Mean z-score", pc1 = "PC1", ssgsea = "ssGSEA", singscore = "singscore")

# ---- 정의 / 경로 -------------------------------------------------------

# cfg$signatures$<SIG> 검사 → list(sig, up, down, source, label, where)
# 중복 유전자, up/down 겹침, 빈 정의, source 없음 → 중단 (조용히 고치지 않음)
sig_def <- function(sig) {
  where <- paste0("cfg$signatures$", sig)
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", sig)) {
    stop("signature 이름 '", sig, "'은 파일명에 쓰므로 영문/숫자/_ . -만 가능 → cfg$signatures 이름 수정")
  }
  if (toupper(sig) %in% toupper(cfg$cancers)) stop("signature 이름 '", sig, "'이 암종 코드와 같음 → ", where, " 이름 수정")
  s <- cfg$signatures[[sig]]
  if (is.null(s)) {
    stop(where, " 없음 (정의된 signature: ",
         if (length(cfg$signatures)) paste(names(cfg$signatures), collapse = ", ") else "없음", ")")
  }
  bad_keys <- setdiff(names(s), c("up", "down", "source", "label"))
  if (length(bad_keys)) stop(where, ": 알 수 없는 항목 ", paste(bad_keys, collapse = ", "), " (up, down, source, label만 사용)")
  get_set <- function(k) {
    x <- s[[k]]
    if (is.null(x)) return(character())
    if (!is.character(x) || anyNA(x)) stop(where, "$", k, ": 유전자 기호 문자 벡터여야 함 (NA 불가)")
    if (any(x != trimws(x) | !nzchar(x))) {
      stop(where, "$", k, ": 공백이 붙었거나 빈 유전자 기호 → ", paste0("\"", x[x != trimws(x) | !nzchar(x)], "\"", collapse = ", "))
    }
    d <- unique(x[duplicated(x)])
    if (length(d)) stop(where, "$", k, ": 중복 유전자 ", paste(d, collapse = ", "))
    x
  }
  up <- get_set("up"); down <- get_set("down")
  if (!length(up) && !length(down)) stop(where, ": up, down 모두 비어 있음")
  both <- intersect(up, down)
  if (length(both)) stop(where, ": up과 down에 모두 있는 유전자 ", paste(both, collapse = ", "))
  if (is.null(s$source) || !is.character(s$source) || !nzchar(trimws(s$source))) {
    stop(where, "$source 필수 (출처 문헌/MSigDB 이름; analysis_info와 Methods에 기록)")
  }
  list(sig = sig, up = up, down = down, source = s$source, label = s$label %||% sig, where = where)
}

# signature 출력 하위 폴더: output/<type>/<C>/signature/<SIG> (GEO: GEO/<GSE>/signature/<SIG>)
sig_subdir <- function(sig, parent = NULL) do.call(file.path, as.list(c(parent, "signature", sig)))

# 점수 파일: data/processed/<prefix>_sigscores_<SIG>.csv (prefix = 암종 또는 GSE)
sig_scores_path <- function(prefix, sig) processed_path(prefix, paste0("sigscores_", sig, ".csv"))

# ---- 발현 --------------------------------------------------------------

# 02a의 log2(TPM+1) 행렬 (database/TCGA_<C>_RNAseq_Expression.csv) → 유전자 × 환자 matrix
sig_tcga_expr_path <- function(cancer) file.path(cfg$clinical_dir, paste0("TCGA_", cancer, "_RNAseq_Expression.csv"))
sig_read_tcga_expr <- function(cancer) {
  f <- sig_tcga_expr_path(cancer)
  if (!file.exists(f)) stop(f, " 없음 → 02a_download_tcga.R 먼저 실행")
  d <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE)
  if (names(d)[1] != "Gene") stop(f, ": 첫 열이 Gene이 아님")
  dup <- unique(d$Gene[duplicated(d$Gene)])
  if (length(dup)) stop(f, ": 중복 유전자 기호 ", paste(head(dup, 10), collapse = ", "), " → 02a 재실행")
  m <- as.matrix(d[, -1])
  rownames(m) <- d$Gene
  m
}

# 저발현 필터 (rule = cfg$signature_expr_filter$tcga/geo; NULL = 분산 0/결측만 제외)
# 반환: list(expr = 통과 유전자 행렬, stats = 모든 유전자의 mean_log2, frac_expressed, sd, kept)
# 통과 행렬 하나를 네 점수와 무작위 gene set 배경에 공통으로 사용
sig_filter_genes <- function(expr, rule) {
  sdv  <- apply(expr, 1, stats::sd, na.rm = TRUE)
  frac <- if (is.null(rule)) rep(NA_real_, nrow(expr)) else rowMeans(expr >= rule$min_log2, na.rm = TRUE)
  keep <- !is.na(sdv) & sdv > 0 & rowSums(is.na(expr)) == 0
  if (!is.null(rule)) keep <- keep & frac >= rule$min_frac
  stats <- data.frame(gene = rownames(expr), mean_log2 = rowMeans(expr, na.rm = TRUE), frac_expressed = frac,
                      sd = sdv, kept = keep, row.names = NULL, stringsAsFactors = FALSE)
  list(expr = expr[keep, , drop = FALSE], stats = stats)
}

sig_filter_text <- function(rule) {
  if (is.null(rule)) "non-zero variance, no missing values"
  else sprintf("log2(TPM+1) >= %g in >= %g%% of samples, non-zero variance", rule$min_log2, 100 * rule$min_frac)
}

# ---- 유전자 매칭 -------------------------------------------------------

# 정의 유전자 → status: matched / <missing_label> (행렬·플랫폼에 없음) / low_expression (필터 탈락)
# 방향별·전체 사용 비율 < cfg$signature_min_match 또는 사용 수 < cfg$signature_min_genes → 중단
# alias 자동 변환 없음 (정확한 기호만)
sig_match <- function(def, all_genes, kept_genes, where, missing_label = "not_in_matrix") {
  tab <- data.frame(gene = c(def$up, def$down),
                    direction = rep(c("up", "down"), c(length(def$up), length(def$down))), stringsAsFactors = FALSE)
  tab$status <- ifelse(!tab$gene %in% all_genes, missing_label,
                       ifelse(tab$gene %in% kept_genes, "matched", "low_expression"))
  used <- tab$status == "matched"
  msgs <- character()
  for (dir in c(unique(tab$direction), if (length(unique(tab$direction)) > 1) "all")) {
    i <- if (dir == "all") rep(TRUE, nrow(tab)) else tab$direction == dir
    r <- mean(used[i])
    if (r < cfg$signature_min_match) {
      msgs <- c(msgs, sprintf("%s: 사용 %d/%d (%.0f%% < %.0f%%)", dir, sum(used[i]), sum(i), 100 * r,
                              100 * cfg$signature_min_match))
    }
    if (dir != "all" && sum(used[i]) < cfg$signature_min_genes) {
      msgs <- c(msgs, sprintf("%s: 사용 유전자 %d개 < cfg$signature_min_genes %d", dir, sum(used[i]), cfg$signature_min_genes))
    }
  }
  if (length(msgs)) {
    miss <- tab[!used, ]
    stop(where, " ", def$sig, ": 매칭 부족 — ", paste(msgs, collapse = "; "),
         "\n  사용 못 한 유전자: ", paste(sprintf("%s (%s, %s)", miss$gene, miss$direction, miss$status), collapse = ", "),
         "\n  → ", def$where, " 유전자 기호 확인 (alias 자동 변환 없음) 또는 cfg$signature_min_match / cfg$signature_min_genes")
  }
  tab
}

# 매칭 표 → 사용 유전자만 남긴 정의
sig_used_def <- function(def, tab) {
  def$up   <- tab$gene[tab$direction == "up" & tab$status == "matched"]
  def$down <- tab$gene[tab$direction == "down" & tab$status == "matched"]
  def
}
sig_sets <- function(def) Filter(length, list(up = def$up, down = def$down))

# up/down 집합 점수 → up − down (한쪽만 있으면 up 또는 −down)
sig_combine <- function(s) (s$up %||% 0) - (s$down %||% 0)

# ---- 점수 4종 (expr = 필터 통과 행렬, 열 = 점수 계산 샘플; def = sig_used_def 결과) ----------

sig_zmat <- function(expr, genes) t(scale(t(expr[genes, , drop = FALSE])))

sig_score_z <- function(expr, def) {
  sig_combine(lapply(sig_sets(def), function(g) colMeans(sig_zmat(expr, g))))
}

# PC1 (표준화 유전자의 prcomp), 집합별로 z-score와 양의 상관이 되도록 부호 맞춤 → up − down
# 반환: list(score, var_explained (집합별), loading (유전자별, 부호 맞춘 값))
sig_score_pc1 <- function(expr, def) {
  res <- lapply(sig_sets(def), function(g) {
    z <- sig_zmat(expr, g)
    if (length(g) == 1) return(list(score = z[1, ], var = 1, loading = setNames(1, g)))
    pr <- stats::prcomp(t(z), center = TRUE, scale. = TRUE)
    pc <- pr$x[, 1]
    flip <- if (stats::cor(pc, colMeans(z)) < 0) -1 else 1
    list(score = flip * pc, var = pr$sdev[1]^2 / sum(pr$sdev^2), loading = flip * pr$rotation[, 1])
  })
  list(score = sig_combine(lapply(res, `[[`, "score")),
       var_explained = vapply(res, `[[`, 0, "var"),
       loading = unlist(unname(lapply(res, `[[`, "loading"))))
}

# ssGSEA (GSVA): up/down을 한 번에 계산 (정규화를 쓰더라도 두 집합에 같은 비율) → up − down
sig_score_ssgsea <- function(expr, def) {
  sets <- sig_sets(def)
  p <- GSVA::ssgseaParam(exprData = expr, geneSets = sets, minSize = 1,
                         normalize = isTRUE(cfg$signature_ssgsea_normalize))
  es <- suppressMessages(GSVA::gsva(p, verbose = FALSE))
  sig_combine(lapply(setNames(names(sets), names(sets)), function(k) es[k, colnames(expr)]))
}

# singscore: 집합마다 simpleScore(upSet = 집합, centerScore = TRUE) → up − down (합의된 정의 그대로)
sig_score_singscore <- function(expr, def, ranks = singscore::rankGenes(expr)) {
  sig_combine(lapply(sig_sets(def), function(g) {
    s <- singscore::simpleScore(ranks, upSet = g, centerScore = TRUE)
    setNames(s$TotalScore, rownames(s))[colnames(expr)]
  }))
}

# 점수 파일 (공통 인터페이스) + PC1 설명력/적재값
# ranks: singscore::rankGenes(expr) (여러 signature에서 재사용; NULL이면 계산)
sig_scores <- function(expr, def, ranks = NULL) {
  pc <- sig_score_pc1(expr, def)
  if (is.null(ranks)) ranks <- singscore::rankGenes(expr)
  scores <- data.frame(sample_id = colnames(expr),
                       z = unname(sig_score_z(expr, def)),
                       pc1 = unname(pc$score),
                       ssgsea = unname(sig_score_ssgsea(expr, def)),
                       singscore = unname(sig_score_singscore(expr, def, ranks)),
                       stringsAsFactors = FALSE)
  bad <- scores$sample_id[!stats::complete.cases(scores)]
  if (length(bad)) stop(def$sig, ": 점수 결측 샘플 ", paste(head(bad, 10), collapse = ", "))
  attr(scores, "pc1") <- pc
  scores
}

# ---- QC (outcome 사용 안 함) --------------------------------------------

# Cronbach's alpha (psych::alpha, check.keys = FALSE): down 유전자는 z × −1
# 반환: list(alpha, items = data.frame(gene, r_drop, alpha_if_dropped)); 유전자 2개 미만이면 NA
sig_alpha <- function(expr, def) {
  genes <- c(def$up, def$down)
  if (length(genes) < 2) return(list(alpha = NA_real_, items = data.frame(gene = genes, r_drop = NA, alpha_if_dropped = NA)))
  z <- sig_zmat(expr, genes)
  z[def$down, ] <- -z[def$down, ]
  a <- suppressMessages(suppressWarnings(psych::alpha(as.data.frame(t(z)), check.keys = FALSE, warnings = FALSE)))
  items <- data.frame(gene = genes, r_drop = unname(a$item.stats[genes, "r.drop"]),
                      alpha_if_dropped = if (length(genes) > 2) unname(a$alpha.drop[genes, "raw_alpha"]) else NA_real_,
                      stringsAsFactors = FALSE)
  list(alpha = unname(a$total$raw_alpha), items = items)
}

# PC1 설명력 (상관행렬 최대 고윳값 / 유전자 수 = prcomp(scale.=TRUE) 설명력과 같음)
sig_pc1_var <- function(expr, genes) {
  if (length(genes) < 2) return(1)
  ev <- eigen(stats::cor(t(expr[genes, , drop = FALSE])), symmetric = TRUE, only.values = TRUE)$values
  ev[1] / length(genes)
}

# 같은 크기 무작위 gene set (배경 = 필터 통과 유전자) 대비 PC1 설명력 백분위
# 반환: list(null, pct = 100 × mean(null < obs), p = (1 + #(null ≥ obs)) / (n + 1))
sig_random_pc1 <- function(expr, size, obs, n = cfg$signature_random_n, seed = cfg$seed) {
  if (size < 2) return(list(null = numeric(), pct = NA_real_, p = NA_real_))
  set.seed(seed)
  bg <- rownames(expr)
  null <- vapply(seq_len(n), function(i) sig_pc1_var(expr, sample(bg, size)), 0)
  list(null = null, pct = 100 * mean(null < obs), p = (1 + sum(null >= obs)) / (n + 1))
}

# 점수 간 Spearman 상관 (long: method1, method2, spearman_rho, n)
sig_score_cor <- function(scores) {
  m <- as.matrix(scores[, sig_methods])
  r <- stats::cor(m, method = "spearman")
  g <- expand.grid(i = seq_along(sig_methods), j = seq_along(sig_methods))
  data.frame(method1 = sig_methods[g$i], method2 = sig_methods[g$j], spearman_rho = r[cbind(g$i, g$j)],
             n = nrow(m), stringsAsFactors = FALSE)
}

# QC 묶음: 유전자 표, 집합별 요약, 상관, 무작위 분포
# tab = sig_match 결과, stats = sig_filter_genes()$stats, expr = 필터 통과 행렬, def = 사용 정의
sig_qc <- function(expr, def, tab, stats, scores) {
  pc <- attr(scores, "pc1")
  al <- sig_alpha(expr, def)
  genes <- merge(tab, stats[, c("gene", "mean_log2", "frac_expressed", "sd")], by = "gene", all.x = TRUE, sort = FALSE)
  genes <- genes[match(tab$gene, genes$gene), ]
  genes$r_drop           <- al$items$r_drop[match(genes$gene, al$items$gene)]
  genes$alpha_if_dropped <- al$items$alpha_if_dropped[match(genes$gene, al$items$gene)]
  genes$pc1_loading      <- unname(pc$loading[genes$gene])
  rownames(genes) <- NULL

  # 집합별 행 (up, down) + 둘 다 있으면 all (alpha만; 결합 PC1은 보고하지 않음)
  sets <- sig_sets(def)
  null <- list()
  count_row <- function(i) data.frame(n_defined = sum(i), n_matched = sum(i & tab$status %in% c("matched", "low_expression")),
                                      n_used = sum(i & tab$status == "matched"), match_ratio = mean(tab$status[i] == "matched"))
  rows <- lapply(names(sets), function(k) {
    obs <- unname(pc$var_explained[[k]])
    rp <- sig_random_pc1(expr, length(sets[[k]]), obs)
    null[[k]] <<- rp$null
    one <- list(up = sets[[k]], down = character())   # 방향 안에서는 부호를 뒤집지 않음
    data.frame(set = k, count_row(tab$direction == k), alpha = sig_alpha(expr, one)$alpha,
               pc1_var_explained = obs, random_pc1_pct = rp$pct, random_pc1_p = rp$p, stringsAsFactors = FALSE)
  })
  if (length(sets) > 1) {
    rows <- c(rows, list(data.frame(set = "all", count_row(rep(TRUE, nrow(tab))), alpha = al$alpha,
                                    pc1_var_explained = NA_real_, random_pc1_pct = NA_real_, random_pc1_p = NA_real_)))
  }
  list(genes = genes, summary = do.call(rbind, rows), cor = sig_score_cor(scores), null = null)
}

# ---- 그림 ---------------------------------------------------------------

# 점수 4종 쌍별 산점도 (Spearman rho 표시)
sig_pairs_plot <- function(scores, title) {
  pr <- utils::combn(sig_methods, 2)
  d <- do.call(rbind, lapply(seq_len(ncol(pr)), function(k) {
    a <- pr[1, k]; b <- pr[2, k]
    rho <- stats::cor(scores[[a]], scores[[b]], method = "spearman")
    data.frame(pair = sprintf("%s vs %s (rho = %.2f)", sig_method_labels[[b]], sig_method_labels[[a]], rho),
               x = scores[[a]], y = scores[[b]])
  }))
  ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.5) +
    ggplot2::facet_wrap(~pair, scales = "free") +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::theme_classic(base_size = 10)
}

# 무작위 gene set PC1 설명력 분포 + 관측값
sig_random_plot <- function(null, summary, title) {
  sets <- names(null)[lengths(null) > 0]
  if (!length(sets)) return(NULL)
  d <- do.call(rbind, lapply(sets, function(k) data.frame(set = k, v = null[[k]])))
  o <- summary[summary$set %in% sets, ]
  o$lab <- sprintf("observed %.3f (percentile %.1f)", o$pc1_var_explained, o$random_pc1_pct)
  ggplot2::ggplot(d, ggplot2::aes(v)) +
    ggplot2::geom_histogram(bins = 40, fill = "grey75", colour = "white") +
    ggplot2::geom_vline(data = o, ggplot2::aes(xintercept = pc1_var_explained), colour = "#C8442F", linewidth = 0.8) +
    ggplot2::geom_text(data = o, ggplot2::aes(x = pc1_var_explained, y = Inf, label = lab), vjust = 1.5, hjust = 1.05,
                       size = 3, colour = "#C8442F") +
    ggplot2::facet_wrap(~set, scales = "free") +
    ggplot2::labs(x = "PC1 variance explained (random gene sets of the same size)", y = "Count", title = title) +
    ggplot2::theme_classic(base_size = 10)
}
