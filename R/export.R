# =============================================================
# R/export.R
# 학생연구자료 내보내기 (90_export.R) 공통 함수
# 사용: source("config.R"); source("R/utils.R"); source("R/export.R")
# 설정: cfg$export (config.R)
# =============================================================

# ---- 경로 --------------------------------------------------------------

# 내보내기 이름: cfg$export$use_custom_filename = TRUE면 custom_filename, 아니면 유전자
# (대상 폴더·생성 파일·복사 파일 이름의 유전자 자리에만 쓰임. 결과 선택·docx 내용은 실제 유전자)
export_name <- function(gene) {
  if (!isTRUE(cfg$export$use_custom_filename)) return(gene)
  nm <- cfg$export$custom_filename
  bad <- function(why) stop("cfg$export$custom_filename ", why, " (use_custom_filename = TRUE)", call. = FALSE)
  if (length(nm) > 1) bad(paste0("은 값 하나여야 함 (현재 ", length(nm), "개)"))
  if (length(nm) == 0 || is.na(nm) || !is.character(nm) || !nzchar(trimws(nm))) bad("이 비어 있음")
  if (nm != trimws(nm)) bad(paste0("앞뒤에 공백이 있음: \"", nm, "\""))
  if (grepl("[<>:\"/\\\\|?*%{}[:cntrl:]]", nm)) {
    bad(paste0("에 파일명에 쓸 수 없는 문자가 있음 (< > : \" / \\ | ? * % { } 제어문자): \"", nm, "\""))
  }
  if (grepl("[. ]$", nm)) bad(paste0("이 마침표/공백으로 끝남: \"", nm, "\""))
  if (grepl("^(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])$", toupper(nm))) bad(paste0("이 Windows 예약 이름: \"", nm, "\""))
  nm
}

# 대상 폴더: <dest_root>/<folder_fmt % 내보내기 이름 (export_name())>
export_dest <- function(name) file.path(cfg$export$dest_root, sprintf(cfg$export$folder_fmt, name))

# output/<type>/<cancer>/<subdir> (out_dir()과 달리 폴더를 만들지 않음 → 존재 확인용)
src_dir <- function(cancer, type, subdir = NULL) {
  do.call(file.path, as.list(c(cfg$output_dir, type, cancer, subdir)))
}

# layout 패턴 확장: {G} → 유전자, {coll} → names(cfg$gsea_collections) 각각
expand_patterns <- function(pats, gene) {
  colls <- names(cfg$gsea_collections)
  unlist(lapply(pats, function(p) {
    p <- gsub("{G}", gene, p, fixed = TRUE)
    if (grepl("{coll}", p, fixed = TRUE)) vapply(colls, function(cl) gsub("{coll}", cl, p, fixed = TRUE), "") else p
  }), use.names = FALSE)
}

# cfg$export$exclude 정규식에 걸리는 이름 (~$ 잠금 파일, _stale, inspect)
is_excluded <- function(name) {
  Reduce(`|`, lapply(cfg$export$exclude, grepl, x = name), rep(FALSE, length(name)))
}

# 폴더에서 <prefix>_<패턴>.<ext> 파일 찾기 → 패턴별 목록 (하위 폴더·제외 대상 무시)
match_outputs <- function(dir, prefix, pats, ext) {
  f <- if (dir.exists(dir)) list.files(dir) else character()
  f <- f[!dir.exists(file.path(dir, f)) & !is_excluded(f)]
  setNames(lapply(pats, function(p) {
    hit <- f[grepl(utils::glob2rx(paste0(prefix, "_", p, ".", ext)), f)]
    if (!isTRUE(cfg$export$include_de)) hit <- hit[!grepl("_de_results\\.", hit)]
    file.path(dir, hit)
  }), pats)
}

# layout 패턴 → 매칭 패턴(pat, {G} = 유전자) + 대상 패턴(dest_pat, {G} = 내보내기 이름)
# 이름 치환은 패턴의 {G} 자리에서만 (예: enrichment_* 의 경로 이름 속 유전자 기호는 바뀌지 않음)
layout_targets <- function(pats, gene, name) {
  raw <- expand_patterns(pats, "{G}")   # {coll}만 확장, {G}는 그대로
  data.frame(raw = raw, gene = gene, name = name,
             pat = gsub("{G}", gene, raw, fixed = TRUE),
             dest_pat = gsub("{G}", name, raw, fixed = TRUE),
             has_g = grepl("{G}", raw, fixed = TRUE), stringsAsFactors = FALSE)
}

rx_escape <- function(x) gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", x)

# 와일드카드(* ?)와 {G}가 함께 있는 패턴: 파일 이름에서 {G} 자리를 찾아 내보내기 이름으로 바꿈
# {G} 위치가 둘 이상으로 해석되면 (최단/최장 일치 결과가 다르면) 중단
glob_dest_name <- function(file, template, gene, name) {
  tok <- regmatches(template, gregexpr("\\{G\\}|[*?]|[^*?{]+|\\{", template))[[1]]
  slot <- tok %in% c("{G}", "*", "?")
  rebuild <- function(lazy) {
    q <- if (lazy) "?" else ""
    rx <- paste0("^", paste(vapply(tok, function(k) switch(k, "{G}" = paste0("(", rx_escape(gene), ")"),
                                                           "*" = paste0("(.*", q, ")"), "?" = "(.)", rx_escape(k)), ""),
                            collapse = ""), "$")
    m <- regmatches(file, regexec(rx, file, perl = TRUE))[[1]]   # TRE는 최단 일치(.*?)를 무시 → perl
    if (!length(m)) stop(file, ": layout 패턴 ", template, "과 일치하지 않음")
    out <- tok
    out[slot] <- m[-1]
    out[tok == "{G}"] <- name
    paste(out, collapse = "")
  }
  a <- rebuild(TRUE); b <- rebuild(FALSE)
  if (a != b) {
    stop(file, ": layout 패턴 ", template, "에서 {G} 위치를 하나로 정할 수 없음 (", a, " / ", b,
         ") → cfg$export$layout 패턴 수정 또는 cfg$export$use_custom_filename = FALSE")
  }
  a
}

# layout_targets() 결과로 파일 찾기 (매칭은 match_outputs()와 동일) + 대상 파일 이름
# 반환: data.frame(source, dest_base, pattern) — 패턴 순서, 패턴 안에서는 match_outputs 순서
match_targets <- function(dir, prefix, tg, ext) {
  hits <- match_outputs(dir, prefix, tg$pat, ext)
  rows <- lapply(seq_len(nrow(tg)), function(i) {
    f <- hits[[i]]
    if (!length(f)) return(NULL)
    b <- basename(f)
    dest <- if (!tg$has_g[i] || tg$pat[i] == tg$dest_pat[i]) {
      b
    } else if (!grepl("[*?]", tg$raw[i])) {
      want <- paste0(prefix, "_", tg$pat[i], ".", ext)
      if (!identical(b, want)) stop(dir, ": 패턴 ", tg$raw[i], "에 예상 밖 파일 ", paste(b, collapse = ", "))
      paste0(prefix, "_", tg$dest_pat[i], ".", ext)
    } else {
      tpl <- paste0(prefix, "_", tg$raw[i], ".", ext)
      vapply(b, glob_dest_name, "", template = tpl, gene = tg$gene[i], name = tg$name[i])
    }
    data.frame(source = f, dest_base = unname(dest), pattern = tg$pat[i], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame(source = character(), dest_base = character(), pattern = character()) else out
}

# 대상 경로 중복 (대소문자 무시: Windows 파일 시스템) → 원본 파일을 나열하고 중단
check_plan_collisions <- function(plan) {
  key <- tolower(plan$rel)
  dup <- unique(key[duplicated(key)])
  if (length(dup)) {
    msg <- vapply(dup, function(k) paste0(plan$rel[key == k][1], " ← ", paste(plan$source[key == k], collapse = ", ")), "")
    stop("대상 경로 중복 (서로 다른 원본이 같은 이름으로 내보내짐):\n  ", paste(msg, collapse = "\n  "),
         if (isTRUE(cfg$export$use_custom_filename)) "\n  → cfg$export$custom_filename 또는 cfg$export$layout 확인" else "",
         call. = FALSE)
  }
  invisible(plan)
}

# ---- 복사 계획 / 실행 ---------------------------------------------------

# 내용 md5: docx는 zip 안의 항목별 md5 (생성 시각이 들어가는 docProps/core.xml 제외)
#           → 같은 내용을 다시 만들어도 "same"
content_md5 <- function(path) {
  if (!grepl("\\.docx$", path, ignore.case = TRUE)) return(unname(tools::md5sum(path)))
  tmp <- tempfile("docx_")
  on.exit(unlink(tmp, recursive = TRUE))
  utils::unzip(path, exdir = tmp)
  f <- sort(list.files(tmp, recursive = TRUE))
  f <- f[f != "docProps/core.xml"]
  h <- unname(tools::md5sum(file.path(tmp, f)))
  digest::digest(paste(f, h, collapse = "\n"), algo = "md5", serialize = FALSE)
}

# 원본 → 대상 상대경로 목록에 상태 부여: new / same / replace
plan_copy <- function(plan, dest_root) {
  if (!nrow(plan)) return(transform(plan, dest = character(), status = character(), md5 = character()))
  plan$dest   <- file.path(dest_root, plan$rel)
  plan$md5    <- vapply(plan$source, content_md5, "")
  exists      <- file.exists(plan$dest)
  dest_md5    <- ifelse(exists, vapply(plan$dest, function(d) if (file.exists(d)) content_md5(d) else "", ""), "")
  plan$status <- ifelse(!exists, "new", ifelse(dest_md5 == plan$md5, "same", "replace"))
  plan
}

# new/replace 실행. replace는 기존 파일을 _old/<ts>/<rel>로 이동 후 복사 (삭제 없음)
# 반환: 실패 메시지 벡터
apply_copy <- function(plan, dest_root, old_dir) {
  fail <- character()
  for (i in which(plan$status %in% c("new", "replace"))) {
    dest <- plan$dest[i]
    ok <- tryCatch({
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      if (plan$status[i] == "replace") {
        bak <- file.path(old_dir, plan$rel[i])
        k <- 1
        while (file.exists(bak)) {   # 같은 분에 두 번 실행한 경우: 덮어쓰지 않고 번호 붙임
          bak <- file.path(old_dir, paste0(plan$rel[i], ".", k))
          k <- k + 1
        }
        dir.create(dirname(bak), recursive = TRUE, showWarnings = FALSE)
        if (!file.rename(dest, bak)) stop("기존 파일 이동 실패 (열려 있는지 확인)")
      }
      if (!file.copy(plan$source[i], dest, overwrite = FALSE, copy.date = TRUE)) stop("복사 실패")
      if (content_md5(dest) != plan$md5[i]) stop("복사 후 md5 불일치")
      TRUE
    }, error = function(e) {
      fail <<- c(fail, paste0(plan$rel[i], ": ", conditionMessage(e)))
      FALSE
    })
    if (!ok) plan$status[i] <- "failed"
  }
  attr(fail, "plan") <- plan
  fail
}

# ---- GEO 결과가 이 유전자의 것인지 ----------------------------------------
# 1) 06 survival_summary_raw.csv의 gene 열 → 2) 없으면 probe_QC.csv의 symbol → 3) 둘 다 없으면 NA
# 반환: list(gene = 확인된 유전자 또는 NA, source = 근거)
geo_result_gene <- function(gse) {
  dir <- src_dir("GEO", "tables", gse)
  raw <- file.path(dir, out_file(gse, "survival_summary_raw.csv"))
  if (file.exists(raw)) {
    r <- read.csv(raw, stringsAsFactors = FALSE, check.names = FALSE)
    if ("gene" %in% names(r)) {
      g <- unique(toupper(r$gene))
      if (length(g) != 1) stop(raw, ": gene 열 값이 하나가 아님: ", paste(g, collapse = ", "))
      return(list(gene = g, source = "survival_summary_raw$gene"))
    }
  }
  pq <- file.path(dir, out_file(gse, "probe_QC.csv"))
  if (file.exists(pq)) {
    p <- read.csv(pq, stringsAsFactors = FALSE)
    if (all(c("symbol", "chosen") %in% names(p))) {
      sym <- unique(trimws(unlist(strsplit(p$symbol[p$chosen %in% TRUE], "///", fixed = TRUE))))
      if (length(sym) == 1) return(list(gene = toupper(sym), source = "probe_QC$symbol"))
    }
  }
  list(gene = NA_character_, source = "none")
}

# ---- 핵심 수치 한 장 (<G>_key_results.csv) -------------------------------

key_cols <- c("cohort", "endpoint", "population", "model", "gene_term", "n", "events",
              "HR", "lo", "hi", "p", "logrank_p", "rmst_diff", "ph_flag", "exploratory",
              "gsea_n_sig", "gsea_top_up", "gsea_top_down")

key_row <- function(...) {
  r <- list(...)
  out <- as.data.frame(setNames(rep(list(NA), length(key_cols)), key_cols))
  for (k in names(r)) out[[k]] <- r[[k]] %||% NA
  out
}

# TCGA: 04 survival_summary_raw.csv의 해당 유전자 행 → 모형별 행
key_rows_tcga <- function(cancer, gene) {
  f <- file.path(src_dir(cancer, "tables"), out_file(cancer, "survival_summary_raw.csv"))
  s <- read.csv(f, stringsAsFactors = FALSE)
  s <- s[s$gene == gene, ]
  if (nrow(s) != 1) stop(f, ": ", gene, " 행이 1개가 아님 (", nrow(s), ")")
  gt <- s$cox_gene_term %||% "group"
  base <- function(model, n, ev, hr, lo, hi, p, term = gt, expl = NA) {
    key_row(cohort = cancer, endpoint = "OS", population = "All", model = model, gene_term = term,
            n = n, events = ev, HR = hr, lo = lo, hi = hi, p = p, logrank_p = s$logrank_p,
            rmst_diff = s$rmst_diff, ph_flag = s$ph_flag, exploratory = expl)
  }
  out <- rbind(
    base("Univariable", s$n, s$events, s$uni_hr, s$uni_lo, s$uni_hi, s$uni_p),
    base("Multivariable", s$multi_n, s$multi_events, s$multi_hr, s$multi_lo, s$multi_hi, s$multi_p,
         expl = s$exploratory),
    base(paste0("Multivariable + strata(", cfg$strata_var, ")"), s$multi_n, s$multi_events,
         s$strata_hr, s$strata_lo, s$strata_hi, s$strata_p, expl = s$exploratory))
  if (!is.na(s$hr_early)) {
    sc <- cfg$ph_split_months
    out <- rbind(out,
      base(sprintf("Multivariable, 0–%d months", sc), s$multi_n, NA, s$hr_early, s$lo_early, s$hi_early, s$p_early),
      base(sprintf("Multivariable, >%d months", sc), NA, NA, s$hr_late, s$lo_late, s$hi_late, s$p_late))
  }
  out
}

# GEO: 06 survival_summary_raw.csv 전체
key_rows_geo <- function(gse) {
  f <- file.path(src_dir("GEO", "tables", gse), out_file(gse, "survival_summary_raw.csv"))
  s <- read.csv(f, stringsAsFactors = FALSE)
  term <- ifelse(s$term == "expr_log2", "continuous", s$term)
  do.call(rbind, lapply(seq_len(nrow(s)), function(i) {
    key_row(cohort = gse, endpoint = s$endpoint[i], population = s$population[i], model = s$model[i],
            gene_term = term[i], n = s$n[i], events = s$events[i], HR = s$hr[i], lo = s$lo[i],
            hi = s$hi[i], p = s$p[i], logrank_p = s$logrank_p[i], rmst_diff = s$rmst_diff[i],
            ph_flag = !is.na(s$ph_p_term[i]) && s$ph_p_term[i] < 0.05, exploratory = s$exploratory[i])
  }))
}

# GSEA: 컬렉션마다 유의 경로 수(padj < cfg$export$gsea_fdr) + NES 상위 경로 (양/음)
key_rows_gsea <- function(cancer, gene) {
  dir <- src_dir(cancer, "tables", file.path("gsea", gene))
  do.call(rbind, lapply(names(cfg$gsea_collections), function(cl) {
    f <- file.path(dir, out_file(cancer, paste0("gsea_", cl, ".csv")))
    if (!file.exists(f)) return(NULL)
    g <- read.csv(f, stringsAsFactors = FALSE)
    sig <- g[!is.na(g$padj) & g$padj < cfg$export$gsea_fdr, ]
    top <- function(x) paste(head(x$pathway, cfg$export$gsea_top), collapse = "; ")
    key_row(cohort = cancer, model = paste("GSEA", cl), gsea_n_sig = nrow(sig),
            gsea_top_up = top(sig[sig$NES > 0, ][order(-sig$NES[sig$NES > 0]), ]),
            gsea_top_down = top(sig[sig$NES < 0, ][order(sig$NES[sig$NES < 0]), ]))
  }))
}

# ---- git 상태 (manifest용) ---------------------------------------------------
git_info <- function() {
  run <- function(...) tryCatch(suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = FALSE)),
                                error = function(e) NA_character_)
  head <- run("rev-parse", "HEAD")
  st   <- run("status", "--porcelain")
  list(head = if (length(head)) head[1] else NA_character_,
       dirty = length(st) > 0 && !all(is.na(st)))
}

# ---- <G>_tables.docx: 양식(cfg$export$template)을 파이프라인 출력으로 채우기 ----------------
# 양식은 읽기만 함. 캡션 문장·글꼴·크기는 양식에서 읽고, 행 구성은 파이프라인 CSV 기준.
# 양식 캡션 순서 = Table 1, 상관 표 × 암종, KM 그림 × 암종, Cox 표 × 암종 (cfg$cancers 순서)

# 양식 읽기: 캡션 (Table/Figure로 시작하는 문단), 표 글꼴/크기, Table 1의 변수 이름
read_template <- function(path) {
  x    <- xml2::read_xml(unz(path, "word/document.xml"))
  txt  <- function(node) paste(xml2::xml_text(xml2::xml_find_all(node, ".//w:t")), collapse = "")
  paras <- vapply(xml2::xml_find_all(x, "/w:document/w:body/w:p"), txt, "")
  caps <- trimws(paras[grepl("^(Table|Figure) [0-9]+\\.", paras)])
  font <- xml2::xml_attr(xml2::xml_find_all(x, "//w:tbl//w:rFonts"), "ascii")
  font <- font[!is.na(font)]
  sz   <- as.numeric(xml2::xml_attr(xml2::xml_find_all(x, "//w:tbl//w:sz"), "val"))
  sz   <- sz[!is.na(sz)]
  # Table 1의 변수 행 = 첫 칸에만 글자가 있는 행 (칸 안 여러 문단은 공백으로 연결)
  cell_txt <- function(tc) trimws(paste(vapply(xml2::xml_find_all(tc, ".//w:p"), txt, ""), collapse = " "))
  rows <- xml2::xml_find_all(xml2::xml_find_first(x, "//w:tbl"), "./w:tr")
  vars <- unlist(lapply(rows, function(tr) {
    v <- vapply(xml2::xml_find_all(tr, "./w:tc"), cell_txt, "")
    if (nzchar(v[1]) && all(!nzchar(v[-1]))) v[1] else NULL
  }))
  list(captions = caps,
       font = if (length(font)) names(which.max(table(font))) else cfg$export$font_name,
       size = if (length(sz)) as.numeric(names(which.max(table(sz)))) / 2 else cfg$export$font_size,  # half-point → pt
       table1_vars = vars)
}

# 캡션: 양식 유전자 → gene, cfg$export$caption_fixes 적용. 반환 list(text, fixed = 고친 항목)
fill_caption <- function(cap, tpl_gene, gene) {
  fixed <- character()
  for (k in names(cfg$export$caption_fixes)) {
    if (grepl(k, cap, fixed = TRUE)) {
      cap <- gsub(k, cfg$export$caption_fixes[[k]], cap, fixed = TRUE)
      fixed <- c(fixed, paste0("\"", k, "\" → \"", cfg$export$caption_fixes[[k]], "\""))
    }
  }
  list(text = gsub(paste0("\\b", tpl_gene, "\\b"), gene, cap), fixed = fixed)
}

# 3선 표 (양식: 머리 위·아래, 본문 아래 0.5 pt 실선)
style_table <- function(ft, tpl, bold_rows = NULL) {
  b <- officer::fp_border(color = "black", width = 0.5)
  ft <- ft |>
    flextable::border_remove() |>
    flextable::hline_top(border = b, part = "header") |>
    flextable::hline_bottom(border = b, part = "header") |>
    flextable::hline_bottom(border = b, part = "body") |>
    flextable::font(fontname = tpl$font, part = "all") |>
    flextable::fontsize(size = tpl$size, part = "all") |>
    flextable::bold(part = "header") |>
    flextable::align(j = -1, align = "center", part = "all") |>
    flextable::padding(padding.top = 1, padding.bottom = 1, part = "all")
  if (length(bold_rows)) ft <- flextable::bold(ft, i = bold_rows, j = 1, part = "body")
  flextable::autofit(ft)
}

# 열 이름이 "<수준>, N = …"인 그룹 열 (03 table1 CSV)
group_col <- function(nm, level) {
  j <- which(startsWith(nm, paste0(level, ",")))
  if (length(j) != 1) stop("table1 CSV에서 ", level, " 열을 찾지 못함: ", paste(nm, collapse = " | "))
  j
}

# Table 1: 암종별 table1 CSV → 변수/수준 행을 합치고 [High | Low | P-value] × 암종
table1_block <- function(gene) {
  lv <- if (isTRUE(cfg$export$table1_alt_first)) rev(cfg$group_levels) else cfg$group_levels
  sep <- "\u001f"
  per <- lapply(cfg$cancers, function(cancer) {
    f <- file.path(src_dir(cancer, "tables"), out_file(cancer, paste0("table1_", gene, ".csv")))
    t <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, na.strings = "", encoding = "UTF-8", colClasses = "character")
    if (!"Test" %in% names(t)) stop(f, ": Test 열 없음 → 03_table1.R 재실행")
    is_var <- !is.na(t$Test) & nzchar(t$Test)
    vname <- t[[1]][is_var][cumsum(is_var)]
    cnt <- function(x) if (isTRUE(cfg$export$table1_counts_only)) sub("\\s*\\(.*\\)$", "", x) else x
    d <- data.frame(key = paste(vname, ifelse(is_var, "", t[[1]]), sep = sep), stringsAsFactors = FALSE)
    for (l in lv) {
      x <- t[[group_col(names(t), l)]]
      d[[l]] <- ifelse(is.na(x), "", ifelse(is_var, x, cnt(x)))   # 연속형 변수 행은 median [IQR] 그대로
    }
    d$p <- ifelse(is_var & !is.na(t[["p-value"]]), t[["p-value"]], "")
    d
  })
  # 변수 순서 = 등장 순서 합집합, 수준도 변수 안에서 등장 순서 합집합
  keys <- unique(unlist(lapply(per, `[[`, "key")))
  kv <- sub(paste0(sep, ".*$"), "", keys)
  kl <- sub(paste0("^.*", sep), "", keys)
  vars <- unique(kv)
  o <- order(match(kv, vars), seq_along(keys))
  keys <- keys[o]; kv <- kv[o]; kl <- kl[o]
  out <- data.frame(Characteristic = ifelse(kl == "", kv, paste0("  ", kl)), stringsAsFactors = FALSE)
  for (k in seq_along(per)) {
    m <- match(keys, per[[k]]$key)
    for (l in lv) out[[paste(cfg$cancers[k], l, sep = "_")]] <- ifelse(is.na(m), "-", per[[k]][[l]][m])
    out[[paste(cfg$cancers[k], "p", sep = "_")]] <- ifelse(is.na(m), "", per[[k]]$p[m])
  }
  list(df = out, var_rows = which(kl == ""), vars = vars, levels = lv)
}

table1_flex <- function(gene, tpl) {
  b <- table1_block(gene)
  ft <- flextable::flextable(b$df)
  hdr <- c(Characteristic = "", setNames(rep(c(b$levels, "P-value"), length(cfg$cancers)), names(b$df)[-1]))
  ft <- flextable::set_header_labels(ft, values = as.list(hdr))
  ft <- flextable::add_header_row(ft, values = c("", cfg$cancers),
                                  colwidths = c(1, rep(length(b$levels) + 1, length(cfg$cancers))))
  list(ft = style_table(ft, tpl, b$var_rows), vars = b$vars)
}

# 상관 표: correlation CSV (long) → 변수 × (R, P) 행렬
cor_flex <- function(cancer, gene, tpl) {
  f <- file.path(src_dir(cancer, "tables"), out_file(cancer, paste0("correlation_", gene, ".csv")))
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  v <- unique(d$var1)
  get <- function(a, b, col) d[[col]][d$var1 == a & d$var2 == b]
  rows <- lapply(v, function(a) {
    r <- vapply(v, function(b) if (a == b) "1" else sprintf("%.3f", get(a, b, "r")), "")
    p <- vapply(v, function(b) if (a == b) "" else fmt_p(get(a, b, "p")), "")
    rbind(c(a, "R", r), c("", "P", p))
  })
  m <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(m) <- c("var", "stat", v)
  ft <- flextable::flextable(m) |> flextable::set_header_labels(values = list(var = "", stat = ""))
  style_table(ft, tpl, which(m$var != ""))
}

# cox CSV (04 format_cox) → 변수·수준별 행: label "<Variable> (<수준> vs <기준>)", hr, p
cox_items <- function(f) {
  t <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, na.strings = "", encoding = "UTF-8", colClasses = "character")
  t$Variable <- Reduce(function(a, b) if (is.na(b)) a else b, t$Variable, accumulate = TRUE)
  t$Level[is.na(t$Level)] <- ""
  is_ref <- grepl(" \\(ref\\)$", t$Level)
  out <- list()
  for (v in unique(t$Variable)) {
    i <- which(t$Variable == v)
    ref <- sub(" \\(ref\\)$", "", t$Level[i][is_ref[i]])
    for (k in i[!is_ref[i]]) {
      lab <- if (length(ref)) paste0(v, " (", t$Level[k], " vs ", ref, ")") else v
      out[[length(out) + 1]] <- data.frame(var = v, label = lab, hr = t[["HR (95% CI)"]][k], p = t$p[k],
                                           stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

# Cox 표: 다변량 모형의 항 순서 (유전자 행은 마지막), 단변량은 cox_uni의 같은 변수·같은 척도 행
cox_flex <- function(cancer, gene, tpl) {
  dir <- src_dir(cancer, "tables")
  mv <- cox_items(file.path(dir, out_file(cancer, paste0("cox_multi_", gene, ".csv"))))
  uv <- cox_items(file.path(dir, out_file(cancer, "cox_uni.csv")))
  is_gene <- startsWith(mv$var, paste0(gene, " expression"))
  mv <- rbind(mv[!is_gene, ], mv[is_gene, ])
  m <- match(mv$label, uv$label)
  if (anyNA(m)) warning(cancer, ": 단변량 표(cox_uni)에 없는 항: ", paste(mv$label[is.na(m)], collapse = ", "), call. = FALSE)
  df <- data.frame(Variable = mv$label,
                   uh = ifelse(is.na(m), "–", uv$hr[m]), up = ifelse(is.na(m), "–", uv$p[m]),
                   mh = mv$hr, mp = mv$p, stringsAsFactors = FALSE)
  df[is.na(df)] <- ""
  ft <- flextable::flextable(df) |>
    flextable::set_header_labels(values = list(Variable = "Variable", uh = "Univariate HR (95% CI)", up = "P value",
                                               mh = "Multivariate HR (95% CI)", mp = "P value"))
  s <- read.csv(file.path(dir, out_file(cancer, "survival_summary_raw.csv")), stringsAsFactors = FALSE)
  s <- s[s$gene == gene, ]
  scale <- if ((s$cox_gene_term %||% "group") == "continuous") {
    paste0(gene, " HR per 1 log2 unit of expression.")
  } else {
    paste0(gene, " HR ", group_contrast_label(), ".")
  }
  note <- sprintf("Multivariate model: n = %d, events = %d, EPV = %.1f%s. %s", s$multi_n, s$multi_events, s$epv,
                  if (isTRUE(s$exploratory)) paste0(" (exploratory, EPV < ", cfg$epv_min, ")") else "", scale)
  list(ft = style_table(ft, tpl), note = note)
}

# TIFF 크기 (픽셀) — magick 없이 첫 IFD의 ImageWidth(256)/ImageLength(257)만 읽음
tiff_dims <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  end <- if (rawToChar(readBin(con, "raw", 2)) == "II") "little" else "big"
  rd <- function(size) {
    v <- readBin(con, "integer", 1, size = size, endian = end, signed = size >= 4)
    if (size == 4 && v < 0) v + 2^32 else v
  }
  rd(2)
  seek(con, rd(4))
  dims <- c()
  for (i in seq_len(rd(2))) {
    tag <- rd(2); type <- rd(2); rd(4)
    pos <- seek(con)
    if (tag %in% c(256, 257)) dims[as.character(tag)] <- if (type == 3) rd(2) else rd(4)
    seek(con, pos + 4)
  }
  if (length(dims) != 2) stop(path, ": TIFF 크기를 읽지 못함")
  list(w = dims[["256"]], h = dims[["257"]])
}

build_tables_docx <- function(gene, out_path) {
  tpl <- read_template(cfg$export$template)
  nc <- length(cfg$cancers)
  if (length(tpl$captions) != 1 + 3 * nc) {
    stop("양식 캡션 ", length(tpl$captions), "개 — Table 1 + 암종(", nc, ")마다 상관 표·KM 그림·Cox 표가 필요 (",
         cfg$export$template, ")")
  }
  tpl_gene <- cfg$export$template_gene %||% sub("^.* (\\S+) expression.*$", "\\1", tpl$captions[1])
  caps <- lapply(tpl$captions, fill_caption, tpl_gene = tpl_gene, gene = gene)
  for (k in seq_len(nc)) for (i in 1 + k + nc * (0:2)) {
    if (!grepl(cfg$cancers[k], caps[[i]]$text, fixed = TRUE)) {
      stop("양식 캡션 순서가 cfg$cancers와 다름: \"", caps[[i]]$text, "\"에 ", cfg$cancers[k], " 없음")
    }
  }
  par_txt <- function(txt, size) {
    officer::fpar(officer::ftext(txt, officer::fp_text(font.family = tpl$font, font.size = size)))
  }

  doc <- officer::read_docx()
  t1 <- table1_flex(gene, tpl)
  doc <- officer::body_add_fpar(doc, par_txt(caps[[1]]$text, tpl$size + 1))
  doc <- flextable::body_add_flextable(doc, t1$ft)
  for (k in seq_len(nc)) {
    doc <- officer::body_add_par(doc, "")
    doc <- officer::body_add_fpar(doc, par_txt(caps[[1 + k]]$text, tpl$size + 1))
    doc <- flextable::body_add_flextable(doc, cor_flex(cfg$cancers[k], gene, tpl))
  }
  for (k in seq_len(nc)) {
    cancer <- cfg$cancers[k]
    img <- file.path(src_dir(cancer, "figures"), out_file(cancer, paste0("km_", gene, ".", cfg$export$fig_ext)))
    td <- tiff_dims(img)
    doc <- officer::body_add_par(doc, "")
    doc <- officer::body_add_img(doc, img, width = cfg$export$fig_width, height = cfg$export$fig_width * td$h / td$w)
    doc <- officer::body_add_fpar(doc, par_txt(caps[[1 + nc + k]]$text, tpl$size + 1))
  }
  for (k in seq_len(nc)) {
    cx <- cox_flex(cfg$cancers[k], gene, tpl)
    doc <- officer::body_add_par(doc, "")
    doc <- officer::body_add_fpar(doc, par_txt(caps[[1 + 2 * nc + k]]$text, tpl$size + 1))
    doc <- flextable::body_add_flextable(doc, cx$ft)
    doc <- officer::body_add_fpar(doc, par_txt(cx$note, tpl$size - 1))
  }
  print(doc, target = out_path)

  # 보고: 고친 양식 오류, 양식에만 있는 Table 1 변수
  fixes <- unique(unlist(lapply(caps, `[[`, "fixed")))
  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  # 출력 라벨·변수 이름과 앞부분이 같거나 cfg$export$template_aliases에 대응이 있으면 "있음"
  nv <- norm(c(t1$vars, intersect(cfg$table1_vars, names(var_labels))))
  tv <- tpl$table1_vars
  al <- cfg$export$template_aliases
  tv_n <- norm(ifelse(tv %in% names(al), al[tv], tv))
  only_tpl <- tv[!vapply(tv_n, function(v) any(startsWith(nv, v) | startsWith(v, nv)), TRUE)]
  cat("\n[", basename(out_path), "] 양식 유전자 ", tpl_gene, " → ", gene, "\n", sep = "")
  if (length(fixes)) cat("  양식 캡션 오타 수정:", paste(fixes, collapse = ", "), "\n")
  cat("  Table 1·Cox 표의 수준 라벨은 파이프라인 CSV 기준 (CEA ≤5 / >5 ng/mL; 양식의 CEA ≤65/>65 행과 \">=5 vs <5\" 라벨은 쓰지 않음)\n")
  if (length(only_tpl)) {
    cat("  양식에 있으나 출력에 없는 변수 (이름이 달라 대응되는 변수가 있을 수 있음):", paste(only_tpl, collapse = ", "), "\n")
  }
  invisible(out_path)
}
