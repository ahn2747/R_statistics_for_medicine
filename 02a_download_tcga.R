# =============================================================
# 02a_download_tcga.R
# GDC에서 TCGA RNA-seq(STAR - Counts) 다운로드 + 전처리 (암종은 config.R의 cfg$cancers)
# 결과:
#   database/TCGA_<COHORT>_RNAseq_Expression.csv  (log2(TPM+1), 유전자 x 환자)
#   data/processed/<COHORT>_counts.rds            (raw counts, GSEA/DESeq2용)
# =============================================================

# ---- 0. 패키지 (최초 1회만) ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("TCGAbiolinks", "SummarizedExperiment")) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE)
}
library(TCGAbiolinks)
library(SummarizedExperiment)

source("config.R")

dir.create("database", showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Windows 경로 길이 문제를 피하려고 짧은 폴더에 받음
gdc_dir <- cfg$gdc_dir

process_cohort <- function(cohort) {
  project <- paste0("TCGA-", cohort)
  message("==== ", project, " ====")

  # ---- 1. 조회 ----
  query <- GDCquery(
    project       = project,
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  # ---- 2. 다운로드 (중간에 끊기면 다시 실행하면 이어서 받음) ----
  GDCdownload(query, method = "api", directory = gdc_dir, files.per.chunk = 20)

  # ---- 3. SummarizedExperiment로 합치기 ----
  se <- GDCprepare(query, directory = gdc_dir,
                   save = TRUE,
                   save.filename = file.path("data/processed", paste0(cohort, "_se.rda")))

  # ---- 4. 종양 샘플(01)만 ----
  bc <- colnames(se)
  se <- se[, substr(bc, 14, 15) == "01"]

  # ---- 5. 환자 단위(12자리)로 중복 정리: 첫 번째 aliquot 유지 ----
  pid <- substr(colnames(se), 1, 12)
  dup <- duplicated(pid)
  if (any(dup)) {
    message("중복 환자 ", sum(dup), "명 → 첫 aliquot만 유지: ",
            paste(unique(pid[dup]), collapse = ", "))
  }
  se <- se[, !dup]
  colnames(se) <- substr(colnames(se), 1, 12)

  # ---- 6. 유전자 ID → Symbol, 중복 Symbol은 평균 발현이 가장 높은 것만 ----
  sym  <- rowData(se)$gene_name
  tpm  <- assay(se, "tpm_unstrand")
  keep <- !is.na(sym) & sym != ""
  tpm  <- tpm[keep, ]; sym <- sym[keep]
  ord  <- order(rowMeans(tpm), decreasing = TRUE)
  first <- ord[!duplicated(sym[ord])]
  tpm  <- tpm[first, ]
  rownames(tpm) <- sym[first]

  cnt <- assay(se, "unstranded")[keep, ][first, ]
  rownames(cnt) <- sym[first]

  # ---- 7. 저장 ----
  expr <- log2(tpm + 1)
  out_csv <- file.path("database", paste0("TCGA_", cohort, "_RNAseq_Expression.csv"))
  write.csv(data.frame(Gene = rownames(expr), expr, check.names = FALSE),
            out_csv, row.names = FALSE)
  saveRDS(cnt, file.path("data/processed", paste0(cohort, "_counts.rds")))

  message(cohort, ": ", nrow(expr), " genes x ", ncol(expr), " patients → ", out_csv)
  invisible(expr)
}

for (cohort in cfg$cancers) process_cohort(cohort)

# ---- 8. 기존 .sav와 겹치는 환자 수 확인 ----
library(haven)
for (cohort in cfg$cancers) {
  sav <- file.path("database", paste0(cohort, ".sav"))
  if (!file.exists(sav)) next
  clin <- read_sav(sav)
  expr <- read.csv(file.path("database", paste0("TCGA_", cohort, "_RNAseq_Expression.csv")),
                   check.names = FALSE, nrows = 1)
  ids  <- trimws(clin$sampleID)
  cat(cohort, ": 임상", length(ids), "명 중 발현 매칭",
      sum(ids %in% colnames(expr)), "명\n")
}
