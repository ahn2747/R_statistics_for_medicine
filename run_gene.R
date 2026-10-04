# =============================================================
# run_gene.R <GENE> [--skip-download] [--from=<step>]
# 유전자 하나: 다운로드 → 02b → 03 → 04 → 05 (암종마다) → 06 → 90_export
#   각 단계는 별도 프로세스, 종료 코드가 0이 아니면 즉시 중단
#   자식 프로세스에 MEDIN_GENE=<GENE> → config.R의 cfg$primary_gene
#   단계 이름: download, 02b, 03, 04, 05, 06, 90  (예: --from=04)
#   로그: output/run_logs/<GENE>_<타임스탬프>.log (단계별 시작/종료/소요 시간)
# 커밋은 하지 않음: 끝에 바뀐 output 파일 목록만 출력
# =============================================================

source("config.R")
source("R/utils.R")

args <- commandArgs(trailingOnly = TRUE)
pos  <- args[!startsWith(args, "--")]
if (length(pos) != 1) stop("사용법: Rscript run_gene.R <GENE> [--skip-download] [--from=<step>]")
gene    <- toupper(pos[1])
skip_dl <- "--skip-download" %in% args
from    <- sub("^--from=", "", args[startsWith(args, "--from=")])

rs <- function(script, ...) list(cmd = cfg$rscript, args = c(script, ...))
steps <- c(
  list(download = list(cmd = cfg$python, args = c("csv_download.py", gene, "--cancers", cfg$cancers)),
       `02b` = rs("02b_merge_genes.R"),
       `03`  = rs("03_table1.R"),
       `04`  = rs("04_survival.R")),
  setNames(lapply(cfg$cancers, function(cancer) rs("05_gsea.R", cancer)), paste0("05_", cfg$cancers)),
  list(`06` = rs("06_external_geo.R"),
       `90` = rs("90_export.R", gene))
)
step_group <- sub("_.*$", "", names(steps))   # 05_COAD → 05
if (length(from)) {
  if (!from %in% step_group) stop("--from=", from, " 없음 (가능: ", paste(unique(step_group), collapse = ", "), ")")
  steps <- steps[seq_along(steps) >= match(from, step_group)]
}
if (skip_dl) steps <- steps[names(steps) != "download"]

# --skip-download / --from 으로 다운로드를 건너뛰면 gene CSV가 있어야 함
if (!"download" %in% names(steps)) {
  f <- list.files(cfg$gene_dir, pattern = cfg$gene_pattern)
  m <- str_match(f, cfg$gene_pattern)
  have <- m[toupper(m[, 4]) == gene, 2]
  miss <- setdiff(cfg$cancers, have)
  if (length(miss)) {
    stop(gene, " gene CSV 없음 (", paste(miss, collapse = ", "), ", ", cfg$gene_dir,
         ") → --skip-download 없이 실행하거나 py csv_download.py ", gene)
  }
}

log_dir <- file.path(cfg$output_dir, "run_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_f <- file.path(log_dir, paste0(gene, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
logw <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = log_f, append = TRUE)
}

cat("\n주의:\n",
    "  - primary gene을 ", gene, "(으)로 바꾸면 04의 BH 보정 대상이 바뀌어 다른 유전자의 q값도 바뀜.\n",
    "    끝에 바뀐 output 파일 목록을 보여 주며, 커밋은 직접 결정 (자동 커밋 없음).\n",
    "  - 05·06은 primary gene만 실행 → 이전 유전자의 GSEA/GEO 결과는 _stale/로 이동되거나 덮어써짐.\n",
    "    유전자마다 90_export까지 끝낸 뒤 다음 유전자로 넘어갈 것.\n\n", sep = "")
logw("run_gene.R ", gene, " 시작: 단계 ", paste(names(steps), collapse = " → "),
     if (skip_dl) " (--skip-download)" else "", if (length(from)) paste0(" (--from=", from, ")") else "")

Sys.setenv(MEDIN_GENE = gene, PYTHONUTF8 = "1")
t_all <- Sys.time()
for (nm in names(steps)) {
  s <- steps[[nm]]
  logw("[", nm, "] 시작: ", basename(s$cmd), " ", paste(s$args, collapse = " "))
  t0 <- Sys.time()
  status <- system2(s$cmd, s$args)
  dt <- format(round(difftime(Sys.time(), t0, units = "mins"), 1))
  if (!identical(as.integer(status), 0L)) {
    logw("[", nm, "] 실패 (종료 코드 ", status, ", ", dt, ") → 중단. 고친 뒤 --from=", sub("_.*$", "", nm), "로 재시작")
    quit(save = "no", status = 1)
  }
  logw("[", nm, "] 완료 (", dt, ")")
}
logw("전체 완료 (", format(round(difftime(Sys.time(), t_all, units = "mins"), 1)), ")")

# 바뀐 output 파일 (커밋은 사용자가 결정)
st <- tryCatch(system2("git", c("status", "--porcelain", "--", cfg$output_dir), stdout = TRUE), error = function(e) NULL)
cat("\ngit status (", cfg$output_dir, "): ", length(st), "개 변경\n", sep = "")
if (length(st)) cat(paste0("  ", st, collapse = "\n"), "\n")
cat("로그:", log_f, "\n")
