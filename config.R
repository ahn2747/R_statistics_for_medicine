# =============================================================
# config.R
# 프로젝트 설정 — 암종/유전자/분석 옵션은 여기서만 수정
# 다른 TCGA 암종 추가: cancers에 코드 추가 + database/<CANCER>.sav (없으면 GDC 임상으로 대체)
# =============================================================

cfg <- list(
  cancers      = c("COAD", "READ"),        # 예: c("COAD", "READ", "LIHC")
  clinical_dir = "database",               # <CANCER>.sav 위치
  gene_dir     = "database/gene_files",    # <CANCER>_<split>_<GENE>.csv 위치
  gene_pattern = "^([A-Za-z]+)_(\\d+_\\d+)_(.+)\\.csv$",   # 암종, split, 유전자
  genes_from_gdc = c(),                    # 02a 발현 행렬에서 추가로 가져올 유전자 (예: c("CD8A"))

  # 생존 엔드포인트 (clean_names 이후 이름)
  surv_time   = "days",
  surv_status = "status",                  # 1 = 사망

  drop_columns = c("rock2g_01"),           # 중복 열 (ROCK2Group과 동일)

  # Table 1 변수 (존재하는 열만 사용)
  table1_vars = c("age", "age_g", "gender", "pathologic_stage",
                  "pathologic_t", "pathologic_n", "pathologic_m",
                  "pathologic_stage_12_34", "cea_g",
                  "lymphovascular_invasion_indicator", "vascular_invasion_indicator",
                  "residual_tumor", "histologic_diagnosis",
                  "anatomic_neoplasm_subdivision", "tumor_status",
                  "kras_gene_analysis_indicator", "braf_gene_analysis_indicator",
                  "mismatch_rep_proteins_tested_by_ihc"),

  exclude_neoadjuvant  = TRUE,            # history_neoadjuvant_treatment == "Yes" 제외
  cox_covariates       = c("age", "gender", "stage"),   # age는 연속형
  epv_min              = 10,               # 다변량 Cox EPV 경고 기준
  median_tie_tolerance = 1,                # median-split 재현 QC 허용 오차 (명)

  # 01 단계 생존 n/사건 수 검증 (NULL이면 검증 생략)
  expected_surv = list(
    COAD = list(excl = c(n = 436, ev = 96), incl = c(n = 439, ev = 97)),
    READ = list(excl = c(n = 156, ev = 25), incl = c(n = 157, ev = 25))
  ),

  gdc_dir       = "D:/GDCdata",            # 02a 다운로드 폴더 (Windows 경로 길이 문제 회피)
  processed_dir = "data/processed",
  output_dir    = "output"
)
