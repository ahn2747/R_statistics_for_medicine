# =============================================================
# config.R
# 프로젝트 설정 — 암종/유전자/분석 옵션은 여기서만 수정
# 새 암종 추가 / .sav 갱신 절차는 CLAUDE.md의 체크리스트 참고
# =============================================================

cfg <- list(
  cancers      = c("COAD", "READ"),        # 예: c("COAD", "READ", "LIHC")
  clinical_dir = "database",               # <CANCER>.sav 위치
  gene_dir     = "database/gene_files",    # <CANCER>_<split>_<GENE>.csv 위치
  gene_pattern = "^([A-Za-z]+)_(\\d+_\\d+)_(.+)\\.csv$",   # 암종, split, 유전자
  genes_from_gdc = c(),                    # 02a 발현 행렬에서 추가로 가져올 유전자 (예: c("CD8A"))
  primary_gene = "MS4A1",                  # 주 가설 유전자: BH 보정에서 제외, raw p 보고

  drop_columns = c("rock2g_01"),           # 중복 열 (ROCK2Group과 동일)

  # Table 1 변수 (존재하는 열만 사용)
  table1_vars = c("age", "age_g", "gender", "pathologic_stage",
                  "pathologic_t", "pathologic_n", "pathologic_m",
                  "pathologic_stage_12_34", "cea_g",
                  "lymphovascular_invasion_indicator", "vascular_invasion_indicator",
                  "residual_tumor", "histologic_diagnosis",
                  "anatomic_neoplasm_subdivision", "tumor_status"),

  exclude_neoadjuvant  = TRUE,             # 신보조요법 환자 제외
  cox_covariates       = c("age", "gender", "stage"),   # age는 연속형 (10세 단위)
  stage_collapsed      = "pathologic_stage_12_34",      # EPV 부족 시 stage 대신 사용 (없으면 stage로 생성)
  epv_min              = 10,               # 다변량 Cox EPV 경고 기준
  min_events           = 10,               # 유전자별 생존분석 최소 사건 수 (미만이면 건너뜀)
  tss_min_n            = 10,               # strata(tss) 민감도 분석: 이 수 미만 기관은 "Other"
  km_times             = c(36, 60),        # 3년/5년 생존율 (개월)
  ph_split_months      = 24,               # PH 위반 유전자: 0–24개월 / >24개월 시간 분할 Cox
  rmst_tau             = 60,               # RMST 차이 계산 시점 (개월)
  median_tie_tolerance = 1,                # median-split 재현 QC 허용 오차 (명)
  group_colors         = c(Low = "#2E6FB7", High = "#C8442F"),   # 모든 그림 공통

  gdc_dir       = "D:/GDCdata",            # 02a 다운로드 폴더 (Windows 경로 길이 문제 회피)
  processed_dir = "data/processed",
  output_dir    = "output",
  manifest      = "data/processed/db_manifest.json"   # .sav 지문 (md5, n, 생존 n/사건 등)
)

# ---- 원본 열 이름 스키마 (.sav의 실제 열 이름) ------------------------------
# 모든 코드는 이 스키마로 원본 열을 찾고 표준 이름(sample_id, surv_time, status,
# pathologic_stage, age, gender, neoadjuvant, last_contact)으로 바꿔서 사용.
# 암종별로 다른 항목만 cfg$schema$<CANCER>에 적음 (modifyList로 병합).
# 선택 항목(neoadjuvant, last_contact)이 없는 암종은 NULL 대신 NA로 지정.
cfg$schema_default <- list(
  id           = "sampleID",
  time         = "Days",
  time_unit    = "days",                   # "days" 또는 "months"
  status       = "Status",
  event_value  = 1,                        # 사건(사망)을 뜻하는 값 (예: 1, "Dead")
  stage        = "pathologic_Stage",
  stage_format = "numeric",                # "numeric" (1–4) 또는 "ajcc_text" ("Stage IIIB")
  age          = "age_at_initial_pathologic_diagnosis",
  sex          = "gender",
  neoadjuvant  = "history_neoadjuvant_treatment",
  last_contact = "last_contact_days_to"
)
cfg$schema <- list(
  COAD = list(),
  READ = list()
  # 예: LIHC = list(stage = "ajcc_pathologic_tumor_stage", stage_format = "ajcc_text")
)

# ---- 암종별 값 라벨 / 라벨 일관성 규칙 재정의 (기본값은 R/utils.R) ----------
cfg$value_labels <- list(COAD = list(), READ = list())
cfg$label_rules  <- list(COAD = list(), READ = list())
