# 학생연구자료 내보내기 자동화 — 구현 명세

> 이 문서는 Claude Code가 구현할 작업 명세다. 시작 전에 CLAUDE.md, docs/outputs.md, docs/methods.md를 먼저 읽을 것.
> CLAUDE.md의 규칙(하드코딩 금지, `config.R`이 유일한 설정 위치, 한국어 메시지/영문 라벨, 불일치 시 `stop()`)을 그대로 따른다.

## 0. 목표

```
Rscript run_gene.R <GENE>          # 유전자 하나: 다운로드 → 분석 → 내보내기 (전체 래퍼)
Rscript 90_export.R [GENE] [--dry-run]   # 이미 돌린 결과만 학생연구자료로 내보내기
```
결과: `C:/Users/안지훈/Documents/0.학생연구자료/<GENE>(COAD&READ)/` 에 표준 구조로 정리된 CSV·TIFF + 양식에 채운 표 docx.

## 1. 확정된 결정 사항 (사용자 합의)

| 항목 | 결정 |
|---|---|
| 그림 | **TIFF만** 복사 (PDF 제외). PDF로만 나오는 `zph_*`는 복사 대상 아님 |
| 표 | **CSV만** 복사 (pipeline이 만드는 개별 .docx 제외) |
| GEO 폴더 이름 | `geo/<GSE>_<validates>/` 예: `geo/GSE39582_COAD/` — `<validates>`는 `cfg$geo_datasets$<GSE>$validates`에서 읽음 |
| 표 양식 채우기 | 이번 작업에 포함 (`table양식.docx` → `<GENE>_tables.docx`) |
| 전체 래퍼 | 이번 작업에 포함 (`run_gene.R`) |

## 2. 대상 폴더 구조

```
0.학생연구자료/<GENE>(COAD&READ)/
├── firstline/      <C>_table1_<G>.csv, <C>_cox_uni.csv, <C>_cox_multi_<G>.csv,
│                   <C>_survival_summary.csv, <C>_km_summary.csv, <C>_ph_tests.csv,
│                   <C>_correlation_<G>.csv        ← 신규 (§5)
├── figure/         <C>_km_<G>.tiff, <C>_forest_multi_<G>.tiff
├── gsea/data/      <C>_gsea_<collection>.csv, <C>_gsea_summary.csv, <C>_gsea_QC.csv
├── gsea/figure/    <C>_nes_<collection>.tiff, <C>_volcano.tiff, <C>_pca.tiff, <C>_enrichment_*.tiff
├── geo/<GSE>_<validates>/data/     <GSE>_cox_*.csv, _km_summary, _survival_summary, _ph_tests, _probe_QC, _pheno_QC
├── geo/<GSE>_<validates>/figure/   <GSE>_km_*.tiff, <GSE>_forest*.tiff
├── <G>_tables.docx                 ← 양식 채운 표 + KM 그림 (§6)
├── <G>_key_results.csv             ← 핵심 수치 한 장 요약 (§4.4)
├── _export_manifest.csv
└── _old/<YYYYMMDD_HHMM>/           ← 덮어쓰기 전 기존 파일 백업
```
- `<C>`는 `cfg$cancers` 순회. 폴더 이름 형식은 config(`folder_fmt`)에서 — `"%s(COAD&READ)"`은 현재 기존 ROCK2 폴더와 맞춘 값이며, 암종 목록에서 자동 생성해도 됨 (`paste(cfg$cancers, collapse="&")`).
- `de_results.csv`(수 MB)는 기본 제외, `include_de = TRUE`일 때만.
- 유전자와 무관한 `table1_overall`, `forest_genes`는 복사하지 않음 (필요하면 config로 켜기).

## 3. config.R 추가 항목

```r
cfg$export <- list(
  dest_root  = "C:/Users/안지훈/Documents/0.학생연구자료",
  folder_fmt = "%s(COAD&READ)",
  template   = "C:/Users/안지훈/Documents/0.학생연구자료/table양식.docx",
  fig_ext    = "tiff",
  table_ext  = "csv",
  include_de = FALSE,
  # 출력 패밀리 → 대상 하위 폴더. {G} = 유전자, {coll} = cfg$gsea_collections
  layout = list(
    tcga = list(
      firstline = c("table1_{G}", "cox_uni", "cox_multi_{G}", "survival_summary",
                    "km_summary", "ph_tests", "correlation_{G}"),
      figure    = c("km_{G}", "forest_multi_{G}")),
    gsea = list(data   = c("gsea_{coll}", "gsea_summary", "gsea_QC"),
                figure = c("nes_{coll}", "volcano", "pca", "enrichment_*")),
    geo  = list(data   = c("cox_*", "km_summary", "survival_summary", "ph_tests", "probe_QC", "pheno_QC"),
                figure = c("km_*", "forest*"))
  )
)
cfg$cor_vars <- list(genes = c("APC", "KRAS", "TP53"),            # §5 상관분석: primary gene + 이 유전자들
                     clinical = c(age = "Age", cea_level_pretreatment = "CEA"))
```
- 경로에 한글이 있으므로 `config.R`이 UTF-8로 저장되는지 확인 (Windows R 4.6은 UTF-8 기본). 문제가 되면 `Sys.getenv("MEDIN_EXPORT_ROOT")` 대체 경로도 허용.
- 새 설정은 docs/covariates.md가 아니라 **docs/outputs.md에 "export" 섹션**으로 문서화.

## 4. `90_export.R [GENE] [--dry-run]`

번호 90: 07(TIMER)·08(proteomics) 자리를 비워 둠.

### 4.1 사전 검증 (실패 시 `stop()`, 무엇을 다시 돌려야 하는지 메시지에 명시)
1. 유전자 = 인자 또는 `cfg$primary_gene`. `detect_genes()`로 merged 데이터에 존재하는지 확인.
2. **03·04 출력 존재**: 암종마다 `layout$tcga`의 모든 파일. 현재 `COAD_table1_MAD2L1.csv`, `READ_table1_MAD2L1.csv`가 **없음** → 03을 MAD2L1 추가 후 다시 안 돌린 상태. 이런 경우 "03_table1.R 재실행 필요"로 중단.
3. **신선도**: 각 출력의 mtime이 `data/processed/<C>_merged.rds`보다 오래됐으면 경고 (재실행 안 한 결과일 가능성).
4. **GSEA**: `output/tables/<C>/gsea/<G>/`가 없으면 경고 후 GSEA만 건너뜀 (`05_gsea.R <C>` 안내; primary가 아니면 `cfg$gsea_genes`). `_stale/<G>/`에만 있으면 그 사실을 알려 주되 자동으로 가져오지는 않음.
5. **GEO가 이 유전자의 결과인지 검증** — 파일명에 유전자가 없어서 다른 유전자 결과가 섞일 위험이 가장 큼.
   - 06이 `survival_summary_raw.csv`에 `gene` 열(과 `cox_gene_term`)을 기록하도록 수정하고, export는 그 값으로 검증. 06 수정은 의도된 출력 변경이므로 새 출력을 baseline으로 커밋.
   - 열이 없는 옛 파일이면 `probe_QC.csv`의 유전자 기호로 대조, 그것도 안 되면 GEO 건너뜀 + 경고.
   - 유전자가 다르면 GEO만 건너뜀 (오류 아님), 콘솔에 "06을 <G>로 재실행" 안내.
   - 06은 primary gene만 돌리므로 GENE ≠ primary면 GEO는 원칙적으로 건너뜀.
6. 제외 패턴: `~$*`(엑셀 잠금 파일; 현재 `output/tables/COAD/~$COAD_cox_multi_MAD2L1.csv` 존재), `_stale/`, `inspect/`.

### 4.2 복사 규칙
- 파일명은 원본 그대로 유지 (`<prefix>_...`).
- 대상에 같은 이름이 있으면 md5 비교: 같으면 건너뜀, 다르면 기존 파일을 `_old/<타임스탬프>/<상대경로>`로 **이동** 후 복사. 삭제는 하지 않음.
- 대상 폴더의 기존 파일 중 이번 매핑에 없는 것(원고 docx, 수동 파일 등)은 건드리지 않음.
- 복사 실패(파일 잠김 등)는 모아서 마지막에 보고하고 0이 아닌 종료 코드.

### 4.3 `--dry-run`
복사/이동 예정 목록(원본 → 대상, 상태: new / same / replace / skip+이유)만 출력. 파일 시스템 변경 없음.

### 4.4 `<G>_key_results.csv`
한 행 = 코호트(COAD, READ, `<GSE>`) × 종점 × 모형. 열: cohort, endpoint, model, gene_term, n, events, HR, lo, hi, p, logrank_p, rmst_diff, ph_flag, exploratory. 출처: 각 `survival_summary_raw.csv`(TCGA는 해당 유전자 행, GEO는 전체). GSEA는 별도 블록 또는 별도 열로 유의 경로 수(FDR < cfg 기준), 상위 NES 경로 3개(양/음).

### 4.5 `_export_manifest.csv`
source, dest, status, md5, size, source_mtime, exported_at, `git rev-parse HEAD`, `git status --porcelain`이 비어 있는지(dirty 여부), gene, cox_gene_term.

## 5. 신규 분석: 상관분석 (양식 Table 2·3용)

양식에는 "Pearson correlation analysis" 표(primary gene, APC, KRAS, TP53, Age, CEA의 R/P 행렬)가 있는데 **현재 파이프라인에는 이 분석이 없다.**
- 03_table1.R(또는 새 helper `cor_matrix()` in R/utils.R)에서 `<C>_correlation_<G>.csv` 생성: 변수 = primary gene + `cfg$cor_vars$genes`의 `<gene>_expression`(log2) + `cfg$cor_vars$clinical`의 연속형 원값. pairwise complete, Pearson; R과 p를 long 형식(var1, var2, r, p, n)으로 저장.
- 행렬에 넣을 유전자가 해당 암종에 없으면 경고 후 제외.
- 방법은 docs/methods.md에 한 줄 추가 (Pearson, pairwise complete obs, 다중비교 미보정 탐색적).
- 기존 출력에는 영향 없음 (새 파일만 추가) — `git diff --stat output/tables/`로 확인.

## 6. 양식 채우기: `<G>_tables.docx`

`table양식.docx` 구조 (현재 ROCK2로 채워진 예시):
| 순서 | 내용 | 원본 |
|---|---|---|
| Table 1 | COAD·READ 나란히, 각 High / Low / P-value (n만) — 55행 | `<C>_table1_<G>.csv` |
| Table 2·3 | COAD·READ Pearson R/P 행렬 | `<C>_correlation_<G>.csv` (§5) |
| Figure 1·2 | COAD·READ KM 곡선 (그림 2개 삽입) | `<C>_km_<G>.tiff` |
| Table 4·5 | COAD·READ 단변량 + 다변량 Cox (Variable / Univariate HR(95% CI) / P / Multivariate HR(95% CI) / P) | `<C>_cox_uni.csv` + `<C>_cox_multi_<G>.csv` |

구현:
- R `officer` + `flextable`로 작성 (R만 쓰는 파이프라인 유지). 양식에서 글꼴·크기·선 스타일을 읽어 재현하고, 캡션은 양식 문장에서 유전자 이름만 바꿔 씀 (`ROCK2` → `<G>`; 하드코딩하지 말고 템플릿 캡션을 읽어 치환하거나 config에 캡션 템플릿을 둠).
- 양식 docx 자체를 직접 수정하지 말고 **새 파일**로 저장. 양식은 읽기 전용.
- 행 구성은 양식의 고정 행이 아니라 **파이프라인 출력(`cfg$table1_vars`) 기준**. 양식에만 있고 파이프라인에 없는 변수(History of colon polyps, Family history, Colon polyps at procurement 등)는 넣지 않고, 콘솔에 "양식에 있으나 출력에 없는 변수" 목록을 출력.
- Table 1의 n/% 표기: 양식은 n만 → `table1` CSV의 "n (x%)"에서 n만 추출하는 옵션(`cfg$export$table1_counts_only = TRUE`).
- Cox 표의 유전자 행 라벨은 `cfg$cox_gene_term`을 따름: 현재 `"continuous"`이므로 다변량은 "MAD2L1 (per 1 log2 unit)"이지 "High vs Low"가 아님. 단변량은 cox_uni의 같은 척도 행을 사용 (group 모형이면 High vs Low 행). 두 척도를 섞지 말 것.
- KM 그림: TIFF를 그대로 넣되 officer가 TIFF를 못 넣으면 `magick`으로 임시 PNG 변환 후 삽입 (임시 파일은 학생연구자료에 남기지 않음).
- 숫자 형식: `fmt_hr()`/`fmt_p()` 결과 그대로 (이미 CSV에 서식 적용됨).

**양식에서 발견된 오류 — 고쳐서 생성하고 사용자에게 보고할 것:**
- Table 1의 CEA 행 수준이 "≤65 / >65"로 되어 있고 값이 Lymphovascular invasion 행과 같음 (복사 실수). 실제 코드는 ≤5 / >5 ng/mL.
- Table 4·5의 CEA 라벨 ">=5 vs <5" — 파이프라인 기준은 ">5 vs ≤5".
- Figure 2 캡션 "foroverall" 띄어쓰기 누락.
- `COAD_cox_uni.csv`에 연령(ageG) 단변량 행이 보이지 않음 (Sex, stage, CEA만 있음). `cfg$cox_covariates`에는 `ageG`가 있으므로 04에서 빠지는 이유를 확인하고, 의도된 것이 아니면 수정 (수정 시 04 결과 baseline 갱신).

## 7. 전체 래퍼: `run_gene.R <GENE> [--skip-download] [--from=<step>]`

순서 (각 단계는 별도 `Rscript` 프로세스, 종료 코드가 0이 아니면 즉시 중단):
```
1. python csv_download.py <GENE>            # OncoLnc → database/gene_files/<C>_50_50_<GENE>.csv
2. 02b_merge_genes.R
3. 03_table1.R                              # + 상관분석 (§5)
4. 04_survival.R
5. 05_gsea.R COAD  /  05_gsea.R READ        # cfg$cancers 순회, 암종마다 별도 프로세스
6. 06_external_geo.R
7. 90_export.R <GENE>
```
필요한 수정:
- **`csv_download.py`**: `TARGET_GENES`, `target1/2` 고정값 제거 → `sys.argv`(유전자 여러 개 허용)로 받고, 암종 목록은 인자 또는 기본값 COAD/READ. 다운로드 실패 시 0이 아닌 종료 코드. `database/`는 읽기 전용 입력 규칙이 있으나 `gene_files/`에 새 CSV를 추가하는 것은 이 스크립트의 원래 역할이므로 허용 — 단 같은 이름이 있으면 덮어쓰지 말고 md5 비교 후 다르면 중단.
- **primary gene 덮어쓰기**: `config.R`에서 `primary_gene = Sys.getenv("MEDIN_GENE", "MAD2L1")`. 래퍼는 자식 프로세스에 `MEDIN_GENE`을 설정. 
- Rscript 경로: `"C:/Program Files/R/R-4.6.1/bin/Rscript.exe"` (config 또는 래퍼 상단 한 곳에). Python 실행 파일도 한 곳에서 지정 (`Sys.which("python")` 우선).
- 단계별 시작/종료 시각, 소요 시간을 `output/run_logs/<GENE>_<타임스탬프>.log`에 기록 (`output/run_logs/`는 .gitignore에 추가).
- `--from=04`처럼 중간부터 재시작 가능. `--skip-download`는 gene CSV가 이미 있을 때.

주의 (래퍼 시작 시 콘솔에 경고로 출력):
- primary gene을 바꾸면 04의 BH 보정 대상이 바뀌어 **다른 유전자의 q값도 바뀐다**. CLAUDE.md의 "baseline 커밋" 규칙과 충돌하므로, 래퍼 완료 후 `git status`로 변경된 CSV 목록을 보여 주고 커밋은 사용자에게 맡긴다 (자동 커밋 금지).
- 05·06이 primary gene만 돌리므로 GSEA/GEO 결과는 마지막에 돌린 유전자의 것만 output에 남고, 이전 유전자 것은 `_stale/`로 이동됨 → 유전자마다 export까지 끝낸 뒤 다음 유전자로 넘어갈 것.

## 8. 기존 학생연구자료 정리 (구현 후 1회, 사용자 확인 후)
- `MAD2L1/` (지금은 하위 폴더 없이 평평): `_old/`로 옮기고 `run_gene.R MAD2L1 --skip-download`(또는 03부터) → `MAD2L1(COAD&READ)/`로 새로 생성.
- `ROCK2(COAD&READ)/`: 현재 output에 ROCK2의 GSEA·GEO 결과가 없음(`_stale` 또는 MAD2L1로 덮임). 자동 재생성하지 말고 그대로 둠. `geo/` 바로 아래(9/28)와 `geo/COAD/`(9/29)의 중복은 사용자에게 확인 후 `geo/GSE39582_COAD/`로 정리.
- `HIST1H2AC/`, `medin_1팀_발표자료/`: 다른 도구로 만든 결과 → 건드리지 않음.

## 9. 검증 (완료 조건)
1. `Rscript 90_export.R MAD2L1 --dry-run` 출력 확인 → 실제 실행.
2. `_export_manifest.csv`의 md5가 원본과 모두 일치.
3. 같은 명령 재실행 시 전부 `same`(변경 0건).
4. GEO 오염 테스트: scratch 폴더에서 06 출력의 gene 값을 다른 유전자로 바꾼 사본으로 export → GEO 건너뜀 확인 (CLAUDE.md의 scratch 테스트 방식 사용, 실제 출력 건드리지 않기).
5. `<G>_tables.docx`를 열어 Table 1 n이 `<C>_table1_<G>.csv`와, Cox HR이 cox CSV와 일치하는지 표본 대조. High/Low n이 CLAUDE.md의 cohort invariant(merge_QC)와 일치.
6. medin_R에서 `git diff --stat output/tables/`: §5 신규 파일과 §6에서 의도적으로 고친 것 외에 변경 없음.
7. CLAUDE.md의 Commands·Architecture에 90_export.R, run_gene.R, `cfg$export`, `cfg$cor_vars` 반영; docs/outputs.md에 export 섹션 추가.

## 10. 작업 순서 제안
1. §3 config 추가 + §4 `90_export.R` (dry-run부터)
2. 06에 gene 열 추가 → GEO 검증 동작 확인
3. §5 상관분석
4. §6 docx 생성
5. §7 `csv_download.py` 인자화 + `run_gene.R`
6. §9 검증 → §8 정리는 사용자 확인 후
