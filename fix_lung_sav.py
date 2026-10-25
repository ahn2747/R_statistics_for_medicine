import pandas as pd
import os
import pyreadstat

# LUAD.sav / LUSC.sav 열 이름 정리 (값은 바꾸지 않음)
# 입력: database/<C>.sav  →  출력: database/temp_gene_files/<C>_fixed.{csv,sav} (database/<C>.sav는 덮어쓰지 않음)
# 규칙: 유전자 열 = <대문자 포함 유전자 기호>Expression + <같은 기호>Group 짝

super_dir = os.path.dirname(os.path.abspath(__file__))
db_dir = os.path.join(super_dir, "database")
temp_genes_files_dir = os.path.join(db_dir, "temp_gene_files")


def same_values(df, a, b):
    # 두 열의 값이 완전히 같은지 (결측 위치 포함)
    x, y = df[a], df[b]
    return bool(((x == y) | (x.isna() & y.isna())).all())


def drop_duplicates(df, pairs):
    # (지울 열, 남길 열): 값이 같을 때만 지움
    for drop, keep in pairs:
        if not same_values(df, drop, keep):
            raise ValueError(f"{drop} 와 {keep} 값이 다름 → 지우지 않음")
        df.drop(columns=[drop], inplace=True)
        print(f"  삭제: {drop} (= {keep})")


def rename(df, mapping):
    miss = [k for k in mapping if k not in df.columns]
    if miss:
        raise ValueError(f"없는 열: {miss}")
    df.rename(columns=mapping, inplace=True)
    for k, v in mapping.items():
        print(f"  이름 변경: {k} → {v}")


def check(df, cancer):
    # 대소문자 무시 중복 이름, Expression/Group 짝, 소문자 기호 확인
    lower = pd.Series(df.columns).str.lower()
    dup = list(df.columns[lower.duplicated(keep=False).values])
    if dup:
        raise ValueError(f"{cancer}: 대소문자만 다른 중복 열 {dup}")
    stems = {}
    for col in df.columns:
        low = col.lower()
        for suf in ("expression", "group"):
            if low.endswith(suf) and len(col) > len(suf):
                stems.setdefault(col[: -len(suf)], set()).add(suf)
    unpaired = {s: sorted(v) for s, v in stems.items() if len(v) == 1}
    lower_stem = [s for s in stems if not any(ch.isupper() for ch in s)]
    print(f"  짝 없는 Expression/Group: {unpaired}")
    print(f"  소문자 기호: {lower_stem}")
    print(f"  유전자 (짝 있음) {sum(len(v) == 2 for v in stems.values())}개")


# ---- LUAD ----
print("== LUAD")
a, meta_a = pyreadstat.read_sav(os.path.join(db_dir, "LUAD.sav"))
drop_duplicates(a, [("TREM1expression_x", "TREM1expression"), ("TREM1group_x", "TREM1group"),
                    ("TREM1expression_y", "TREM1expression"), ("TREM1group_y", "TREM1group"),
                    ("BRCA1expression", "BRCA1Expression"), ("BRCA1group", "BRCA1Group")])
rename(a, {"SFTPBEexpression": "SFTPBExpression",
           "c10orf91group": "C10orf91Group",
           "gas6group": "GAS6Group"})
check(a, "LUAD")
a.to_csv(os.path.join(temp_genes_files_dir, "LUAD_fixed.csv"), index=False)
pyreadstat.write_sav(a, os.path.join(temp_genes_files_dir, "LUAD_fixed.sav"))

# ---- LUSC ----
print("== LUSC")
b, meta_b = pyreadstat.read_sav(os.path.join(db_dir, "LUSC.sav"))
drop_duplicates(b, [("brca1expression", "BRCA1expression"), ("brca1group", "BRCA1group")])
lower_genes = ["cdk1", "cdk2", "ccna2", "dlgap5", "e2f1", "ccnb1", "plk1", "myc", "smad3", "vim",
               "ocln", "esr1", "pparg", "scd", "kras", "hspa5", "tmem125", "linc00261"]
mapping = {"RRM1Expresion": "RRM1Expression", "c": "RRM2Expression",
           "c10orf91expression": "C10orf91Expression", "c10orf91group": "C10orf91Group"}
for g in lower_genes:
    mapping[g + "expression"] = g.upper() + "Expression"
    mapping[g + "group"] = g.upper() + "Group"
rename(b, mapping)
check(b, "LUSC")
b.to_csv(os.path.join(temp_genes_files_dir, "LUSC_fixed.csv"), index=False)
pyreadstat.write_sav(b, os.path.join(temp_genes_files_dir, "LUSC_fixed.sav"))
