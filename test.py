import pandas as pd
import os
import pyreadstat

super_dir = os.path.dirname(os.path.abspath(__file__))
db_dir = os.path.join(super_dir, "database")
genes_files_dir = os.path.join(db_dir, "gene_files")
temp_genes_files_dir = os.path.join(db_dir, "temp_gene_files")

df1, meta1  = pyreadstat.read_sav(os.path.join(db_dir, "LUAD.sav"))
df2, meta2  = pyreadstat.read_sav(os.path.join(db_dir, "LUSC.sav"))
adam10 = pd.read_csv(os.path.join(temp_genes_files_dir, "LUAD_50_50_ADAM10.csv"))
adam10.rename(columns = {"Patient" : "sampleID"}, inplace = True)
ect2 = pd.read_csv(os.path.join(temp_genes_files_dir, "LUAD_50_50_ECT2.csv"))
ect2.rename(columns = {"Patient" : "sampleID"}, inplace = True)


# a = df1
# b = adam10[["sampleID", "Expression", "Group"]]
# b.rename(columns = {"Expression" : "ADAM10expression", "Group" : "ADAM10group"}, inplace = True)
# c = ect2[["sampleID", "Expression", "Group"]]
# c.rename(columns = {"Expression" : "ECT2expression", "Group" : "ECT2group"}, inplace = True)
# a.drop(columns=["ADAM10group", "ADAM10expression"], inplace=True)
# a.drop(columns=["ECT2group", "ECT2expression", "ECT2group_x", "ECT2expression_y"], inplace=True)
# merged_df = pd.merge(a, b, on="sampleID", how="left")
# merged_df = pd.merge(merged_df, c, on="sampleID", how="left")
# merged_df.to_csv(os.path.join(temp_genes_files_dir, "LUAD_50_50_ADAM10_merged.csv"), index=False)
# pyreadstat.write_sav(merged_df, os.path.join(temp_genes_files_dir, "LUAD.sav"))

# a = df2
# ect2_lusc = pd.read_csv(os.path.join(temp_genes_files_dir, "LUSC_50_50_ECT2.csv"))
# b = ect2_lusc[["Patient", "Expression", "Group"]]
# b.rename(columns = {"Patient": "sampleID", "Expression" : "ECT2expression", "Group" : "ECT2group"}, inplace = True)
# a.drop(columns=["ECT2group", "ECT2expression"], inplace=True)
# merged_df = pd.merge(a, b, on="sampleID", how="left")
# merged_df.to_csv(os.path.join(temp_genes_files_dir, "LUSC_50_50_ECT2_merged.csv"), index=False)
# pyreadstat.write_sav(merged_df, os.path.join(temp_genes_files_dir, "LUSC.sav"))

a = df1
a.loc[a["sampleID"] == "TCGA-49-6742", "OS.time"] = 1678
a.loc[a["sampleID"] == "TCGA-49-6742", "OS"] = 0
print(a.loc[a["sampleID"] == "TCGA-49-6742", "OS.time"])
pyreadstat.write_sav(a, os.path.join(temp_genes_files_dir, "LUAD_os_changed.sav"))

