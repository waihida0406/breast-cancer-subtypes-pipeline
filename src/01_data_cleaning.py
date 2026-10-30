"""
Step 1 - Data cleaning and exploratory data analysis (Objective 1)

Input : data/raw/METABRIC_RNA_Mutation.csv, built by Step 0 from cBioPortal
        or downloaded from Kaggle ("Breast Cancer Gene Expression Profiles
        (METABRIC)"); 1,904 patients, clinical attributes, mRNA z-scores,
        mutations
Output: data/interim/patients_clean.csv      clinical + derived variables
        data/interim/expression.csv          patient x gene mRNA z-scores
        data/interim/mutations.csv           patient x gene 0/1 mutation flags
        outputs/tables/cleaning_report.csv, table1_cohort_characteristics.csv
        outputs/figures/fig1-fig3

Why each cleaning step matters is written next to it, so the same
reasoning can go into the Materials & Methods chapter.
"""
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from config import (METABRIC_CSV, INTERIM, TAB, FIG, SUBTYPES, SUBTYPE_ORDER,
                    CLINICAL_RENAME)

COLORS = {"LumA": "#2b6cb0", "LumB": "#63b3ed", "HER2E": "#c05621",
          "Basal": "#9b2c2c", "ClaudinLow": "#b7791f", "Normal": "#718096"}


def normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = [re.sub(r"\s+", "_", c.strip().lower()) for c in df.columns]
    return df


def split_blocks(df: pd.DataFrame):
    """METABRIC mixes three blocks in one file: clinical, mRNA, mutation."""
    pam = next(c for c in df.columns if "pam50" in c)
    mut_cols = [c for c in df.columns if c.endswith("_mut")]
    known = set(CLINICAL_RENAME) | {pam, "patient_id", "overall_survival",
                                    "er_status_measured_by_ihc",
                                    "her2_status_measured_by_snp6",
                                    "primary_tumor_laterality", "oncotree_code",
                                    "cancer_type", "tumor_other_histologic_subtype",
                                    "3-gene_classifier_subtype"}
    expr_cols = [c for c in df.columns
                 if c not in known and c not in mut_cols
                 and pd.api.types.is_numeric_dtype(df[c])]
    return pam, expr_cols, mut_cols


def pos_neg(s: pd.Series) -> pd.Series:
    """Receptor status: handles 'Positive', 'Positve' (typo in the source), 'Neg'."""
    s = s.astype(str).str.strip().str.lower()
    return np.where(s.str.startswith("pos"), "Positive",
                    np.where(s.str.startswith("neg"), "Negative", None))


def clean(df: pd.DataFrame):
    log, n0 = [], len(df)
    df = normalise_columns(df)
    pam, expr_cols, mut_cols = split_blocks(df)

    # 1. Subtype label: keep the six PAM50 + claudin-low classes; drop 'NC'
    #    (not classified) because it has too few patients to analyse.
    df[pam] = df[pam].astype(str).str.strip()
    keep = df[pam].isin(SUBTYPES)
    log.append({"step": "drop unclassified / missing PAM50 subtype",
                "rows_removed": int((~keep).sum())})
    df = df[keep].copy()
    df["subtype_code"] = df[pam].map(lambda s: SUBTYPES[s][0])

    # 2. Duplicate patients would be counted twice in every test
    dup = df["patient_id"].duplicated()
    log.append({"step": "drop duplicate patient_id", "rows_removed": int(dup.sum())})
    df = df[~dup]

    # 3. Plausible age; missing values in other fields are kept and handled
    #    per analysis (complete-case for Cox, imputation inside the ML pipeline)
    df["age_at_diagnosis"] = pd.to_numeric(df["age_at_diagnosis"], errors="coerce")
    bad_age = ~df["age_at_diagnosis"].between(18, 100)
    log.append({"step": "drop missing / implausible age", "rows_removed": int(bad_age.sum())})
    df = df[~bad_age]

    # 4. Survival outcome must be present for survival analysis
    df["overall_survival_months"] = pd.to_numeric(df["overall_survival_months"], errors="coerce")
    no_surv = df["overall_survival_months"].isna() | (df["overall_survival_months"] < 0)
    log.append({"step": "drop missing survival time", "rows_removed": int(no_surv.sum())})
    df = df[~no_surv]

    # 5. Build the clinical table with tidy names and derived variables
    c = df[["patient_id"] + list(CLINICAL_RENAME)].rename(columns=CLINICAL_RENAME)
    c["patient_id"] = "MB-" + c["patient_id"].astype(str).str.replace("MB-", "", regex=False).str.zfill(4)
    c["subtype_code"] = df["subtype_code"].values
    for col in ["er_status", "pr_status", "her2_status"]:
        c[col] = pos_neg(c[col])
    for col in ["tumor_size_mm", "tumor_stage", "grade", "nodes_positive", "npi",
                "mutation_count", "chemotherapy", "hormone_therapy", "radiotherapy"]:
        c[col] = pd.to_numeric(c[col], errors="coerce")
    c["hr_positive"] = ((c.er_status == "Positive") | (c.pr_status == "Positive")).astype(int)
    c["her2_positive"] = (c.her2_status == "Positive").astype(int)
    c["triple_negative"] = ((c.er_status == "Negative") & (c.pr_status == "Negative")
                            & (c.her2_status == "Negative")).astype(int)
    c["vital_status"] = c["vital_status"].fillna("Unknown").astype(str).str.strip()
    c["os_event"] = c.vital_status.str.startswith("Died").astype(int)
    c["dss_event"] = (c.vital_status == "Died of Disease").astype(int)
    c["age_band"] = pd.cut(c.age, [0, 39, 49, 59, 69, 120],
                           labels=["<40", "40-49", "50-59", "60-69", "70+"]).astype(str)
    c["histology"] = c["histology"].fillna("Unknown").astype(str).str.replace(
        "Breast ", "", regex=False)
    log.append({"step": "receptor status harmonised (incl. 'Positve' typo)", "rows_removed": 0})

    # 6. Expression block: z-scores; drop genes with >20% missing
    expr = df[expr_cols].apply(pd.to_numeric, errors="coerce")
    gmiss = expr.isna().mean()
    expr = expr.loc[:, gmiss <= 0.2]
    log.append({"step": f"expression genes kept (of {len(expr_cols)})",
                "rows_removed": 0, "note": expr.shape[1]})
    expr = pd.concat([pd.Series(c["patient_id"].values, name="patient_id", index=expr.index), expr], axis=1)

    # 7. Mutation block: protein-change strings -> 0/1 present flags
    mut = df[mut_cols].astype(str).apply(lambda s: s.str.strip())
    mut = (~mut.isin(["0", "0.0", "nan", "", "None"])).astype(int)
    mut.columns = [m[:-4] for m in mut.columns]
    c["pik3ca_mut"] = mut["pik3ca"].values if "pik3ca" in mut else 0
    c["brca_mut"] = (mut.get("brca1", 0) + mut.get("brca2", 0) > 0).astype(int).values \
        if ("brca1" in mut or "brca2" in mut) else 0
    mut = pd.concat([pd.Series(c["patient_id"].values, name="patient_id", index=mut.index), mut], axis=1)

    log.append({"step": "final patients", "rows_removed": n0 - len(c), "note": len(c)})
    return c.reset_index(drop=True), expr.reset_index(drop=True), mut.reset_index(drop=True), log


def cohort_table(c: pd.DataFrame) -> pd.DataFrame:
    """Table 1 of the thesis: cohort characteristics by subtype."""
    g = c.groupby("subtype_code")
    t = pd.DataFrame({
        "n": g.size(),
        "percent": (100 * g.size() / len(c)).round(1),
        "age_mean": g.age.mean().round(1),
        "age_sd": g.age.std().round(1),
        "tumor_size_median_mm": g.tumor_size_mm.median(),
        "grade3_pct": (100 * g.grade.apply(lambda s: (s == 3).mean())).round(1),
        "node_positive_pct": (100 * g.nodes_positive.apply(lambda s: (s > 0).mean())).round(1),
        "er_pos_pct": (100 * g.er_status.apply(lambda s: (s == "Positive").mean())).round(1),
        "her2_pos_pct": (100 * g.her2_positive.mean()).round(1),
        "tnbc_pct": (100 * g.triple_negative.mean()).round(1),
        "deaths_from_disease_pct": (100 * g.dss_event.mean()).round(1),
    }).reindex(SUBTYPE_ORDER)
    return t


def eda_figures(c: pd.DataFrame):
    order = [s for s in SUBTYPE_ORDER if s in set(c.subtype_code)]

    # Fig 1 - class distribution (motivates imbalance handling in Step 7)
    counts = c.subtype_code.value_counts().reindex(order)
    fig, ax = plt.subplots(figsize=(7.5, 4))
    ax.bar(order, counts.values, color=[COLORS[s] for s in order])
    for i, v in enumerate(counts.values):
        ax.text(i, v + 8, f"{v}\n({v / len(c):.0%})", ha="center", fontsize=8)
    ax.set_ylim(0, counts.max() * 1.2)
    ax.set_ylabel("Number of patients"); ax.set_xlabel("PAM50 / claudin-low subtype")
    ax.set_title("Figure 1. Distribution of molecular subtypes")
    fig.tight_layout(); fig.savefig(FIG / "fig1_subtype_distribution.png", dpi=200); plt.close(fig)

    # Fig 2 - age and tumour size by subtype
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, col, lab in [(axes[0], "age", "Age at diagnosis (years)"),
                         (axes[1], "tumor_size_mm", "Tumour size (mm)")]:
        data = [c.loc[c.subtype_code == s, col].dropna() for s in order]
        bp = ax.boxplot(data, patch_artist=True, showfliers=False)
        for patch, s in zip(bp["boxes"], order):
            patch.set_facecolor(COLORS[s]); patch.set_alpha(0.7)
        ax.set_xticks(range(1, len(order) + 1), order, rotation=20)
        ax.set_ylabel(lab)
    fig.suptitle("Figure 2. Age and tumour size by subtype")
    fig.tight_layout(); fig.savefig(FIG / "fig2_age_tumour_size.png", dpi=200); plt.close(fig)

    # Fig 3 - receptor status by subtype (IHC/SNP6 vs molecular subtype)
    rec = c.groupby("subtype_code").agg(
        ER=("er_status", lambda s: (s == "Positive").mean()),
        PR=("pr_status", lambda s: (s == "Positive").mean()),
        HER2=("her2_positive", "mean"),
        TNBC=("triple_negative", "mean")).reindex(order)
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(order)); w = 0.2
    for i, (col, colr) in enumerate(zip(rec.columns, ["#2b6cb0", "#63b3ed", "#c05621", "#9b2c2c"])):
        ax.bar(x + (i - 1.5) * w, rec[col] * 100, w, label=f"{col}{'' if col == 'TNBC' else '+'}", color=colr)
    ax.set_xticks(x, order); ax.set_ylim(0, 105)
    ax.set_ylabel("% of patients"); ax.legend(ncol=4, fontsize=8, loc="upper right")
    ax.set_title("Figure 3. Receptor status by molecular subtype")
    fig.tight_layout(); fig.savefig(FIG / "fig3_receptor_status.png", dpi=200); plt.close(fig)


def main():
    raw = pd.read_csv(METABRIC_CSV, low_memory=False)
    print(f"[step 1] raw file: {raw.shape[0]} rows x {raw.shape[1]} columns")
    c, expr, mut, log = clean(raw)
    c.to_csv(INTERIM / "patients_clean.csv", index=False)
    expr.to_csv(INTERIM / "expression.csv", index=False)
    mut.to_csv(INTERIM / "mutations.csv", index=False)
    pd.DataFrame(log).to_csv(TAB / "cleaning_report.csv", index=False)
    t1 = cohort_table(c)
    t1.to_csv(TAB / "table1_cohort_characteristics.csv")
    eda_figures(c)
    print(f"[step 1] {len(c)} patients, {expr.shape[1] - 1} expression genes, "
          f"{mut.shape[1] - 1} mutation genes")
    print(t1.to_string())


if __name__ == "__main__":
    main()
