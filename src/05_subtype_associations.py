"""
Step 5 - Clinical and genomic features by molecular subtype (Objective 2)

Categorical features (receptor status, grade, stage, histology, treatment,
mutations in recurrently mutated genes):
  H0: the feature is independent of PAM50 subtype
  Test: chi-square test of independence, with a 5,000-permutation Monte
  Carlo p-value (small cells make the asymptotic p unreliable),
  Benjamini-Hochberg FDR across all features, bias-corrected Cramer's V.
Continuous features (age, tumour size, positive nodes, NPI, mutation count):
  Kruskal-Wallis test (non-normal, unequal variances), epsilon-squared effect.

Outputs: outputs/tables/table2_chi_square.csv, table3_kruskal_wallis.csv,
         mutation_prevalence_by_subtype.csv, marker_expression_by_subtype.csv
         outputs/figures/fig4_mutation_heatmap.png, fig5_marker_expression.png,
         fig6_effect_sizes.png
"""
import sqlite3
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency, kruskal
from config import (DB_PATH, INTERIM, TAB, FIG, RANDOM_STATE, SUBTYPE_ORDER,
                    MARKER_GENES)

N_PERM = 5000
MIN_MUT_FREQ = 0.02     # test genes mutated in at least 2% of patients

CATEGORICAL = {
    "er_status": "ER status", "pr_status": "PR status", "her2_status": "HER2 status",
    "grade": "Histologic grade", "tumor_stage": "Tumour stage",
    "histology": "Histological type", "menopausal_state": "Menopausal state",
    "cellularity": "Cellularity", "chemotherapy": "Received chemotherapy",
    "hormone_therapy": "Received hormone therapy", "radiotherapy": "Received radiotherapy",
}
CONTINUOUS = {"age": "Age", "tumor_size_mm": "Tumour size (mm)",
              "nodes_positive": "Positive lymph nodes", "npi": "Nottingham Prognostic Index",
              "mutation_count": "Mutation count"}


# ---------------------------------------------------------------- statistics
def _chi2_stat(t: np.ndarray) -> float:
    e = t.sum(1, keepdims=True) * t.sum(0, keepdims=True) / t.sum()
    return float(((t - e) ** 2 / e).sum())


def perm_pvalue(x, y, observed, rng) -> float:
    """Monte Carlo p-value: shuffle subtype labels, recompute chi-square."""
    xi = pd.factorize(x)[0]; yi = pd.factorize(y)[0]
    r, k = xi.max() + 1, yi.max() + 1
    hits = 0
    for _ in range(N_PERM):
        t = np.bincount(xi * k + rng.permutation(yi), minlength=r * k).reshape(r, k)
        hits += _chi2_stat(t) >= observed - 1e-9
    return (hits + 1) / (N_PERM + 1)


def cramers_v_corrected(table: np.ndarray) -> float:
    chi2 = _chi2_stat(table)
    n = table.sum(); r, k = table.shape
    phi2 = max(0, chi2 / n - (k - 1) * (r - 1) / (n - 1))
    rc, kc = r - (r - 1) ** 2 / (n - 1), k - (k - 1) ** 2 / (n - 1)
    return float(np.sqrt(phi2 / max(1e-12, min(kc - 1, rc - 1))))


def bh_fdr(p: np.ndarray) -> np.ndarray:
    n = len(p); order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(q, 1)
    return out


def chi_square_table(df: pd.DataFrame, features: dict) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_STATE)
    rows = []
    for col, label in features.items():
        d = df[[col, "subtype_code"]].dropna()
        tab = pd.crosstab(d[col], d.subtype_code)
        if tab.shape[0] < 2:
            continue
        chi2, p, dof, exp = chi2_contingency(tab.values, correction=False)
        rows.append({
            "feature": col, "label": label,
            "type": "Mutation" if col.endswith("_mut") else "Clinical",
            "n": len(d), "chi2": round(chi2, 2), "dof": dof, "p_asymptotic": p,
            "p_permutation": perm_pvalue(d[col].values, d.subtype_code.values, chi2, rng),
            "pct_cells_expected_lt5": round(100 * (exp < 5).mean(), 1),
            "cramers_v": round(cramers_v_corrected(tab.values), 3)})
    res = pd.DataFrame(rows)
    res["q_fdr"] = bh_fdr(res.p_permutation.values)
    res["significant_q05"] = res.q_fdr < 0.05
    res["effect"] = pd.cut(res.cramers_v, [-1, 0.1, 0.3, 0.5, 2],
                           labels=["negligible", "small", "medium", "large"])
    return res.sort_values("cramers_v", ascending=False)


def kruskal_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for col, label in CONTINUOUS.items():
        d = df[[col, "subtype_code"]].dropna()
        groups = [g[col].values for _, g in d.groupby("subtype_code")]
        h, p = kruskal(*groups)
        n, k = len(d), len(groups)
        rows.append({"feature": col, "label": label, "n": n, "H": round(h, 2),
                     "dof": k - 1, "p_value": p,
                     "epsilon_squared": round(max(0, (h - k + 1) / (n - k)), 3),
                     **{f"median_{s}": d.loc[d.subtype_code == s, col].median()
                        for s in SUBTYPE_ORDER}})
    res = pd.DataFrame(rows)
    res["q_fdr"] = bh_fdr(res.p_value.values)
    return res


# ------------------------------------------------------------------ figures
def heatmap(mat: pd.DataFrame, title: str, path, cmap, vmin, vmax, fmt, cbar_label):
    fig, ax = plt.subplots(figsize=(7.5, 0.38 * len(mat) + 1.8))
    im = ax.imshow(mat.values, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(mat.shape[1]), mat.columns)
    ax.set_yticks(range(mat.shape[0]), mat.index)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            strong = abs(v - (vmin + vmax) / 2) > 0.3 * (vmax - vmin)
            ax.text(j, i, format(v, fmt), ha="center", va="center", fontsize=7,
                    color="white" if strong else "black")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.04, label=cbar_label)
    fig.tight_layout(); fig.savefig(path, dpi=200); plt.close(fig)


def effect_plot(res: pd.DataFrame):
    r = res.sort_values("cramers_v").tail(25)
    colors = ["#2b6cb0" if s else "#a0aec0" for s in r.significant_q05]
    fig, ax = plt.subplots(figsize=(7.5, 0.3 * len(r) + 1.5))
    ax.barh(r.label, r.cramers_v, color=colors)
    for x, lab in [(0.1, "small"), (0.3, "medium"), (0.5, "large")]:
        ax.axvline(x, ls="--", lw=0.8, color="grey")
        ax.text(x, len(r) - 0.4, lab, fontsize=7, color="grey", ha="center")
    ax.set_xlabel("Bias-corrected Cramer's V")
    ax.set_title("Figure 6. Association of clinical and mutation features with subtype\n"
                 "(blue = significant after FDR correction, q < 0.05)")
    fig.tight_layout(); fig.savefig(FIG / "fig6_effect_sizes.png", dpi=200); plt.close(fig)


def main():
    pts = pd.read_csv(INTERIM / "patients_clean.csv")
    mut = pd.read_csv(INTERIM / "mutations.csv")
    expr = pd.read_csv(INTERIM / "expression.csv")
    order = [s for s in SUBTYPE_ORDER if s in set(pts.subtype_code)]

    # Recurrently mutated genes join the categorical test set
    freq = mut.drop(columns="patient_id").mean().sort_values(ascending=False)
    recurrent = freq[freq >= MIN_MUT_FREQ].index.tolist()
    df = pts.drop(columns=["pik3ca_mut", "brca_mut"]).merge(mut[["patient_id"] + recurrent].rename(
        columns={g: f"{g}_mut" for g in recurrent}), on="patient_id")
    features = {**CATEGORICAL, **{f"{g}_mut": f"{g.upper()} mutation" for g in recurrent}}

    chi = chi_square_table(df, features)
    chi.to_csv(TAB / "table2_chi_square.csv", index=False)
    kw = kruskal_table(pts)
    kw.to_csv(TAB / "table3_kruskal_wallis.csv", index=False)
    effect_plot(chi)

    # Fig 4 - mutation landscape: top 15 genes x subtype
    prev = df.groupby("subtype_code")[[f"{g}_mut" for g in recurrent[:15]]].mean().T
    prev.index = [i[:-4].upper() for i in prev.index]
    prev = prev[order]
    prev.round(4).to_csv(TAB / "mutation_prevalence_by_subtype.csv")
    heatmap(prev * 100, "Figure 4. Mutation prevalence (%) of the 15 most mutated genes",
            FIG / "fig4_mutation_heatmap.png", "Blues", 0, max(60, prev.values.max() * 100),
            ".0f", "% of patients")

    # Fig 5 - marker gene expression validates the PAM50 labels
    markers = [g for g in MARKER_GENES if g in expr.columns]
    me = expr[["patient_id"] + markers].merge(pts[["patient_id", "subtype_code"]]) \
        .groupby("subtype_code")[markers].mean().T[order]
    me.index = [g.upper() for g in me.index]
    me.round(3).to_csv(TAB / "marker_expression_by_subtype.csv")
    heatmap(me, "Figure 5. Mean mRNA z-score of marker genes by subtype",
            FIG / "fig5_marker_expression.png", "RdBu_r", -1.5, 1.5, ".1f", "Mean z-score")

    # Cross-check: mutation prevalence from the SQL view equals pandas
    with sqlite3.connect(DB_PATH) as con:
        v = pd.read_sql("SELECT subtype_code, gene_symbol, prevalence "
                        "FROM v_mutation_prevalence_by_subtype", con)
    chk = v.pivot(index="gene_symbol", columns="subtype_code", values="prevalence") \
           .reindex(index=prev.index, columns=prev.columns)
    print(f"[step 5] SQL vs pandas max mutation-prevalence diff: {(chk - prev).abs().max().max():.4f}")
    print(f"[step 5] {int(chi.significant_q05.sum())}/{len(chi)} categorical features significant (q<0.05)")
    print(chi[["label", "chi2", "cramers_v", "q_fdr"]].head(10).to_string(index=False))
    print(kw[["label", "H", "epsilon_squared", "q_fdr"]].to_string(index=False))


if __name__ == "__main__":
    main()
