"""
Step 7 - Predicting molecular subtype (Objective 4)

Research question: how well can routine clinical pathology (age, size,
grade, nodes, ER/PR/HER2) recover the PAM50 molecular subtype, compared
with a gene-expression panel? This mirrors a real clinical gap: most
hospitals report IHC receptor status, not PAM50.

Feature sets
  A. Clinical + IHC  - routinely available variables
  B. Gene expression - mRNA z-scores of the METABRIC gene panel
     (PAM50 is itself defined by expression, so B is expected to score
      higher; the useful finding is the size of the gap and which subtypes
      clinical data cannot separate)

Class-imbalance handling
  * Stratified 80/20 hold-out + repeated stratified 5-fold CV on train only
  * class_weight="balanced" (LogReg) / "balanced_subsample" (RF)
  * SMOTE inside an imblearn pipeline when imbalanced-learn is installed,
    so synthetic samples never leak into validation folds
  * Imputation and scaling inside the pipeline (fitted per fold, no leakage)
  * Balanced accuracy, macro-F1 and per-class recall; majority baseline

Outputs: outputs/tables/table6_model_comparison.csv,
         table7_classification_report_<set>.csv, permutation_importance_<set>.csv
         outputs/figures/fig9_model_comparison.png, fig10_confusion_matrices.png,
         fig11_feature_importance.png
"""
import json
import joblib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, RepeatedStratifiedKFold, cross_validate
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.metrics import (classification_report, ConfusionMatrixDisplay,
                             balanced_accuracy_score, f1_score)
from sklearn.inspection import permutation_importance
from config import INTERIM, PROCESSED, TAB, FIG, MODELS, RANDOM_STATE, SUBTYPE_ORDER

try:
    from imblearn.pipeline import Pipeline as ImbPipeline
    from imblearn.over_sampling import SMOTE
    HAS_IMBLEARN = True
except ImportError:
    HAS_IMBLEARN = False

SCORING = {"balanced_accuracy": "balanced_accuracy", "f1_macro": "f1_macro",
           "accuracy": "accuracy"}


def clinical_features(p: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({
        "age": p.age,
        "tumor_size_mm": p.tumor_size_mm,
        "grade": p.grade,
        "nodes_positive": p.nodes_positive,
        "er_positive": (p.er_status == "Positive").astype(float).where(p.er_status.notna()),
        "pr_positive": (p.pr_status == "Positive").astype(float).where(p.pr_status.notna()),
        "her2_positive": p.her2_positive.astype(float),
        "postmenopausal": (p.menopausal_state.astype(str).str.lower() == "post").astype(float),
        "lobular_histology": p.histology.astype(str).str.contains("Lobular", case=False).astype(float),
        "high_cellularity": (p.cellularity.astype(str).str.lower() == "high").astype(float),
    })


def models(n_min: int) -> dict:
    pre = [("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())]
    m = {
        "baseline_majority": DummyClassifier(strategy="most_frequent"),
        "logreg_balanced": Pipeline(pre + [("clf", LogisticRegression(
            class_weight="balanced", max_iter=5000, C=0.5))]),
        "rf_balanced": Pipeline([("impute", SimpleImputer(strategy="median")),
                                 ("clf", RandomForestClassifier(
                                     n_estimators=400, class_weight="balanced_subsample",
                                     min_samples_leaf=2, random_state=RANDOM_STATE, n_jobs=-1))]),
    }
    if HAS_IMBLEARN:
        k = max(1, min(5, int(n_min * 0.8) - 1))
        m["logreg_smote"] = ImbPipeline(pre + [
            ("smote", SMOTE(k_neighbors=k, random_state=RANDOM_STATE)),
            ("clf", LogisticRegression(max_iter=5000, C=0.5))])
    return m


def run_feature_set(name, X, y, ids):
    X_tr, X_te, y_tr, y_te, id_tr, id_te = train_test_split(
        X, y, ids, test_size=0.2, stratify=y, random_state=RANDOM_STATE)
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=RANDOM_STATE)
    cand = models(y_tr.value_counts().min())
    rows = []
    for mname, model in cand.items():
        s = cross_validate(model, X_tr, y_tr, cv=cv, scoring=SCORING, n_jobs=-1)
        rows.append({"feature_set": name, "model": mname,
                     **{f"{k}_mean": s[f"test_{k}"].mean() for k in SCORING},
                     **{f"{k}_sd": s[f"test_{k}"].std() for k in SCORING}})
        print(f"[step 7] {name:16s} {mname:18s} bal_acc={s['test_balanced_accuracy'].mean():.3f} "
              f"macroF1={s['test_f1_macro'].mean():.3f}")
    cvres = pd.DataFrame(rows)
    best_name = cvres[cvres.model != "baseline_majority"].sort_values(
        "f1_macro_mean", ascending=False).iloc[0].model
    best = cand[best_name].fit(X_tr, y_tr)
    pred = best.predict(X_te)
    slug = name.lower().replace(" + ", "_").replace(" ", "_")
    pd.DataFrame(classification_report(y_te, pred, output_dict=True, zero_division=0)).T \
      .round(3).to_csv(TAB / f"table7_classification_report_{slug}.csv")
    metrics = {"feature_set": name, "best_model": best_name,
               "test_balanced_accuracy": round(balanced_accuracy_score(y_te, pred), 4),
               "test_f1_macro": round(f1_score(y_te, pred, average="macro"), 4)}
    imp = permutation_importance(best, X_te, y_te, scoring="f1_macro",
                                 n_repeats=10 if X.shape[1] < 50 else 3,
                                 random_state=RANDOM_STATE, n_jobs=-1)
    impdf = pd.DataFrame({"feature": X.columns, "importance_mean": imp.importances_mean,
                          "importance_sd": imp.importances_std}) \
        .sort_values("importance_mean", ascending=False).round(4)
    impdf.to_csv(TAB / f"permutation_importance_{slug}.csv", index=False)
    joblib.dump(best, MODELS / f"best_model_{slug}.joblib")
    preds = pd.DataFrame({"patient_id": id_te, "feature_set": name, "actual": y_te,
                          "predicted": pred, "correct": (y_te == pred).astype(int)})
    return cvres, metrics, (y_te, pred, best_name), impdf, preds


def main():
    pts = pd.read_csv(INTERIM / "patients_clean.csv")
    expr = pd.read_csv(INTERIM / "expression.csv")
    df = pts.merge(expr, on="patient_id")
    y = df.subtype_code
    print("[step 7] class distribution:", y.value_counts().to_dict(),
          "| SMOTE available:", HAS_IMBLEARN)

    sets = {"Clinical + IHC": clinical_features(df),
            "Gene expression": df[expr.columns[1:]]}
    all_cv, all_metrics, cms, imps, preds = [], [], {}, {}, []
    for name, X in sets.items():
        cvres, met, cm, imp, pr = run_feature_set(name, X, y, df.patient_id)
        all_cv.append(cvres); all_metrics.append(met); cms[name] = cm; imps[name] = imp
        preds.append(pr)
    cvall = pd.concat(all_cv).round(4)
    cvall.to_csv(TAB / "table6_model_comparison.csv", index=False)
    (TAB / "test_metrics.json").write_text(json.dumps(all_metrics, indent=2))
    pd.concat(preds).to_csv(PROCESSED / "predictions.csv", index=False)
    print("[step 7] hold-out:", all_metrics)

    # Fig 9 - CV comparison, both feature sets
    fig, ax = plt.subplots(figsize=(9, 4.2))
    labels = [f"{r.feature_set}\n{r.model}" for r in cvall.itertuples()]
    x = np.arange(len(cvall)); w = 0.38
    ax.bar(x - w / 2, cvall.balanced_accuracy_mean, w, yerr=cvall.balanced_accuracy_sd,
           color="#2b6cb0", capsize=3, label="Balanced accuracy")
    ax.bar(x + w / 2, cvall.f1_macro_mean, w, yerr=cvall.f1_macro_sd,
           color="#90cdf4", capsize=3, label="Macro F1")
    ax.set_xticks(x, labels, fontsize=7); ax.set_ylim(0, 1.05); ax.legend(fontsize=8)
    ax.set_ylabel("Score (mean +/- SD, 5x3 CV)")
    ax.set_title("Figure 9. Subtype prediction: clinical vs gene-expression features")
    fig.tight_layout(); fig.savefig(FIG / "fig9_model_comparison.png", dpi=200); plt.close(fig)

    # Fig 10 - confusion matrices side by side
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, (name, (yt, yp, mname)) in zip(axes, cms.items()):
        labs = [s for s in SUBTYPE_ORDER if s in set(yt)]
        ConfusionMatrixDisplay.from_predictions(yt, yp, labels=labs, ax=ax, cmap="Blues",
                                                normalize="true", values_format=".2f",
                                                colorbar=False)
        ax.set_title(f"{name} ({mname})", fontsize=10)
        ax.tick_params(axis="x", rotation=30)
    fig.suptitle("Figure 10. Hold-out confusion matrices (row-normalised recall)")
    fig.tight_layout(); fig.savefig(FIG / "fig10_confusion_matrices.png", dpi=200); plt.close(fig)

    # Fig 11 - top features
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for ax, (name, imp) in zip(axes, imps.items()):
        top = imp.head(10).iloc[::-1]
        lab = top.feature.str.upper() if name == "Gene expression" else top.feature
        ax.barh(lab, top.importance_mean, xerr=top.importance_sd, color="#2b6cb0", capsize=3)
        ax.set_title(name, fontsize=10); ax.set_xlabel("Drop in macro-F1 when shuffled")
    fig.suptitle("Figure 11. Top 10 predictive features (permutation importance)")
    fig.tight_layout(); fig.savefig(FIG / "fig11_feature_importance.png", dpi=200); plt.close(fig)


if __name__ == "__main__":
    main()
