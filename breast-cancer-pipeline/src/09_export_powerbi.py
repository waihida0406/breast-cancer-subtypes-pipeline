"""
Step 9 (Objective 6) - Export star-schema CSVs for Power BI from the SQLite database.

Power BI: Get Data > Text/CSV (or Folder) on powerbi/data/, then build the
relationships listed in powerbi/DASHBOARD_SPEC.md.
"""
import sqlite3
import pandas as pd
from config import DB_PATH, POWERBI_DIR, PROCESSED, TAB

EXPORTS = {
    # dimensions
    "dim_subtype": "SELECT * FROM subtypes",
    "dim_gene": "SELECT * FROM genes",
    "dim_drug": "SELECT drug_id, generic_name, drug_class, target_group, metabric_therapy FROM drugs",
    "dim_trial": "SELECT * FROM trials",
    # facts / bridges
    "fact_patient": "SELECT * FROM patients",
    "fact_patient_mutation": "SELECT * FROM patient_mutations",
    # expression is large: export only the marker genes used on the dashboard
    "fact_marker_expression": """SELECT e.* FROM patient_expression e JOIN genes g USING (gene_id)
                                 WHERE g.gene_symbol IN ('ESR1','PGR','ERBB2','GRB7','FOXA1','GATA3',
                                 'KRT5','KRT14','KRT17','EGFR','MKI67','CCNB1','AURKA','FOXC1')""",
    "fact_patient_treatment": "SELECT * FROM patient_treatments",
    "bridge_patient_drug_eligibility": "SELECT * FROM patient_drug_eligibility",
    "fact_drug_safety": "SELECT * FROM v_drug_safety",
    "fact_adverse_event": "SELECT drug_id, reaction, report_count FROM adverse_events",
    "fact_fda_yearly": "SELECT * FROM fda_yearly",
    "bridge_trial_drug": "SELECT * FROM trial_interventions",
}


def main():
    with sqlite3.connect(DB_PATH) as con:
        for name, q in EXPORTS.items():
            pd.read_sql(q, con).to_csv(POWERBI_DIR / f"{name}.csv", index=False)
    # Analysis outputs as extra tables
    for f in ["table2_chi_square.csv", "table4_survival_by_subtype.csv",
              "table5_cox_model.csv", "km_curves.csv", "table6_model_comparison.csv",
              "permutation_importance_clinical_ihc.csv",
              "permutation_importance_gene_expression.csv"]:
        if (TAB / f).exists():
            pd.read_csv(TAB / f).to_csv(POWERBI_DIR / f"analysis_{f}", index=False)
    if (PROCESSED / "predictions.csv").exists():
        pd.read_csv(PROCESSED / "predictions.csv").to_csv(
            POWERBI_DIR / "fact_prediction.csv", index=False)
    print(f"[step 9] wrote {len(list(POWERBI_DIR.glob('*.csv')))} CSVs to {POWERBI_DIR}")


if __name__ == "__main__":
    main()
