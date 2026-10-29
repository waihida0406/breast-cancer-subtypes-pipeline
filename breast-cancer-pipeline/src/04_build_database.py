"""
Step 4 - Integrate all sources into a relational SQLite database

Tables (sql/schema.sql): subtypes, patients, genes, patient_mutations,
patient_expression, drugs, patient_treatments, patient_drug_eligibility,
adverse_events, fda_yearly, trials, trial_interventions.
Views (sql/views.sql) feed Steps 5 and 8 and the Power BI export.

Eligibility rules (target_group in config.DRUGS), applied to each patient's
receptor and mutation status:
  HR+                   ER+ or PR+
  HR+/HER2-             ER+ or PR+, and HER2-
  PIK3CA-mut HR+/HER2-  as above, plus a PIK3CA mutation
  HER2+                 HER2+
  TNBC                  ER-, PR-, HER2-
  BRCA-mut HER2-        BRCA1 or BRCA2 mutation, and HER2-
  Any                   every patient
"""
import sqlite3
import pandas as pd
from config import INTERIM, DB_PATH, SQL_DIR, SUBTYPES, DRUGS

RULES = {
    "HR+":                  lambda p: p.hr_positive == 1,
    "HR+/HER2-":            lambda p: (p.hr_positive == 1) & (p.her2_positive == 0),
    "PIK3CA-mut HR+/HER2-": lambda p: (p.hr_positive == 1) & (p.her2_positive == 0) & (p.pik3ca_mut == 1),
    "HER2+":                lambda p: p.her2_positive == 1,
    "TNBC":                 lambda p: p.triple_negative == 1,
    "BRCA-mut HER2-":       lambda p: (p.brca_mut == 1) & (p.her2_positive == 0),
    "Any":                  lambda p: pd.Series(True, index=p.index),
}


def build():
    pts = pd.read_csv(INTERIM / "patients_clean.csv")
    expr = pd.read_csv(INTERIM / "expression.csv")
    mut = pd.read_csv(INTERIM / "mutations.csv")
    fda = pd.read_csv(INTERIM / "fda_drug_summary.csv")
    rx = pd.read_csv(INTERIM / "fda_reactions.csv")
    yr = pd.read_csv(INTERIM / "fda_reports_by_year.csv")
    trials = pd.read_csv(INTERIM / "trials.csv")
    inters = pd.read_csv(INTERIM / "trial_interventions.csv")

    con = sqlite3.connect(DB_PATH)
    con.executescript((SQL_DIR / "schema.sql").read_text())

    # subtypes
    pd.DataFrame([{"subtype_code": c, "subtype_name": n, "description": d}
                  for c, n, d in SUBTYPES.values()]) \
      .to_sql("subtypes", con, if_exists="append", index=False)

    # patients
    cols = ["patient_id", "subtype_code", "age", "age_band", "menopausal_state",
            "histology", "cellularity", "tumor_size_mm", "tumor_stage", "grade",
            "nodes_positive", "npi", "er_status", "pr_status", "her2_status",
            "hr_positive", "her2_positive", "triple_negative", "pik3ca_mut", "brca_mut",
            "integrative_cluster", "mutation_count", "surgery", "os_months",
            "os_event", "dss_event", "vital_status", "cohort"]
    pts[cols].to_sql("patients", con, if_exists="append", index=False)

    # genes (union of expression and mutation panels)
    e_genes = [g.upper() for g in expr.columns[1:]]
    m_genes = [g.upper() for g in mut.columns[1:]]
    genes = pd.DataFrame({"gene_symbol": sorted(set(e_genes) | set(m_genes))})
    genes["gene_id"] = range(1, len(genes) + 1)
    genes["in_expression_panel"] = genes.gene_symbol.isin(e_genes).astype(int)
    genes["in_mutation_panel"] = genes.gene_symbol.isin(m_genes).astype(int)
    genes.to_sql("genes", con, if_exists="append", index=False)
    gid = dict(zip(genes.gene_symbol, genes.gene_id))

    # mutations (present only) and expression, both in long format
    ml = mut.melt(id_vars="patient_id", var_name="gene", value_name="present")
    ml = ml[ml.present == 1]
    ml.assign(gene_id=ml.gene.str.upper().map(gid))[["patient_id", "gene_id"]] \
      .to_sql("patient_mutations", con, if_exists="append", index=False)
    el = expr.melt(id_vars="patient_id", var_name="gene", value_name="z_score").dropna()
    el.assign(gene_id=el.gene.str.upper().map(gid))[["patient_id", "gene_id", "z_score"]] \
      .to_sql("patient_expression", con, if_exists="append", index=False, chunksize=50000)

    # drugs (+ FAERS totals)
    drugs = pd.DataFrame([d[:4] for d in DRUGS],
                         columns=["generic_name", "drug_class", "target_group", "metabric_therapy"])
    drugs.insert(0, "drug_id", range(1, len(drugs) + 1))
    drugs = drugs.merge(fda, on="generic_name", how="left")
    drugs.to_sql("drugs", con, if_exists="append", index=False)
    did = dict(zip(drugs.generic_name, drugs.drug_id))

    # treatments recorded in METABRIC
    tr = []
    for flag, name in [("chemotherapy", "Chemotherapy"), ("hormone_therapy", "Hormone therapy"),
                       ("radiotherapy", "Radiotherapy")]:
        ids = pts.loc[pts[flag] == 1, "patient_id"]
        tr.append(pd.DataFrame({"patient_id": ids, "therapy": name}))
    pd.concat(tr).to_sql("patient_treatments", con, if_exists="append", index=False)

    # eligibility for current drugs
    el_rows = []
    for d in drugs.itertuples():
        ids = pts.loc[RULES[d.target_group](pts), "patient_id"]
        el_rows.append(pd.DataFrame({"patient_id": ids, "drug_id": d.drug_id}))
    pd.concat(el_rows).to_sql("patient_drug_eligibility", con, if_exists="append", index=False)

    # adverse events / yearly report counts
    rx = rx.assign(drug_id=rx.generic_name.map(did)).dropna(subset=["drug_id", "reaction"])
    rx.groupby(["drug_id", "reaction"], as_index=False)["report_count"].sum() \
      .to_sql("adverse_events", con, if_exists="append", index=False)
    yr = yr.assign(drug_id=yr.generic_name.map(did)).dropna(subset=["drug_id"])
    yr[["drug_id", "year", "report_count"]].to_sql("fda_yearly", con, if_exists="append", index=False)

    # trials
    trials.to_sql("trials", con, if_exists="append", index=False)
    inters = inters.assign(drug_id=inters.generic_name.map(did))
    inters[inters.nct_id.isin(trials.nct_id)][
        ["nct_id", "intervention_type", "intervention_name", "drug_id"]
    ].to_sql("trial_interventions", con, if_exists="append", index=False)

    con.executescript((SQL_DIR / "views.sql").read_text())
    con.commit()
    for t in ["patients", "genes", "patient_mutations", "patient_expression", "drugs",
              "patient_treatments", "patient_drug_eligibility", "adverse_events",
              "trials", "trial_interventions"]:
        n = con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        print(f"[step 4] {t:26s} {n:>9,}")
    con.close()


if __name__ == "__main__":
    build()
