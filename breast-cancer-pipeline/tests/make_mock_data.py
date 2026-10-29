"""Generate SYNTHETIC test fixtures so the pipeline can be smoke-tested offline.

Writes a fake METABRIC-shaped CSV and fake API responses (in the JSON shape
of openFDA and ClinicalTrials.gov v2) into the raw/cache folders.
Numbers are random with rough subtype patterns: never report them.

Usage (on a throwaway copy of the project):
    python tests/make_mock_data.py && python run_pipeline.py
"""
import json
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, "src")
from config import RAW, METABRIC_CSV, DRUGS

rng = np.random.default_rng(0)
# approximate subtype sizes in the real METABRIC Kaggle file (n=1,904)
SIZES = {"LumA": 679, "LumB": 461, "Her2": 220, "Basal": 199,
         "claudin-low": 199, "Normal": 140, "NC": 6}
# P(ER+), P(PR+), P(HER2+), relative hazard, P(TP53 mut), P(PIK3CA mut)
PROFILE = {"LumA": (.98, .75, .05, 1.0, .12, .45), "LumB": (.95, .55, .12, 1.8, .30, .35),
           "Her2": (.45, .25, .65, 2.4, .55, .30), "Basal": (.12, .06, .08, 2.6, .80, .08),
           "claudin-low": (.30, .20, .08, 2.0, .45, .15), "Normal": (.85, .65, .08, 1.3, .20, .35),
           "NC": (.6, .4, .1, 1.5, .3, .3)}
EXPR_SHIFT = {  # gene -> subtype mean z-score
    "esr1": {"LumA": 1, "LumB": .9, "Her2": -.4, "Basal": -1.5, "claudin-low": -.8, "Normal": .5},
    "pgr": {"LumA": .9, "LumB": .3, "Her2": -.5, "Basal": -1, "claudin-low": -.5, "Normal": .5},
    "erbb2": {"Her2": 2.2, "LumB": .2},
    "grb7": {"Her2": 2.0},
    "foxa1": {"LumA": .8, "LumB": .8, "Her2": .2, "Basal": -1.6, "claudin-low": -.9},
    "gata3": {"LumA": .9, "LumB": .7, "Basal": -1.4, "claudin-low": -.9},
    "krt5": {"Basal": 1.8, "claudin-low": .8, "Normal": .9},
    "krt14": {"Basal": 1.6, "claudin-low": .9, "Normal": 1.0},
    "krt17": {"Basal": 1.9, "claudin-low": .7},
    "egfr": {"Basal": 1.2, "claudin-low": 1.0, "Normal": .6},
    "foxc1": {"Basal": 1.9, "claudin-low": .6},
    "mki67": {"LumB": .9, "Her2": 1.0, "Basal": 1.4, "LumA": -.6, "Normal": -.7},
    "ccnb1": {"LumB": 1.0, "Her2": 1.0, "Basal": 1.5, "LumA": -.7, "Normal": -.8},
    "aurka": {"LumB": .9, "Her2": .9, "Basal": 1.4, "LumA": -.6},
    "cldn3": {"claudin-low": -1.8, "LumA": .3},
    "cldn7": {"claudin-low": -1.9, "LumA": .3},
    "vim": {"claudin-low": 1.6, "Basal": .5},
}
OTHER_GENES = ["brca1", "brca2", "palb2", "pten", "tp53", "atm", "cdh1", "chek2", "nbn",
               "nf1", "stk11", "bard1", "mlh1", "msh2", "msh6", "pms2", "epcam", "rad51c",
               "rad51d", "rad50", "rb1", "akt1", "akt2", "pik3ca", "map3k1", "mtor", "notch1",
               "cdk4", "cdk6", "ccnd1", "myc", "kit", "bcl2", "bax", "casp8"]
MUT_GENES = ["pik3ca", "tp53", "muc16", "ahnak2", "kmt2c", "syne1", "gata3", "map3k1",
             "ahnak", "dnah11", "cdh1", "dnah2", "kmt2d", "ush2a", "ryr2", "dnah5", "herc2",
             "pde4dip", "akap9", "tg", "birc6", "utrn", "tbx3", "col6a3", "arid1a", "lama2",
             "dnah7", "ncor1", "cbfb", "pten", "brca1", "brca2", "akt1"]
BASE_MUT = {g: f for g, f in zip(MUT_GENES, np.linspace(.12, .01, len(MUT_GENES)))}


def metabric():
    rows, pid = [], 0
    for st, n in SIZES.items():
        er, pr, her2, hz, p53, pik = PROFILE[st]
        for _ in range(n):
            age = float(np.clip(rng.normal(61 if st.startswith("Lum") else 53, 12), 22, 96))
            size = float(np.clip(rng.lognormal(3.1 + .15 * (hz > 2), .45), 2, 180))
            grade = int(rng.choice([1, 2, 3], p=[.25, .55, .2] if st == "LumA" else [.03, .27, .7]))
            nodes = int(rng.poisson(1.2 * hz))
            npi = round(.2 * size / 10 + grade + (1 if nodes == 0 else 2 if nodes <= 3 else 3), 2)
            t_death = rng.exponential(300 / hz * (1 + .01 * (60 - age)))
            t_other = rng.exponential(250 * np.exp(-.04 * (age - 60)))
            t_cens = rng.uniform(10, 355)
            t = min(t_death, t_other, t_cens)
            vital = ("Died of Disease" if t == t_death else
                     "Died of Other Causes" if t == t_other else "Living")
            r = {"patient_id": pid, "age_at_diagnosis": round(age, 2),
                 "type_of_breast_surgery": rng.choice(["MASTECTOMY", "BREAST CONSERVING"]),
                 "cancer_type": "Breast Cancer",
                 "cancer_type_detailed": rng.choice(
                     ["Breast Invasive Ductal Carcinoma", "Breast Invasive Lobular Carcinoma",
                      "Breast Mixed Ductal and Lobular Carcinoma",
                      "Breast Invasive Mixed Mucinous Carcinoma"],
                     p=[.78, .12 if st == "LumA" else .05, .07, .03 if st == "LumA" else .10]),
                 "cellularity": rng.choice(["High", "Moderate", "Low"], p=[.5, .38, .12]),
                 "chemotherapy": int(rng.random() < (.1 if st.startswith("Lum") else .45)),
                 "pam50_+_claudin-low_subtype": st, "cohort": int(rng.integers(1, 6)),
                 "er_status_measured_by_ihc": "Positve" if rng.random() < er else "Negative",
                 "er_status": "Positive" if rng.random() < er else "Negative",
                 "neoplasm_histologic_grade": grade,
                 "her2_status_measured_by_snp6": "NEUTRAL",
                 "her2_status": "Positive" if rng.random() < her2 else "Negative",
                 "tumor_other_histologic_subtype": "Ductal/NST",
                 "hormone_therapy": int(rng.random() < er * .75),
                 "inferred_menopausal_state": "Post" if age > 50 else "Pre",
                 "integrative_cluster": str(rng.integers(1, 11)),
                 "primary_tumor_laterality": rng.choice(["Left", "Right"]),
                 "lymph_nodes_examined_positive": nodes, "mutation_count": int(rng.poisson(6)),
                 "nottingham_prognostic_index": npi, "oncotree_code": "IDC",
                 "overall_survival_months": round(t, 2),
                 "overall_survival": int(vital == "Living"),
                 "pr_status": "Positive" if rng.random() < pr else "Negative",
                 "radio_therapy": int(rng.random() < .6),
                 "3-gene_classifier_subtype": "ER+/HER2- Low Prolif",
                 "tumor_size": round(size, 0), "tumor_stage": int(rng.choice([1, 2, 3], p=[.35, .55, .1])),
                 "death_from_cancer": vital}
            for g, shifts in EXPR_SHIFT.items():
                r[g] = round(rng.normal(shifts.get(st, 0), .8), 4)
            for g in OTHER_GENES:
                r[g] = round(rng.normal(0, 1), 4)
            for g in MUT_GENES:
                p = {"tp53": p53, "pik3ca": pik}.get(g, BASE_MUT[g])
                p = p * (1.8 if g in ("gata3", "map3k1", "cdh1") and st == "LumA" else 1)
                r[f"{g}_mut"] = rng.choice(["E545K", "R175H", "H1047R", "X123_splice"]) \
                    if rng.random() < p else "0"
            rows.append(r); pid += 1
    df = pd.DataFrame(rows)
    # dirty values the cleaner must catch
    df.loc[5, "age_at_diagnosis"] = 999
    df.loc[9, "overall_survival_months"] = np.nan
    df.loc[15, "tumor_size"] = np.nan
    df.loc[20, "er_status"] = "Positve"
    df.loc[25, "pam50_+_claudin-low_subtype"] = " LumA "
    df = pd.concat([df, df.iloc[[30]]])          # one duplicated patient
    df.to_csv(METABRIC_CSV, index=False)


def openfda():
    d = RAW / "openfda"; d.mkdir(parents=True, exist_ok=True)
    rx = ["NEUTROPENIA", "NAUSEA", "FATIGUE", "DIARRHOEA", "ALOPECIA", "ARTHRALGIA",
          "HOT FLUSH", "EJECTION FRACTION DECREASED", "INTERSTITIAL LUNG DISEASE",
          "HYPERGLYCAEMIA", "DRUG INEFFECTIVE", "DEATH"]
    for g, *_ in DRUGS:
        s = g.replace(" ", "_"); tot = int(rng.integers(300, 30000))
        meta = lambda n: {"meta": {"results": {"total": n}}, "results": []}
        (d / f"{s}_total.json").write_text(json.dumps(meta(tot)))
        (d / f"{s}_serious.json").write_text(json.dumps(meta(int(tot * rng.uniform(.3, .8)))))
        (d / f"{s}_death.json").write_text(json.dumps(meta(int(tot * rng.uniform(.01, .12)))))
        (d / f"{s}_reactions.json").write_text(json.dumps({"results": [
            {"term": r, "count": int(rng.integers(10, max(11, tot // 5)))} for r in rx]}))
        (d / f"{s}_dates.json").write_text(json.dumps({"results": [
            {"time": f"{y}{m:02d}15", "count": int(rng.integers(1, 300))}
            for y in range(2010, 2026) for m in (3, 9)]}))


def trials():
    d = RAW / "clinicaltrials"; d.mkdir(parents=True, exist_ok=True)
    names = [g for g, *_ in DRUGS] + ["Herceptin", "T-DM1", "Trastuzumab deruxtecan",
                                      "Datopotamab deruxtecan", "Placebo", "Radiotherapy",
                                      "Exercise program", "Nab-paclitaxel", "Ibrance"]
    studies = []
    for i in range(1200):
        y = int(rng.integers(1999, 2026))
        studies.append({"hasResults": bool(rng.random() < .35), "protocolSection": {
            "identificationModule": {"nctId": f"NCT{i:08d}", "briefTitle": f"Mock breast cancer trial {i}"},
            "statusModule": {"overallStatus": str(rng.choice(["COMPLETED", "RECRUITING", "TERMINATED", "ACTIVE_NOT_RECRUITING", "UNKNOWN"])),
                             "startDateStruct": {"date": f"{y}-0{rng.integers(1, 9)}"},
                             "completionDateStruct": {"date": f"{y + 3}-01"}},
            "designModule": {"studyType": str(rng.choice(["INTERVENTIONAL", "OBSERVATIONAL"], p=[.75, .25])),
                             "phases": [str(rng.choice(["PHASE1", "PHASE2", "PHASE3", "PHASE4"]))],
                             "enrollmentInfo": {"count": int(rng.integers(20, 3000))}},
            "sponsorCollaboratorsModule": {"leadSponsor": {"name": "Mock Sponsor",
                                                           "class": str(rng.choice(["INDUSTRY", "OTHER", "NIH"]))}},
            "armsInterventionsModule": {"interventions": [
                {"type": "DRUG", "name": str(n)} for n in rng.choice(names, 2, replace=False)]},
        }})
    for p in range(2):   # two pages, to exercise pagination
        page = {"studies": studies[p * 600:(p + 1) * 600]}
        if p == 0:
            page["nextPageToken"] = "mock-token"
        (d / f"page_{p:03d}.json").write_text(json.dumps(page))


if __name__ == "__main__":
    metabric(); openfda(); trials()
    print("Synthetic fixtures written (NOT real data).")
