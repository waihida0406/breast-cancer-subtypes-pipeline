"""
Project : Integrative Analysis of Breast Cancer Molecular Subtypes, Drug Safety
          and Clinical Trial Trends using Genomic Data, Machine Learning and
          Public Health APIs
Author  : Waihida Yashmin, M.Sc. Biotechnology
          SRM Institute of Science and Technology, Kattankulathur
Module  : config.py - shared paths, constants and curated reference lists
"""
from pathlib import Path
import os

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RAW, INTERIM, PROCESSED = DATA / "raw", DATA / "interim", DATA / "processed"
OUT = ROOT / "outputs"
FIG, TAB, MODELS = OUT / "figures", OUT / "tables", OUT / "models"
SQL_DIR = ROOT / "sql"
DB_PATH = PROCESSED / "breast_cancer.db"
POWERBI_DIR = ROOT / "powerbi" / "data"

for p in (RAW, INTERIM, PROCESSED, FIG, TAB, MODELS, POWERBI_DIR):
    p.mkdir(parents=True, exist_ok=True)

# Kaggle "Breast Cancer Gene Expression Profiles (METABRIC)"
METABRIC_CSV = RAW / "METABRIC_RNA_Mutation.csv"
OPENFDA_KEY = os.getenv("OPENFDA_API_KEY")   # optional; raises the daily limit
RANDOM_STATE = 42

# ---------------------------------------------------------------- subtypes
# PAM50 + claudin-low labels as written in METABRIC -> short codes
SUBTYPES = {
    "LumA":        ("LumA",   "Luminal A",
                    "ER+, low proliferation; best prognosis"),
    "LumB":        ("LumB",   "Luminal B",
                    "ER+, high proliferation"),
    "Her2":        ("HER2E",  "HER2-enriched",
                    "ERBB2 amplicon highly expressed"),
    "Basal":       ("Basal",  "Basal-like",
                    "Mostly triple-negative; basal keratins high"),
    "claudin-low": ("ClaudinLow", "Claudin-low",
                    "Mostly triple-negative; low cell-adhesion genes"),
    "Normal":      ("Normal", "Normal-like",
                    "Resembles normal breast tissue"),
}
SUBTYPE_ORDER = ["LumA", "LumB", "HER2E", "Basal", "ClaudinLow", "Normal"]

# Clinical columns used from the METABRIC file (names after normalisation)
CLINICAL_RENAME = {
    "age_at_diagnosis": "age",
    "tumor_size": "tumor_size_mm",
    "tumor_stage": "tumor_stage",
    "neoplasm_histologic_grade": "grade",
    "lymph_nodes_examined_positive": "nodes_positive",
    "nottingham_prognostic_index": "npi",
    "er_status": "er_status",
    "pr_status": "pr_status",
    "her2_status": "her2_status",
    "inferred_menopausal_state": "menopausal_state",
    "cellularity": "cellularity",
    "cancer_type_detailed": "histology",
    "integrative_cluster": "integrative_cluster",
    "mutation_count": "mutation_count",
    "chemotherapy": "chemotherapy",
    "hormone_therapy": "hormone_therapy",
    "radio_therapy": "radiotherapy",
    "type_of_breast_surgery": "surgery",
    "overall_survival_months": "os_months",
    "death_from_cancer": "vital_status",
    "cohort": "cohort",
}

# Marker genes for the expression heatmap (used if present in the panel)
MARKER_GENES = ["esr1", "pgr", "foxa1", "gata3", "erbb2", "grb7", "krt5",
                "krt14", "krt17", "egfr", "foxc1", "mki67", "ccnb1", "aurka",
                "cldn3", "cldn7", "vim"]

# ------------------------------------------------------------------ drugs
# generic_name, drug_class, target_group, METABRIC therapy flag it maps to,
# openFDA search override (None = phrase search on the generic name).
# target_group drives patient eligibility (see 04_build_database.py).
# Curated from FDA labels; check current labels before citing.
DRUGS = [
    ("tamoxifen",          "Endocrine - SERM",        "HR+",        "Hormone therapy", None),
    ("anastrozole",        "Endocrine - AI",          "HR+",        "Hormone therapy", None),
    ("letrozole",          "Endocrine - AI",          "HR+",        "Hormone therapy", None),
    ("exemestane",         "Endocrine - AI",          "HR+",        "Hormone therapy", None),
    ("fulvestrant",        "Endocrine - SERD",        "HR+",        "Hormone therapy", None),
    ("palbociclib",        "CDK4/6 inhibitor",        "HR+/HER2-",  None, None),
    ("ribociclib",         "CDK4/6 inhibitor",        "HR+/HER2-",  None, None),
    ("abemaciclib",        "CDK4/6 inhibitor",        "HR+/HER2-",  None, None),
    ("alpelisib",          "PI3K inhibitor",          "PIK3CA-mut HR+/HER2-", None, None),
    ("trastuzumab",        "HER2 antibody",           "HER2+",      None,
     'patient.drug.openfda.generic_name.exact:"TRASTUZUMAB"'),
    ("pertuzumab",         "HER2 antibody",           "HER2+",      None, None),
    ("lapatinib",          "HER2 TKI",                "HER2+",      None, None),
    ("neratinib",          "HER2 TKI",                "HER2+",      None, None),
    ("tucatinib",          "HER2 TKI",                "HER2+",      None, None),
    ("ado-trastuzumab emtansine",  "HER2 ADC",        "HER2+",      None, None),
    ("fam-trastuzumab deruxtecan", "HER2 ADC",        "HER2+",      None, None),
    ("pembrolizumab",      "Immune checkpoint inhibitor", "TNBC",   None, None),
    ("sacituzumab govitecan", "TROP2 ADC",            "TNBC",       None, None),
    ("olaparib",           "PARP inhibitor",          "BRCA-mut HER2-", None, None),
    ("talazoparib",        "PARP inhibitor",          "BRCA-mut HER2-", None, None),
    ("cyclophosphamide",   "Chemotherapy",            "Any",        "Chemotherapy", None),
    ("doxorubicin",        "Chemotherapy",            "Any",        "Chemotherapy", None),
    ("paclitaxel",         "Chemotherapy",            "Any",        "Chemotherapy", None),
    ("docetaxel",          "Chemotherapy",            "Any",        "Chemotherapy", None),
    ("capecitabine",       "Chemotherapy",            "Any",        "Chemotherapy", None),
]

# Restrict openFDA reports to those with a breast-cancer indication, so
# drugs used in many cancers (pembrolizumab, chemotherapy) are comparable.
FDA_BREAST_ONLY = True
