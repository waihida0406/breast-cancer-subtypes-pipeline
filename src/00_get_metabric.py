"""
Step 0 - Download METABRIC from cBioPortal and build the analysis file

The Kaggle file "Breast Cancer Gene Expression Profiles (METABRIC)" is a
derived subset of the cBioPortal study `brca_metabric`. This step rebuilds
an equivalent file directly from the original source, so the project needs
no Kaggle account and the provenance is clear:

  Source  : cBioPortal DataHub (github.com/cBioPortal/datahub,
            public/brca_metabric), Curtis et al. 2012; Pereira et al. 2016
  Cohort  : patients with BOTH mRNA expression and targeted sequencing
            (1,904 patients, the same set as the Kaggle file)
  Output  : data/raw/METABRIC_RNA_Mutation.csv, in the Kaggle column format
            (clinical columns, then mRNA z-scores, then <gene>_mut columns)

Expression panel: the 173 sequenced cancer genes plus the breast cancer
marker genes in config.MARKER_GENES (the Kaggle file used a similar
~330-gene cancer panel). Mutation columns hold the protein change of the
first non-silent mutation, or "0" if none.

If you already downloaded the Kaggle CSV, put it in data/raw/ and this step
is skipped.
"""
import io
import requests
import numpy as np
import pandas as pd
from config import RAW, METABRIC_CSV, MARKER_GENES

LFS = "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_metabric/"
RAWGH = "https://raw.githubusercontent.com/cBioPortal/datahub/master/public/brca_metabric/"
FILES = {  # file -> base URL (large files live in Git LFS)
    "data_clinical_patient.txt": LFS,
    "data_clinical_sample.txt": LFS,
    "data_mutations.txt": LFS,
    "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt": LFS,
    "case_lists/cases_sequenced.txt": RAWGH,
}
CACHE = RAW / "cbioportal"
NON_SILENT = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
              "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation",
              "Translation_Start_Site"}


def download():
    CACHE.mkdir(parents=True, exist_ok=True)
    for name, base in FILES.items():
        path = CACHE / name.split("/")[-1]
        if path.exists() and path.stat().st_size > 1000:
            continue
        print(f"[step 0] downloading {name} ...")
        with requests.get(base + name, stream=True, timeout=120) as r:
            r.raise_for_status()
            with open(path, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
        if path.read_bytes()[:40].startswith(b"version https://git-lfs"):
            raise RuntimeError(f"{name}: got a Git LFS pointer, not the data; check the URL")


def yes_no(s):
    return s.map({"YES": 1, "NO": 0})


def build() -> pd.DataFrame:
    read = lambda f: pd.read_csv(CACHE / f, sep="\t", comment="#", low_memory=False)
    pat, smp = read("data_clinical_patient.txt"), read("data_clinical_sample.txt")
    seq_ids = (CACHE / "cases_sequenced.txt").read_text().split("case_list_ids:")[1].split()

    # Expression: average duplicate probes per gene; transpose to patient x gene
    expr = pd.read_csv(CACHE / "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt",
                       sep="\t", low_memory=False)
    expr = expr.dropna(subset=["Hugo_Symbol"]).drop(columns="Entrez_Gene_Id")
    expr = expr.groupby("Hugo_Symbol").mean()

    # Cohort = patients with expression AND sequencing (as in the Kaggle file)
    cohort = sorted(set(expr.columns) & set(seq_ids))

    mut = read("data_mutations.txt")
    mut = mut[mut.Variant_Classification.isin(NON_SILENT)
              & mut.Tumor_Sample_Barcode.isin(cohort)]
    panel = sorted(read("data_mutations.txt").Hugo_Symbol.dropna().unique())

    genes = [g for g in dict.fromkeys(panel + [m.upper() for m in MARKER_GENES]) if g in expr.index]
    e = expr.loc[genes, cohort].T
    e.columns = [g.lower() for g in e.columns]

    # Clinical block in the Kaggle column names
    c = pat.merge(smp, on="PATIENT_ID").set_index("PATIENT_ID").loc[cohort]
    clin = pd.DataFrame({
        "patient_id": c.index.str.replace("MB-", "", regex=False),
        "age_at_diagnosis": c.AGE_AT_DIAGNOSIS,
        "type_of_breast_surgery": c.BREAST_SURGERY,
        "cancer_type": c.CANCER_TYPE,
        "cancer_type_detailed": c.CANCER_TYPE_DETAILED,
        "cellularity": c.CELLULARITY,
        "chemotherapy": yes_no(c.CHEMOTHERAPY),
        "pam50_+_claudin-low_subtype": c.CLAUDIN_SUBTYPE,
        "cohort": c.COHORT,
        "er_status_measured_by_ihc": c.ER_IHC,
        "er_status": c.ER_STATUS,
        "neoplasm_histologic_grade": c.GRADE,
        "her2_status_measured_by_snp6": c.HER2_SNP6,
        "her2_status": c.HER2_STATUS,
        "tumor_other_histologic_subtype": c.HISTOLOGICAL_SUBTYPE,
        "hormone_therapy": yes_no(c.HORMONE_THERAPY),
        "inferred_menopausal_state": c.INFERRED_MENOPAUSAL_STATE,
        "integrative_cluster": c.INTCLUST,
        "primary_tumor_laterality": c.LATERALITY,
        "lymph_nodes_examined_positive": c.LYMPH_NODES_EXAMINED_POSITIVE,
        "mutation_count": mut.groupby("Tumor_Sample_Barcode").size().reindex(cohort).fillna(0).astype(int).values,
        "nottingham_prognostic_index": c.NPI,
        "oncotree_code": c.ONCOTREE_CODE,
        "overall_survival_months": c.OS_MONTHS,
        "overall_survival": (c.OS_STATUS == "0:LIVING").astype(int),   # 1 = alive
        "pr_status": c.PR_STATUS,
        "radio_therapy": yes_no(c.RADIO_THERAPY),
        "3-gene_classifier_subtype": c.THREEGENE,
        "tumor_size": c.TUMOR_SIZE,
        "tumor_stage": c.TUMOR_STAGE,
        "death_from_cancer": c.VITAL_STATUS,
    }, index=c.index)

    # Mutation block: first non-silent protein change per patient x gene
    mut["change"] = mut.HGVSp_Short.fillna(mut.Variant_Classification).astype(str).str.replace("p.", "", regex=False)
    mm = mut.drop_duplicates(["Tumor_Sample_Barcode", "Hugo_Symbol"]) \
            .pivot(index="Tumor_Sample_Barcode", columns="Hugo_Symbol", values="change") \
            .reindex(index=cohort, columns=panel).fillna("0")
    mm.columns = [f"{g.lower()}_mut" for g in mm.columns]

    out = pd.concat([clin, e, mm], axis=1).reset_index(drop=True)
    return out


def main(force: bool = False):
    if METABRIC_CSV.exists() and not force:
        print(f"[step 0] {METABRIC_CSV.name} already present - skipping download")
        return
    download()
    df = build()
    df.to_csv(METABRIC_CSV, index=False)
    n_expr = sum(1 for c in df.columns if not c.endswith("_mut")) - 31
    print(f"[step 0] wrote {METABRIC_CSV.name}: {len(df)} patients, "
          f"{n_expr} expression genes, {sum(c.endswith('_mut') for c in df.columns)} mutation genes")


if __name__ == "__main__":
    import sys
    main(force="--force" in sys.argv)
