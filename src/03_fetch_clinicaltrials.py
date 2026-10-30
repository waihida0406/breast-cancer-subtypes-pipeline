"""
Step 3 (Objective 4) - Breast cancer trials from the ClinicalTrials.gov v2 API.

Endpoint: https://clinicaltrials.gov/api/v2/studies?query.cond=breast cancer
Outputs (data/interim/):
  trials.csv               one row per study
  trial_interventions.csv  one row per study x intervention, matched to
                           our drug list where the name contains a generic
"""
import json
import requests
import pandas as pd
from config import RAW, INTERIM, DRUGS

BASE = "https://clinicaltrials.gov/api/v2/studies"
CACHE = RAW / "clinicaltrials"
CACHE.mkdir(parents=True, exist_ok=True)

# Search patterns per drug (generic, brand and code names). Order matters:
# the trastuzumab ADCs are checked before plain trastuzumab, and
# "trastuzumab deruxtecan" is used rather than "deruxtecan" so that
# datopotamab deruxtecan (a TROP2 ADC) is not mis-assigned.
PATTERNS = [
    ("ado-trastuzumab emtansine", ["emtansine", "t-dm1", "kadcyla"]),
    ("fam-trastuzumab deruxtecan", ["trastuzumab deruxtecan", "t-dxd", "ds-8201", "enhertu"]),
    ("trastuzumab", ["trastuzumab", "herceptin"]),
    ("pertuzumab", ["pertuzumab", "perjeta", "phesgo"]),
    ("tamoxifen", ["tamoxifen", "nolvadex"]),
    ("anastrozole", ["anastrozole", "arimidex"]),
    ("letrozole", ["letrozole", "femara"]),
    ("exemestane", ["exemestane", "aromasin"]),
    ("fulvestrant", ["fulvestrant", "faslodex"]),
    ("palbociclib", ["palbociclib", "ibrance"]),
    ("ribociclib", ["ribociclib", "kisqali"]),
    ("abemaciclib", ["abemaciclib", "verzenio"]),
    ("alpelisib", ["alpelisib", "piqray"]),
    ("lapatinib", ["lapatinib", "tykerb"]),
    ("neratinib", ["neratinib", "nerlynx"]),
    ("tucatinib", ["tucatinib", "tukysa"]),
    ("pembrolizumab", ["pembrolizumab", "keytruda", "mk-3475"]),
    ("sacituzumab govitecan", ["sacituzumab", "trodelvy"]),
    ("olaparib", ["olaparib", "lynparza"]),
    ("talazoparib", ["talazoparib", "talzenna"]),
    ("cyclophosphamide", ["cyclophosphamide", "cytoxan"]),
    ("doxorubicin", ["doxorubicin", "adriamycin"]),
    ("paclitaxel", ["paclitaxel", "taxol", "abraxane"]),
    ("docetaxel", ["docetaxel", "taxotere"]),
    ("capecitabine", ["capecitabine", "xeloda"]),
]
assert {g for g, _ in PATTERNS} == {d[0] for d in DRUGS}, "PATTERNS out of sync with DRUGS"


def fetch_all(condition: str = "breast cancer", page_size: int = 1000) -> list[dict]:
    studies, token, page = [], None, 0
    while True:
        cache = CACHE / f"page_{page:03d}.json"
        if cache.exists():
            data = json.loads(cache.read_text())
        else:
            params = {"query.cond": condition, "pageSize": page_size, "format": "json"}
            if token:
                params["pageToken"] = token
            r = requests.get(BASE, params=params, timeout=60)
            r.raise_for_status()
            data = r.json()
            cache.write_text(json.dumps(data))
        studies += data.get("studies", [])
        token = data.get("nextPageToken")
        page += 1
        if not token:
            return studies


def _year(date_str):
    try:
        return int(str(date_str)[:4])
    except (TypeError, ValueError):
        return None


def match_drug(name: str) -> str | None:
    n = (name or "").lower()
    for generic, pats in PATTERNS:
        if any(p in n for p in pats):
            return generic
    return None


def parse(studies: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    trials, inters = [], []
    for s in studies:
        p = s.get("protocolSection", {})
        ident, status = p.get("identificationModule", {}), p.get("statusModule", {})
        design = p.get("designModule", {})
        sponsor = p.get("sponsorCollaboratorsModule", {}).get("leadSponsor", {})
        nct = ident.get("nctId")
        start = status.get("startDateStruct", {}).get("date")
        phases = design.get("phases") or ["NA"]
        trials.append({
            "nct_id": nct,
            "title": ident.get("briefTitle"),
            "status": status.get("overallStatus"),
            "study_type": design.get("studyType"),
            "phase": "/".join(phases),
            "start_date": start,
            "start_year": _year(start),
            "completion_date": status.get("completionDateStruct", {}).get("date"),
            "enrollment": design.get("enrollmentInfo", {}).get("count"),
            "sponsor": sponsor.get("name"),
            "sponsor_class": sponsor.get("class"),
            "has_results": int(bool(s.get("hasResults"))),
        })
        for iv in p.get("armsInterventionsModule", {}).get("interventions", []) or []:
            inters.append({"nct_id": nct, "intervention_type": iv.get("type"),
                           "intervention_name": iv.get("name"),
                           "generic_name": match_drug(iv.get("name"))})
    t = pd.DataFrame(trials).drop_duplicates("nct_id")
    i = pd.DataFrame(inters, columns=["nct_id", "intervention_type",
                                      "intervention_name", "generic_name"]).drop_duplicates()
    return t, i


def main():
    studies = fetch_all()
    trials, inters = parse(studies)
    trials.to_csv(INTERIM / "trials.csv", index=False)
    inters.to_csv(INTERIM / "trial_interventions.csv", index=False)
    print(f"[trials] {len(trials)} studies, "
          f"{inters['generic_name'].notna().sum()} interventions matched to drug list")


if __name__ == "__main__":
    main()
