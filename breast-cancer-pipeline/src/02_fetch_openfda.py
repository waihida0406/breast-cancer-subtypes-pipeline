"""
Step 2 (Objective 4) - Adverse-event reports for breast cancer drugs (openFDA FAERS)

25 drugs across endocrine therapy, CDK4/6 and PI3K inhibitors, HER2-targeted
agents, PARP inhibitors, immunotherapy, TROP2 ADC and chemotherapy.
With FDA_BREAST_ONLY, a report is counted only if it lists a breast-cancer
indication, so multi-cancer drugs (pembrolizumab, taxanes) are comparable.

Endpoint: https://api.fda.gov/drug/event.json
Outputs (data/interim/):
  fda_drug_summary.csv   one row per drug: total, serious, death reports
  fda_reactions.csv      top MedDRA reactions per drug
  fda_reports_by_year.csv report counts per drug per year
Raw JSON responses are cached in data/raw/openfda/ so re-runs are free.

FAERS is spontaneous reporting: counts reflect reporting, not incidence,
and cannot be read as risk rates without a denominator.
"""
import json
import time
import requests
import pandas as pd
from config import RAW, INTERIM, OPENFDA_KEY, DRUGS, FDA_BREAST_ONLY

BASE = "https://api.fda.gov/drug/event.json"
CACHE = RAW / "openfda"
CACHE.mkdir(parents=True, exist_ok=True)
PAUSE = 0.3  # stays well under 240 requests/minute


def _get(params: dict, cache_name: str) -> dict:
    path = CACHE / f"{cache_name}.json"
    if path.exists():
        return json.loads(path.read_text())
    if OPENFDA_KEY:
        params = {**params, "api_key": OPENFDA_KEY}
    for attempt in range(4):
        r = requests.get(BASE, params=params, timeout=30)
        if r.status_code == 404:          # openFDA returns 404 for "no matches"
            data = {"results": [], "meta": {"results": {"total": 0}}}
            break
        if r.status_code == 429:
            time.sleep(2 ** attempt * 5)
            continue
        r.raise_for_status()
        data = r.json()
        break
    else:
        raise RuntimeError(f"openFDA rate-limited on {cache_name}")
    path.write_text(json.dumps(data))
    time.sleep(PAUSE)
    return data


BREAST = 'patient.drug.drugindication:"breast"'


def drug_query(generic: str, extra: str = "", override: str | None = None) -> str:
    """Build an openFDA search string. Spaces (not '+') join terms because
    requests URL-encodes '+'; the API decodes spaces correctly."""
    parts = [override or f'patient.drug.openfda.generic_name:"{generic}"']
    if FDA_BREAST_ONLY:
        parts.append(BREAST)
    if extra:
        parts.append(extra)
    return " AND ".join(parts)


def total(data: dict) -> int:
    return int(data.get("meta", {}).get("results", {}).get("total", 0))


def parse_counts(data: dict, key_name: str) -> pd.DataFrame:
    # count=<field>.exact returns {"term","count"}; count=receivedate returns {"time","count"}
    return pd.DataFrame(data.get("results", [])).rename(
        columns={"term": key_name, "time": key_name, "count": "report_count"})


def fetch_drug(generic: str, override=None) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    slug = generic.replace(" ", "_")
    q = lambda extra="": drug_query(generic, extra, override)
    all_ = _get({"search": q(), "limit": 1}, f"{slug}_total")
    ser = _get({"search": q("serious:1"), "limit": 1}, f"{slug}_serious")
    death = _get({"search": q("seriousnessdeath:1"), "limit": 1}, f"{slug}_death")
    rx = _get({"search": q(),
               "count": "patient.reaction.reactionmeddrapt.exact", "limit": 50},
              f"{slug}_reactions")
    yr = _get({"search": q(), "count": "receivedate"}, f"{slug}_dates")

    summary = {"generic_name": generic, "total_reports": total(all_),
               "serious_reports": total(ser), "death_reports": total(death)}

    reactions = parse_counts(rx, "reaction")
    reactions["generic_name"] = generic
    reactions["reaction"] = reactions.get("reaction", pd.Series(dtype=str)).str.title()

    by_day = parse_counts(yr, "receivedate")
    if not by_day.empty:
        by_day["year"] = by_day["receivedate"].astype(str).str[:4].astype(int)
        by_year = by_day.groupby("year", as_index=False)["report_count"].sum()
    else:
        by_year = pd.DataFrame(columns=["year", "report_count"])
    by_year["generic_name"] = generic
    return summary, reactions, by_year


def main():
    summaries, reactions, years = [], [], []
    for generic, _cls, _target, _therapy, override in DRUGS:
        try:
            s, r, y = fetch_drug(generic, override)
        except Exception as e:  # keep going if one drug fails
            print(f"[openFDA] {generic}: FAILED ({e})")
            continue
        summaries.append(s); reactions.append(r); years.append(y)
        print(f"[openFDA] {generic}: {s['total_reports']:,} reports")

    pd.DataFrame(summaries).to_csv(INTERIM / "fda_drug_summary.csv", index=False)
    pd.concat(reactions).to_csv(INTERIM / "fda_reactions.csv", index=False)
    pd.concat(years).to_csv(INTERIM / "fda_reports_by_year.csv", index=False)


if __name__ == "__main__":
    main()
