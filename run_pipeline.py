"""
Run the complete project pipeline.

    python run_pipeline.py              # all steps, including API downloads
    python run_pipeline.py --skip-fetch # reuse cached API data (offline)
"""
import runpy
import sys
from pathlib import Path

SRC = Path(__file__).parent / "src"
sys.path.insert(0, str(SRC))

STEPS = [
    "00_get_metabric.py",        # skipped if data/raw/METABRIC_RNA_Mutation.csv exists
    "01_data_cleaning.py",
    "02_fetch_openfda.py",
    "03_fetch_clinicaltrials.py",
    "04_build_database.py",
    "05_subtype_associations.py",
    "06_survival_analysis.py",
    "07_subtype_classifier.py",
    "08_drug_safety_trials.py",
    "09_export_powerbi.py",
]

if __name__ == "__main__":
    skip = "--skip-fetch" in sys.argv
    for step in STEPS:
        if skip and "fetch" in step:
            print(f"--- skipping {step} (using cached data)")
            continue
        print(f"\n=== {step} ===")
        runpy.run_path(str(SRC / step), run_name="__main__")
    print("\nAll steps complete. Figures: outputs/figures | Tables: outputs/tables "
          "| Power BI: powerbi/data")
