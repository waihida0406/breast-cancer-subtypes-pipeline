"""
Step 8 - Drug safety and clinical-trial trend analysis (Objective 5)

All numbers here come from SQL queries on the integrated database, which
shows the value of linking the sources in Step 4.

Outputs: outputs/tables/table8_drug_safety.csv, table9_trials_by_class.csv
         outputs/figures/fig12_serious_share_by_class.png
         outputs/figures/fig13_trials_over_time.png
         outputs/figures/fig14_trials_by_drug_class.png

Interpretation note for the Discussion chapter: FAERS is a spontaneous
reporting system. Counts reflect reporting behaviour and market size,
not incidence, so they cannot be compared as risk rates.
"""
import sqlite3
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from config import DB_PATH, TAB, FIG

# Light-to-dark blues with a visible lightest shade
SHADES = ["#90cdf4", "#4299e1", "#2b6cb0", "#1a365d", "#63b3ed", "#2c5282"]


def main():
    con = sqlite3.connect(DB_PATH)

    # Table 8 - safety profile per drug (SQL view v_drug_safety)
    safety = pd.read_sql("""
        SELECT generic_name, drug_class, target_group, total_reports,
               serious_reports, death_reports, serious_share, n_trials,
               n_eligible_patients
        FROM v_drug_safety ORDER BY drug_class, total_reports DESC""", con)
    safety.to_csv(TAB / "table8_drug_safety.csv", index=False)

    by_class = pd.read_sql("""
        SELECT drug_class,
               1.0 * SUM(serious_reports) / NULLIF(SUM(total_reports), 0) AS serious_share,
               SUM(total_reports) AS total_reports
        FROM drugs GROUP BY drug_class ORDER BY serious_share""", con)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(by_class.drug_class, by_class.serious_share, color="#c53030")
    for i, (s, n) in enumerate(zip(by_class.serious_share, by_class.total_reports)):
        ax.text(s + 0.01, i, f"{s:.0%}  (n={n:,})", va="center", fontsize=7)
    ax.set_xlim(0, 1); ax.set_xlabel("Share of FAERS reports classed as serious")
    ax.set_title("Figure 12. Serious adverse-event report share by drug class\n(FAERS reports with a breast-cancer indication)")
    fig.tight_layout(); fig.savefig(FIG / "fig12_serious_share_by_class.png", dpi=200); plt.close(fig)

    # Trials over time by phase
    trials = pd.read_sql("""
        SELECT start_year, phase, COUNT(*) AS n FROM trials
        WHERE start_year BETWEEN 2000 AND strftime('%Y','now')
        GROUP BY start_year, phase""", con)
    piv = trials.pivot(index="start_year", columns="phase", values="n").fillna(0)
    fig, ax = plt.subplots(figsize=(8, 4))
    piv.plot(kind="bar", stacked=True, ax=ax, width=0.85, color=SHADES[:piv.shape[1]])
    ax.set_xlabel("Start year"); ax.set_ylabel("Number of trials")
    ax.set_title("Figure 13. Breast cancer clinical trials started per year, by phase")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(FIG / "fig13_trials_over_time.png", dpi=200); plt.close(fig)

    # Table 9 / Fig 14 - trial activity per drug class over time
    cls = pd.read_sql("SELECT * FROM v_trials_by_drug_class", con)
    cls.to_csv(TAB / "table9_trials_by_class.csv", index=False)
    cls["period"] = pd.cut(cls.start_year, [1990, 2009, 2014, 2019, 2100],
                           labels=["<=2009", "2010-14", "2015-19", "2020+"])
    pv = cls.groupby(["drug_class", "period"], observed=False)["n_trials"].sum().unstack()
    fig, ax = plt.subplots(figsize=(10, 4.5))
    pv.plot(kind="bar", ax=ax, color=SHADES[:pv.shape[1]])
    ax.set_xlabel(""); ax.set_ylabel("Number of trials"); ax.tick_params(axis="x", rotation=35, labelsize=8)
    plt.setp(ax.get_xticklabels(), ha="right")
    ax.set_title("Figure 14. Trials per drug class across time periods")
    fig.tight_layout(); fig.savefig(FIG / "fig14_trials_by_drug_class.png", dpi=200); plt.close(fig)

    con.close()
    print(f"[step 8] {len(safety)} drugs profiled; figures 12-14 saved")


if __name__ == "__main__":
    main()
