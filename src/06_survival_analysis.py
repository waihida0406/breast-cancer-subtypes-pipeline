"""
Step 6 - Survival analysis by molecular subtype (Objective 3)

Endpoints
  OS  : overall survival, event = death from any cause
  DSS : disease-specific survival, event = death from breast cancer
        (other deaths censored; METABRIC follows patients for up to ~30 years
        and many die of other causes, so DSS is the cleaner endpoint)

Methods (implemented directly in NumPy so every formula can be shown in
the Methods chapter; results match standard packages such as lifelines):
  * Kaplan-Meier estimator with Greenwood 95% CI
  * Multi-group log-rank test
  * Cox proportional-hazards regression (Breslow ties, Newton-Raphson),
    adjusted for age, tumour size, grade and nodal status, LumA as reference

Outputs: outputs/tables/table4_survival_by_subtype.csv, table5_cox_model.csv,
         km_curves.csv (for Power BI)
         outputs/figures/fig7_km_curves.png, fig8_cox_forest.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2_dist, norm
from config import INTERIM, TAB, FIG, SUBTYPE_ORDER

COLORS = {"LumA": "#2b6cb0", "LumB": "#63b3ed", "HER2E": "#c05621",
          "Basal": "#9b2c2c", "ClaudinLow": "#b7791f", "Normal": "#718096"}


# ------------------------------------------------------------ Kaplan-Meier
def kaplan_meier(time, event):
    """Return event times, S(t), lower/upper 95% CI (Greenwood, log-log)."""
    time, event = np.asarray(time, float), np.asarray(event, int)
    ts = np.unique(time[event == 1])
    n_at_risk = np.array([(time >= t).sum() for t in ts])
    d = np.array([((time == t) & (event == 1)).sum() for t in ts])
    s = np.cumprod(1 - d / n_at_risk)
    with np.errstate(divide="ignore", invalid="ignore"):
        var = np.cumsum(d / (n_at_risk * (n_at_risk - d)))
        se_loglog = np.sqrt(var) / np.abs(np.log(s))
        lo = s ** np.exp(1.96 * se_loglog)
        hi = s ** np.exp(-1.96 * se_loglog)
    t = np.r_[0, ts]; s = np.r_[1, s]
    return t, s, np.r_[1, np.nan_to_num(lo, nan=0)], np.r_[1, np.nan_to_num(hi, nan=1)]


def surv_at(t, s, x):
    return s[np.searchsorted(t, x, side="right") - 1]


def median_survival(t, s):
    below = np.where(s <= 0.5)[0]
    return t[below[0]] if len(below) else np.nan


# --------------------------------------------------------- log-rank test
def logrank(time, event, group):
    time, event = np.asarray(time, float), np.asarray(event, int)
    labels = pd.unique(group); g = pd.Categorical(group, categories=labels).codes
    k = len(labels); O = np.zeros(k); E = np.zeros(k); V = np.zeros((k, k))
    for t in np.unique(time[event == 1]):
        at_risk = time >= t
        n_j = np.bincount(g[at_risk], minlength=k).astype(float)
        d_j = np.bincount(g[(time == t) & (event == 1)], minlength=k).astype(float)
        n, d = n_j.sum(), d_j.sum()
        O += d_j; E += d * n_j / n
        if n > 1:
            V += d * (n - d) / (n - 1) * (np.diag(n_j / n) - np.outer(n_j, n_j) / n ** 2)
    diff = (O - E)[:-1]
    stat = float(diff @ np.linalg.solve(V[:-1, :-1], diff))
    return stat, k - 1, float(chi2_dist.sf(stat, k - 1))


# --------------------------------------------------------- Cox regression
def cox_ph(X: np.ndarray, time, event, max_iter=50, tol=1e-9):
    """Cox PH with Breslow ties via Newton-Raphson. Returns beta, SE."""
    time, event = np.asarray(time, float), np.asarray(event, int)
    ev = np.where(event == 1)[0]
    R = (time[None, :] >= time[ev][:, None]).astype(float)   # risk sets
    beta = np.zeros(X.shape[1])

    def loglik_grad_hess(b):
        w = np.exp(np.clip(X @ b, -50, 50))
        Rw = R * w
        s0 = Rw.sum(1)
        s1 = Rw @ X
        s2 = np.einsum("ij,jk,jl->ikl", Rw, X, X)
        ll = float((X[ev] @ b - np.log(s0)).sum())
        grad = (X[ev] - s1 / s0[:, None]).sum(0)
        m = s1 / s0[:, None]
        hess = -(s2 / s0[:, None, None] - np.einsum("ik,il->ikl", m, m)).sum(0)
        return ll, grad, hess

    ll, g, H = loglik_grad_hess(beta)
    for _ in range(max_iter):
        step = np.linalg.solve(-H, g)
        new = beta + step
        ll_new, g_new, H_new = loglik_grad_hess(new)
        halv = 0
        while ll_new < ll - 1e-12 and halv < 20:      # step halving
            step /= 2; new = beta + step
            ll_new, g_new, H_new = loglik_grad_hess(new); halv += 1
        converged = abs(ll_new - ll) < tol
        beta, ll, g, H = new, ll_new, g_new, H_new
        if converged:
            break
    se = np.sqrt(np.diag(np.linalg.inv(-H)))
    return beta, se, ll


def cox_table(pts: pd.DataFrame) -> pd.DataFrame:
    d = pts[["os_months", "dss_event", "age", "tumor_size_mm", "grade",
             "nodes_positive", "subtype_code"]].dropna().copy()
    X = pd.DataFrame({
        "Age (per 10 years)": d.age / 10,
        "Tumour size (per 10 mm)": d.tumor_size_mm / 10,
        "Grade 3 vs 1-2": (d.grade == 3).astype(float),
        "Node-positive": (d.nodes_positive > 0).astype(float),
    })
    for s in SUBTYPE_ORDER[1:]:
        if (d.subtype_code == s).any():
            X[f"{s} vs LumA"] = (d.subtype_code == s).astype(float)
    beta, se, _ = cox_ph(X.values, d.os_months.values, d.dss_event.values)
    z = beta / se
    return pd.DataFrame({
        "variable": X.columns, "coef": beta.round(4), "se": se.round(4),
        "hazard_ratio": np.exp(beta).round(3),
        "hr_ci_lower": np.exp(beta - 1.96 * se).round(3),
        "hr_ci_upper": np.exp(beta + 1.96 * se).round(3),
        "z": z.round(2), "p_value": 2 * norm.sf(np.abs(z)),
        "n": len(d), "events": int(d.dss_event.sum())})


# ------------------------------------------------------------------- main
def main():
    pts = pd.read_csv(INTERIM / "patients_clean.csv")
    order = [s for s in SUBTYPE_ORDER if s in set(pts.subtype_code)]
    endpoints = {"OS": "os_event", "DSS": "dss_event"}

    rows, curves = [], []
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for ax, (ep, col) in zip(axes, endpoints.items()):
        stat, dof, p = logrank(pts.os_months, pts[col], pts.subtype_code)
        for s in order:
            g = pts[pts.subtype_code == s]
            t, sv, lo, hi = kaplan_meier(g.os_months, g[col])
            ax.step(t / 12, sv, where="post", color=COLORS[s], label=f"{s} (n={len(g)})", lw=1.6)
            ax.fill_between(t / 12, lo, hi, step="post", color=COLORS[s], alpha=0.08)
            rows.append({"endpoint": ep, "subtype_code": s, "n": len(g),
                         "events": int(g[col].sum()),
                         "median_months": median_survival(t, sv),
                         "surv_5y": round(surv_at(t, sv, 60), 3),
                         "surv_10y": round(surv_at(t, sv, 120), 3),
                         "logrank_chi2": round(stat, 2), "logrank_dof": dof,
                         "logrank_p": p})
            curves.append(pd.DataFrame({"endpoint": ep, "subtype_code": s,
                                        "months": t, "survival": sv,
                                        "ci_lower": lo, "ci_upper": hi}))
        ax.set_title(f"{'Overall' if ep == 'OS' else 'Disease-specific'} survival\n"
                     f"log-rank chi2 = {stat:.1f}, df = {dof}, p = {p:.1e}", fontsize=10)
        ax.set_xlabel("Years since diagnosis"); ax.set_ylim(0, 1.02); ax.set_xlim(0, 25)
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("Survival probability"); axes[1].legend(fontsize=8, loc="lower left")
    fig.suptitle("Figure 7. Kaplan-Meier survival by molecular subtype")
    fig.tight_layout(); fig.savefig(FIG / "fig7_km_curves.png", dpi=200); plt.close(fig)

    surv = pd.DataFrame(rows)
    surv.to_csv(TAB / "table4_survival_by_subtype.csv", index=False)
    pd.concat(curves).round(4).to_csv(TAB / "km_curves.csv", index=False)

    cox = cox_table(pts)
    cox.to_csv(TAB / "table5_cox_model.csv", index=False)
    c = cox.iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.5, 0.42 * len(c) + 1.4))
    y = np.arange(len(c))
    ax.errorbar(c.hazard_ratio, y, xerr=[c.hazard_ratio - c.hr_ci_lower, c.hr_ci_upper - c.hazard_ratio],
                fmt="s", color="#1a365d", ecolor="#4a5568", capsize=3, ms=5)
    ax.axvline(1, ls="--", color="grey", lw=0.8)
    ax.set_xscale("log"); ax.set_yticks(y, c.variable)
    lo_lim, hi_lim = c.hr_ci_lower.min() * 0.9, c.hr_ci_upper.max() * 1.1
    ticks = [t for t in [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8] if lo_lim <= t <= hi_lim] or [1]
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks]); ax.minorticks_off()
    for yi, (hr, lo, hi, p) in enumerate(zip(c.hazard_ratio, c.hr_ci_lower, c.hr_ci_upper, c.p_value)):
        ax.text(1.02, yi, f"{hr:.2f} ({lo:.2f}-{hi:.2f})  p={p:.1e}",
                transform=ax.get_yaxis_transform(), va="center", fontsize=7.5)
    ax.set_xlabel("Hazard ratio for breast-cancer death (log scale, 95% CI)")
    ax.set_title(f"Figure 8. Multivariable Cox model (n={cox.n[0]}, events={cox.events[0]})")
    fig.tight_layout(); fig.savefig(FIG / "fig8_cox_forest.png", dpi=200, bbox_inches="tight"); plt.close(fig)

    print("[step 6]", surv[surv.endpoint == "DSS"][["subtype_code", "n", "events", "surv_10y"]]
          .to_string(index=False))
    print(f"[step 6] DSS log-rank p = {surv[surv.endpoint == 'DSS'].logrank_p.iloc[0]:.2e}")
    print(cox[["variable", "hazard_ratio", "hr_ci_lower", "hr_ci_upper", "p_value"]].to_string(index=False))


if __name__ == "__main__":
    main()
