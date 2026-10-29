# Power BI dashboard build guide (Objective 6)

Load every CSV in `powerbi/data/` (Get Data → Folder, or Text/CSV one by one).

## Data model (relationships, single direction 1 → *)

| From (1) | To (*) | Key |
|---|---|---|
| dim_subtype | fact_patient | subtype_code |
| dim_subtype | analysis_km_curves | subtype_code |
| dim_subtype | analysis_table4_survival_by_subtype | subtype_code |
| fact_patient | fact_patient_mutation | patient_id |
| fact_patient | fact_marker_expression | patient_id |
| fact_patient | fact_patient_treatment | patient_id |
| fact_patient | bridge_patient_drug_eligibility | patient_id |
| dim_gene | fact_patient_mutation | gene_id |
| dim_gene | fact_marker_expression | gene_id |
| dim_drug | bridge_patient_drug_eligibility | drug_id |
| dim_drug | fact_drug_safety | drug_id |
| dim_drug | fact_adverse_event | drug_id |
| dim_drug | fact_fda_yearly | drug_id |
| dim_drug | bridge_trial_drug | drug_id |
| dim_trial | bridge_trial_drug | nct_id |

Set the cross-filter direction to **Both** on `bridge_patient_drug_eligibility` → `dim_drug`, so a subtype slicer filters drug visuals. Create a `dim_year` table and relate it to `fact_fda_yearly[year]` and `dim_trial[start_year]`.

## DAX measures

```DAX
dim_year = GENERATESERIES(1995, YEAR(TODAY()), 1)          -- rename column to Year

Patients            = DISTINCTCOUNT(fact_patient[patient_id])
Median Age          = MEDIAN(fact_patient[age])
% HER2+             = DIVIDE(SUM(fact_patient[her2_positive]), [Patients])
% Triple-negative   = DIVIDE(SUM(fact_patient[triple_negative]), [Patients])
Deaths from Disease = SUM(fact_patient[dss_event])
% Died of Disease   = DIVIDE([Deaths from Disease], [Patients])
Mutated Patients    = DISTINCTCOUNT(fact_patient_mutation[patient_id])
Mutation Prevalence = DIVIDE([Mutated Patients], [Patients])
Mean Expression     = AVERAGE(fact_marker_expression[z_score])
Eligible Patients   = DISTINCTCOUNT(bridge_patient_drug_eligibility[patient_id])
% Eligible          = DIVIDE([Eligible Patients], CALCULATE([Patients], ALL(dim_drug)))

FAERS Reports       = SUM(fact_drug_safety[total_reports])
Serious Reports     = SUM(fact_drug_safety[serious_reports])
Serious Share       = DIVIDE([Serious Reports], [FAERS Reports])
Reaction Reports    = SUM(fact_adverse_event[report_count])
Yearly Reports      = SUM(fact_fda_yearly[report_count])

Trials              = DISTINCTCOUNT(dim_trial[nct_id])
Trials Recruiting   = CALCULATE([Trials], dim_trial[status] = "RECRUITING")
Trials w/ Results % = DIVIDE(CALCULATE([Trials], dim_trial[has_results] = 1), [Trials])
Total Enrollment    = SUM(dim_trial[enrollment])
```

## Page 1 — Patients, Genomics & Survival
- KPI cards: Patients, Median Age, % HER2+, % Triple-negative, % Died of Disease
- Column chart: Patients by `subtype_name` (shows the class imbalance)
- Clustered bar: ER+/PR+/HER2+/TNBC share by subtype (from `fact_patient`)
- Matrix heatmap: rows `dim_gene[gene_symbol]` (top 15 by Mutated Patients), columns subtype, value Mutation Prevalence, background colour scale
- Matrix heatmap: marker genes × subtype, value Mean Expression, diverging colour scale (blue → white → red)
- Line chart: `analysis_km_curves` with x = months, y = survival, legend = subtype_code, filter endpoint = "DSS" (set the line to stepped)
- Table: `analysis_table5_cox_model` (variable, hazard_ratio, CI, p_value)
- Slicers: subtype, age band, endpoint

## Page 2 — Drug Safety (openFDA FAERS)
- KPI cards: FAERS Reports, Serious Share, death reports
- Bar: Serious Share by `drug_class`, drill down to `generic_name`
- Matrix: subtype × drug_class, value % Eligible (shows which subtypes each drug class serves)
- Treemap: top 15 `reaction` by Reaction Reports, filtered by a drug slicer
- Line: Yearly Reports by `dim_year[Year]`, legend `drug_class`
- Footnote text box: "FAERS counts are spontaneous reports with a breast-cancer indication. They depend on market size and time on market, and are not incidence rates."

## Page 3 — Clinical Trial Trends (ClinicalTrials.gov)
- KPI cards: Trials, Trials Recruiting, Total Enrollment, Trials w/ Results %
- Stacked column: Trials by `dim_year[Year]`, legend `phase`
- Line: Trials by year, legend `dim_drug[drug_class]` (CDK4/6 and ADC growth)
- Donut: Trials by `status`; bar: Trials by `sponsor_class`
- Scatter: x = Trials, y = FAERS Reports, size = Eligible Patients, details = `generic_name`
- Table: nct_id, title, phase, status, start_year (drill-through target)

Theme: one colour per subtype, used on all pages (LumA dark blue, LumB light blue, HER2E orange, Basal dark red, Claudin-low ochre, Normal grey); red only for serious/death measures.
