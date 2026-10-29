# Project report outline

Chapter structure for the project report, with each figure and table mapped to the chapter it belongs in. After running the pipeline, figures are in `outputs/figures/` and tables in `outputs/tables/`.

## Front matter
Title page · Declaration of originality · Acknowledgement · Abstract (≈250 words, from README) · Table of contents · List of figures · List of tables · Abbreviations (PAM50, ER, PR, HER2, HR, TNBC, IHC, NPI, OS, DSS, HR (hazard ratio), CI, FAERS, MedDRA, ADC, CDK4/6, PARP, SMOTE, FDR, CV)

## Chapter 1 — Introduction (5–7 pages)
- 1.1 Breast cancer burden, globally and in India [Ref 1]
- 1.2 Clinical classification: ER, PR, HER2 by IHC; grade; Nottingham Prognostic Index [Ref 9]
- 1.3 Molecular subtypes: PAM50 and claudin-low [Ref 2, 3, 4]
- 1.4 HER2 biology and HER2-targeted therapy [Ref 7, 10]
- 1.5 Treatment by subtype: endocrine, CDK4/6, PI3K, HER2 agents, PARP, immunotherapy, ADCs
- 1.6 Real-world data: METABRIC, cBioPortal, FAERS, trial registries [Ref 5, 6, 11, 12]
- 1.7 Rationale, aim, objectives, research questions

## Chapter 2 — Review of literature (8–10 pages)
- Development and prognostic value of intrinsic subtypes
- METABRIC findings: integrative clusters, mutation landscape
- Concordance between IHC surrogates and PAM50 [Ref 8]
- Machine learning for breast cancer subtype prediction
- Survival differences between subtypes
- FAERS pharmacovigilance studies of HER2 agents, CDK4/6 inhibitors and ADCs
- Breast cancer trial landscape analyses
- Research gap: few studies connect subtype biology, outcomes, drug safety and trial activity in one framework

## Chapter 3 — Materials and methods
- 3.1 Workflow — **Fig 0** (`docs/methodology_workflow.png`)
- 3.2 Data sources and variables (data dictionary table); METABRIC retrieval from cBioPortal and cohort selection (Step 0: 2,509 → 1,904 patients with expression and sequencing → 1,898 after removing 6 unclassified)
- 3.3 Data cleaning and derived variables (Step 1) — `cleaning_report.csv`
- 3.4 API retrieval: endpoints, query fields, indication filter, caching (Steps 2–3)
- 3.5 Database schema and eligibility rules — draw an ER diagram from `sql/schema.sql` (dbdiagram.io or draw.io)
- 3.6 Statistical analysis: chi-square, permutation test, BH-FDR, Cramér's V, Kruskal–Wallis (Step 5) [Ref 16, 17]
- 3.7 Survival analysis: KM, log-rank, Cox PH, with the formulas (Step 6) [Ref 13, 14, 15]
- 3.8 Machine learning: feature sets, split, CV design, imbalance handling, metrics (Step 7) [Ref 18, 19]
- 3.9 Drug-safety and trial-trend analysis (Step 8)
- 3.10 Dashboard development (Step 9, `powerbi/DASHBOARD_SPEC.md`)
- 3.11 Tools: Python (pandas, NumPy, SciPy, scikit-learn, matplotlib), SQLite, Power BI

## Chapter 4 — Results
| Section | Figures / tables |
|---|---|
| 4.1 Data quality and cohort | Table 1, Fig 1, 2 |
| 4.2 Receptor status vs molecular subtype | Fig 3, SQL query 4 (discordance) |
| 4.3 Mutation landscape | Fig 4 |
| 4.4 Marker gene expression (label validation) | Fig 5 |
| 4.5 Statistical associations | Tables 2–3, Fig 6 |
| 4.6 Survival by subtype | Table 4, Fig 7 |
| 4.7 Multivariable Cox model | Table 5, Fig 8 |
| 4.8 Subtype prediction | Tables 6–7, Fig 9, 10 |
| 4.9 Predictive features | Fig 11 |
| 4.10 Drug eligibility and safety | Table 8, Fig 12, `v_eligibility_by_subtype` |
| 4.11 Clinical-trial trends | Table 9, Fig 13, 14 |
| 4.12 Dashboard | Screenshots of the 3 Power BI pages |

Describe what each figure shows here; keep interpretation for Chapter 5.

## Chapter 5 — Discussion
- Do mutation patterns match published METABRIC findings (TP53 in Basal/HER2E; PIK3CA, GATA3, MAP3K1, CDH1 in luminal tumours)?
- Does marker expression (Fig 5) confirm the subtype labels? ERBB2/GRB7 should be high in HER2E and KRT5/KRT17/FOXC1 in Basal.
- Which subtypes are worst after adjustment (Fig 8), and why HER2E risk is overstated in this pre-trastuzumab cohort
- What the clinical-vs-expression gap (Fig 9) means for hospitals that report only IHC; why LumA vs LumB and Normal-like are hardest to separate (Ki-67 and proliferation are not in routine data)
- Discordant cases: HER2E tumours that are HER2-negative by IHC/SNP6
- FAERS results in context: reporting bias, market size and time on market; the indication filter
- Trial trends: the shift toward ADCs and CDK4/6 inhibitors
- Relevance to Indian patients (younger age, more TNBC) and limitations

## Chapter 6 — Summary and conclusion
One paragraph per objective, answering the matching research question, followed by future work.

## References
Vancouver style. The README list is the starting set; add the papers you cite in Chapter 2.

## Appendix
- A: SQL schema and key queries (`sql/`)
- B: Data dictionary
- C: Eligibility rules table
- D: Code repository link

## Likely viva questions
- What is PAM50, and why does claudin-low sit outside it?
- Why do IHC subtypes and PAM50 subtypes disagree for some patients?
- Why use disease-specific survival rather than overall survival in METABRIC?
- What does a hazard ratio of 2 mean? What is the proportional-hazards assumption, and how would you check it?
- Why use a log-rank test rather than comparing 5-year survival percentages?
- Why a permutation p-value? When does the chi-square approximation fail?
- Why correct for multiple testing, and why FDR rather than Bonferroni?
- Why is accuracy misleading for these subtypes? What is macro-F1?
- Why must SMOTE, imputation and scaling be fitted inside each CV fold?
- Isn't predicting PAM50 from expression circular? What does the clinical model add?
- Can FAERS show that one drug class is more toxic than another? Why not?
- Why are HER2-enriched patients' outcomes in METABRIC worse than you would expect today?
