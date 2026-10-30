-- Analytical views (used by Steps 5 and 8 and exported for Power BI)

DROP VIEW IF EXISTS v_subtype_size;
CREATE VIEW v_subtype_size AS
SELECT subtype_code, COUNT(*) AS n_patients FROM patients GROUP BY subtype_code;

-- Share of each subtype carrying a mutation in each gene
DROP VIEW IF EXISTS v_mutation_prevalence_by_subtype;
CREATE VIEW v_mutation_prevalence_by_subtype AS
SELECT s.subtype_code, g.gene_symbol,
       COUNT(pm.patient_id)                                AS n_mutated,
       s.n_patients,
       ROUND(1.0 * COUNT(pm.patient_id) / s.n_patients, 4) AS prevalence
FROM v_subtype_size s
CROSS JOIN genes g
LEFT JOIN patients p          ON p.subtype_code = s.subtype_code
LEFT JOIN patient_mutations pm ON pm.patient_id = p.patient_id AND pm.gene_id = g.gene_id
WHERE g.in_mutation_panel = 1
GROUP BY s.subtype_code, g.gene_id;

-- Receptor profile (IHC / SNP6) per molecular subtype
DROP VIEW IF EXISTS v_receptor_profile_by_subtype;
CREATE VIEW v_receptor_profile_by_subtype AS
SELECT subtype_code, COUNT(*) AS n,
       ROUND(AVG(er_status = 'Positive'), 4) AS er_pos,
       ROUND(AVG(pr_status = 'Positive'), 4) AS pr_pos,
       ROUND(AVG(her2_positive), 4)          AS her2_pos,
       ROUND(AVG(triple_negative), 4)        AS tnbc,
       ROUND(AVG(dss_event), 4)              AS died_of_disease
FROM patients GROUP BY subtype_code;

-- Treatments received in the METABRIC cohort, by subtype
DROP VIEW IF EXISTS v_treatment_by_subtype;
CREATE VIEW v_treatment_by_subtype AS
SELECT s.subtype_code, t.therapy, COUNT(pt.patient_id) AS n_treated, s.n_patients,
       ROUND(1.0 * COUNT(pt.patient_id) / s.n_patients, 4) AS share_treated
FROM v_subtype_size s
CROSS JOIN (SELECT 'Chemotherapy' AS therapy UNION ALL SELECT 'Hormone therapy'
            UNION ALL SELECT 'Radiotherapy') t
LEFT JOIN patients p            ON p.subtype_code = s.subtype_code
LEFT JOIN patient_treatments pt ON pt.patient_id = p.patient_id AND pt.therapy = t.therapy
GROUP BY s.subtype_code, t.therapy;

-- Share of each subtype that would be a candidate for each drug class today
DROP VIEW IF EXISTS v_eligibility_by_subtype;
CREATE VIEW v_eligibility_by_subtype AS
SELECT s.subtype_code, d.drug_class,
       COUNT(DISTINCT e.patient_id) AS n_eligible, s.n_patients,
       ROUND(1.0 * COUNT(DISTINCT e.patient_id) / s.n_patients, 4) AS share_eligible
FROM v_subtype_size s
CROSS JOIN (SELECT DISTINCT drug_class FROM drugs) d
LEFT JOIN patients p ON p.subtype_code = s.subtype_code
LEFT JOIN patient_drug_eligibility e ON e.patient_id = p.patient_id
     AND e.drug_id IN (SELECT drug_id FROM drugs WHERE drug_class = d.drug_class)
GROUP BY s.subtype_code, d.drug_class;

DROP VIEW IF EXISTS v_drug_safety;
CREATE VIEW v_drug_safety AS
SELECT d.drug_id, d.generic_name, d.drug_class, d.target_group,
       d.total_reports, d.serious_reports, d.death_reports,
       ROUND(1.0 * d.serious_reports / NULLIF(d.total_reports, 0), 4) AS serious_share,
       (SELECT COUNT(DISTINCT ti.nct_id) FROM trial_interventions ti
         WHERE ti.drug_id = d.drug_id)                                AS n_trials,
       (SELECT COUNT(*) FROM patient_drug_eligibility e
         WHERE e.drug_id = d.drug_id)                                 AS n_eligible_patients
FROM drugs d;

-- Reporting share of each reaction within a drug (not an incidence rate)
DROP VIEW IF EXISTS v_top_reactions;
CREATE VIEW v_top_reactions AS
SELECT d.generic_name, d.drug_class, ae.reaction, ae.report_count,
       ROUND(1.0 * ae.report_count / NULLIF(d.total_reports, 0), 4) AS share_of_reports
FROM adverse_events ae JOIN drugs d USING (drug_id);

DROP VIEW IF EXISTS v_trials_by_year;
CREATE VIEW v_trials_by_year AS
SELECT start_year, phase, status, COUNT(*) AS n_trials,
       SUM(COALESCE(enrollment, 0)) AS total_enrollment
FROM trials WHERE start_year IS NOT NULL
GROUP BY start_year, phase, status;

DROP VIEW IF EXISTS v_trials_by_drug_class;
CREATE VIEW v_trials_by_drug_class AS
SELECT d.drug_class, t.start_year, COUNT(DISTINCT t.nct_id) AS n_trials
FROM trial_interventions ti
JOIN drugs d  USING (drug_id)
JOIN trials t USING (nct_id)
GROUP BY d.drug_class, t.start_year;
