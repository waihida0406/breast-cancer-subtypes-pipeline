-- Example queries that exercise every link in the schema.

-- 1. Cohort by subtype: size, age, tumour size, disease-specific deaths
SELECT st.subtype_name, COUNT(*) AS n,
       ROUND(AVG(p.age), 1) AS mean_age,
       ROUND(AVG(p.tumor_size_mm), 1) AS mean_size_mm,
       ROUND(100.0 * AVG(p.dss_event), 1) AS pct_died_of_disease
FROM patients p JOIN subtypes st USING (subtype_code)
GROUP BY st.subtype_name ORDER BY n DESC;

-- 2. Three most frequently mutated genes in each subtype
WITH ranked AS (
  SELECT subtype_code, gene_symbol, prevalence,
         ROW_NUMBER() OVER (PARTITION BY subtype_code ORDER BY prevalence DESC) AS rk
  FROM v_mutation_prevalence_by_subtype
)
SELECT * FROM ranked WHERE rk <= 3 ORDER BY subtype_code, rk;

-- 3. Mean ERBB2 and ESR1 expression by subtype (validates the PAM50 labels)
SELECT p.subtype_code,
       ROUND(AVG(CASE WHEN g.gene_symbol = 'ERBB2' THEN e.z_score END), 2) AS erbb2_z,
       ROUND(AVG(CASE WHEN g.gene_symbol = 'ESR1'  THEN e.z_score END), 2) AS esr1_z
FROM patient_expression e
JOIN genes g    USING (gene_id)
JOIN patients p USING (patient_id)
WHERE g.gene_symbol IN ('ERBB2', 'ESR1')
GROUP BY p.subtype_code ORDER BY erbb2_z DESC;

-- 4. Discordance: patients whose clinical HER2 status disagrees with the HER2-enriched label
SELECT subtype_code, her2_status, COUNT(*) AS n
FROM patients
WHERE (subtype_code = 'HER2E' AND her2_status = 'Negative')
   OR (subtype_code <> 'HER2E' AND her2_status = 'Positive')
GROUP BY subtype_code, her2_status ORDER BY n DESC;

-- 5. Patients eligible for each modern drug class, with that class's FAERS profile
SELECT d.drug_class,
       COUNT(DISTINCT e.patient_id) AS eligible_patients,
       (SELECT SUM(d2.total_reports) FROM drugs d2
         WHERE d2.drug_class = d.drug_class) AS faers_reports
FROM drugs d LEFT JOIN patient_drug_eligibility e USING (drug_id)
GROUP BY d.drug_class ORDER BY eligible_patients DESC;

-- 6. Serious-report share by drug class
SELECT drug_class, SUM(serious_reports) AS serious, SUM(total_reports) AS total,
       ROUND(1.0 * SUM(serious_reports) / NULLIF(SUM(total_reports), 0), 3) AS serious_share
FROM drugs GROUP BY drug_class ORDER BY serious_share DESC;

-- 7. Trial activity per drug class since 2015
SELECT drug_class, SUM(n_trials) AS trials_since_2015
FROM v_trials_by_drug_class WHERE start_year >= 2015
GROUP BY drug_class ORDER BY trials_since_2015 DESC;

-- 8. HER2+ patients in METABRIC never received HER2 therapy (pre-trastuzumab era):
--    how many would be eligible today, and how did they fare?
SELECT COUNT(*) AS her2_pos_patients,
       ROUND(100.0 * AVG(dss_event), 1) AS pct_died_of_disease,
       ROUND(AVG(os_months), 1) AS mean_follow_up_months
FROM patients WHERE her2_positive = 1;
