-- Breast cancer integrative database (SQLite; portable to PostgreSQL / SQL Server)
--
-- Link paths
--   patients -> subtypes
--   patients -> patient_mutations / patient_expression -> genes
--   patients -> patient_treatments -> drugs (by therapy type, METABRIC era)
--   patients -> patient_drug_eligibility -> drugs (by receptor / mutation status)
--   drugs -> adverse_events, fda_yearly          (openFDA FAERS)
--   drugs -> trial_interventions -> trials       (ClinicalTrials.gov)

PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS trial_interventions;
DROP TABLE IF EXISTS trials;
DROP TABLE IF EXISTS fda_yearly;
DROP TABLE IF EXISTS adverse_events;
DROP TABLE IF EXISTS patient_drug_eligibility;
DROP TABLE IF EXISTS patient_treatments;
DROP TABLE IF EXISTS drugs;
DROP TABLE IF EXISTS patient_expression;
DROP TABLE IF EXISTS patient_mutations;
DROP TABLE IF EXISTS genes;
DROP TABLE IF EXISTS patients;
DROP TABLE IF EXISTS subtypes;

CREATE TABLE subtypes (
    subtype_code TEXT PRIMARY KEY,          -- LumA, LumB, HER2E, Basal, ClaudinLow, Normal
    subtype_name TEXT NOT NULL UNIQUE,
    description  TEXT
);

CREATE TABLE patients (
    patient_id        TEXT PRIMARY KEY,
    subtype_code      TEXT NOT NULL REFERENCES subtypes(subtype_code),
    age               REAL NOT NULL CHECK (age BETWEEN 18 AND 100),
    age_band          TEXT,
    menopausal_state  TEXT,
    histology         TEXT,
    cellularity       TEXT,
    tumor_size_mm     REAL,
    tumor_stage       INTEGER,
    grade             INTEGER,
    nodes_positive    INTEGER,
    npi               REAL,
    er_status         TEXT CHECK (er_status IN ('Positive','Negative') OR er_status IS NULL),
    pr_status         TEXT CHECK (pr_status IN ('Positive','Negative') OR pr_status IS NULL),
    her2_status       TEXT CHECK (her2_status IN ('Positive','Negative') OR her2_status IS NULL),
    hr_positive       INTEGER NOT NULL,
    her2_positive     INTEGER NOT NULL,
    triple_negative   INTEGER NOT NULL,
    pik3ca_mut        INTEGER NOT NULL DEFAULT 0,
    brca_mut          INTEGER NOT NULL DEFAULT 0,
    integrative_cluster TEXT,
    mutation_count    INTEGER,
    surgery           TEXT,
    os_months         REAL NOT NULL,
    os_event          INTEGER NOT NULL,      -- 1 = died (any cause)
    dss_event         INTEGER NOT NULL,      -- 1 = died of breast cancer
    vital_status      TEXT,
    cohort            INTEGER
);

CREATE TABLE genes (
    gene_id     INTEGER PRIMARY KEY,
    gene_symbol TEXT NOT NULL UNIQUE,        -- upper-case HGNC symbol
    in_expression_panel INTEGER NOT NULL DEFAULT 0,
    in_mutation_panel   INTEGER NOT NULL DEFAULT 0
);

-- Long (tidy) bridge: one row per patient x mutated gene (present only)
CREATE TABLE patient_mutations (
    patient_id TEXT    NOT NULL REFERENCES patients(patient_id),
    gene_id    INTEGER NOT NULL REFERENCES genes(gene_id),
    PRIMARY KEY (patient_id, gene_id)
);

-- Long table of mRNA z-scores
CREATE TABLE patient_expression (
    patient_id TEXT    NOT NULL REFERENCES patients(patient_id),
    gene_id    INTEGER NOT NULL REFERENCES genes(gene_id),
    z_score    REAL,
    PRIMARY KEY (patient_id, gene_id)
);

CREATE TABLE drugs (
    drug_id          INTEGER PRIMARY KEY,
    generic_name     TEXT NOT NULL UNIQUE,
    drug_class       TEXT NOT NULL,
    target_group     TEXT NOT NULL,          -- HR+, HER2+, TNBC, PIK3CA-mut ..., Any
    metabric_therapy TEXT,                   -- therapy flag in METABRIC it belongs to
    total_reports    INTEGER,
    serious_reports  INTEGER,
    death_reports    INTEGER
);

-- Treatments actually recorded in METABRIC (yes/no by therapy type)
CREATE TABLE patient_treatments (
    patient_id TEXT NOT NULL REFERENCES patients(patient_id),
    therapy    TEXT NOT NULL CHECK (therapy IN ('Chemotherapy','Hormone therapy','Radiotherapy')),
    PRIMARY KEY (patient_id, therapy)
);

-- Which current drugs each patient's tumour profile would make them a candidate for
CREATE TABLE patient_drug_eligibility (
    patient_id TEXT    NOT NULL REFERENCES patients(patient_id),
    drug_id    INTEGER NOT NULL REFERENCES drugs(drug_id),
    PRIMARY KEY (patient_id, drug_id)
);

CREATE TABLE adverse_events (
    drug_id      INTEGER NOT NULL REFERENCES drugs(drug_id),
    reaction     TEXT    NOT NULL,
    report_count INTEGER NOT NULL,
    PRIMARY KEY (drug_id, reaction)
);

CREATE TABLE fda_yearly (
    drug_id      INTEGER NOT NULL REFERENCES drugs(drug_id),
    year         INTEGER NOT NULL,
    report_count INTEGER NOT NULL,
    PRIMARY KEY (drug_id, year)
);

CREATE TABLE trials (
    nct_id          TEXT PRIMARY KEY,
    title           TEXT,
    status          TEXT,
    study_type      TEXT,
    phase           TEXT,
    start_date      TEXT,
    start_year      INTEGER,
    completion_date TEXT,
    enrollment      INTEGER,
    sponsor         TEXT,
    sponsor_class   TEXT,
    has_results     INTEGER
);

CREATE TABLE trial_interventions (
    nct_id            TEXT NOT NULL REFERENCES trials(nct_id),
    intervention_type TEXT,
    intervention_name TEXT,
    drug_id           INTEGER REFERENCES drugs(drug_id)   -- NULL if not in the drug list
);

CREATE INDEX ix_patients_subtype ON patients(subtype_code);
CREATE INDEX ix_pm_gene          ON patient_mutations(gene_id);
CREATE INDEX ix_pe_gene          ON patient_expression(gene_id);
CREATE INDEX ix_elig_drug        ON patient_drug_eligibility(drug_id);
CREATE INDEX ix_ti_drug          ON trial_interventions(drug_id);
CREATE INDEX ix_trials_year      ON trials(start_year);
