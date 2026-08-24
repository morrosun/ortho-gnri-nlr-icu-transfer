-- Merge all three databases into a single dataset for analysis
-- This script should be run after extracting data from each database
-- Output: merged_icu_cohort.csv

-- Step 1: Create a unified schema
DROP TABLE IF EXISTS icu_cohort_merged;
CREATE TABLE icu_cohort_merged (
    db VARCHAR(20),
    subject_id VARCHAR(50),
    hadm_id VARCHAR(50),
    age NUMERIC,
    sex VARCHAR(10),
    bmi NUMERIC,
    weight_kg NUMERIC,
    height_cm NUMERIC,
    albumin NUMERIC,
    lymph_abs NUMERIC,
    neutrophil_abs NUMERIC,
    platelet NUMERIC,
    hb NUMERIC,
    crp NUMERIC,
    nlr NUMERIC,
    gnri NUMERIC,
    charlson NUMERIC,
    hip_fracture INTEGER,
    icu_admit INTEGER,
    inhosp_death INTEGER,
    death_30d INTEGER,
    los_days NUMERIC,
    asa INTEGER,
    emop INTEGER,
    antype INTEGER,
    total_protein NUMERIC,
    height_raw NUMERIC,
    icu_los_days NUMERIC
);

-- Step 2: Insert data from each database
-- Note: Run the individual extraction scripts first and export to CSV
-- Then import using COPY command or similar

-- For MIMIC-IV:
-- COPY icu_cohort_merged FROM '/path/to/mimic_icu_export.csv' WITH CSV HEADER;

-- For INSPIRE:
-- INSERT INTO icu_cohort_merged (db, subject_id, hadm_id, age, sex, bmi, weight_kg, height_cm,
--       albumin, lymph_abs, neutrophil_abs, platelet, hb, crp, nlr, gnri,
--       charlson, hip_fracture, icu_admit, inhosp_death, death_30d, los_days,
--       asa, emop, antype, total_protein)
-- SELECT 'inspire', subject_id, hadm_id, age, sex, bmi, weight_kg, height_cm,
--       albumin, lymph_abs, neutrophil_abs, platelet, hb, crp, nlr, gnri,
--       charlson, hip_fracture, icu_admit, inhosp_death, death_30d, los_days,
--       asa, emop, antype, total_protein
-- FROM inspire_export;

-- For eICU:
-- INSERT INTO icu_cohort_merged (db, subject_id, hadm_id, age, sex, bmi, weight_kg, height_cm,
--       albumin, lymph_abs, neutrophil_abs, platelet, hb, crp, nlr, gnri,
--       charlson, hip_fracture, icu_admit, inhosp_death, death_30d, los_days,
--       height_raw, icu_los_days)
-- SELECT 'eicu', subject_id, hadm_id, age, sex, bmi, weight_kg, height_cm,
--       albumin, lymph_abs, neutrophil_abs, platelet, hb, crp, nlr, gnri,
--       charlson, hip_fracture, icu_admit, inhosp_death, death_30d, los_days,
--       height_raw, icu_los_days
-- FROM eicu_export;

-- Step 3: Data quality checks
SELECT db, 
       COUNT(*) AS total_patients,
       COUNT(*) FILTER (WHERE gnri IS NOT NULL) AS has_gnri,
       COUNT(*) FILTER (WHERE nlr IS NOT NULL) AS has_nlr,
       COUNT(*) FILTER (WHERE icu_admit = 1) AS icu_transfers,
       COUNT(*) FILTER (WHERE inhosp_death = 1) AS hospital_deaths,
       ROUND(AVG(age), 1) AS mean_age,
       ROUND(AVG(bmi), 1) AS mean_bmi
FROM icu_cohort_merged
GROUP BY db
ORDER BY db;

-- Step 4: Export to CSV for R analysis
-- COPY (SELECT * FROM icu_cohort_merged) TO '/path/to/ICU/data/icu_cohort_merged.csv' WITH CSV HEADER;