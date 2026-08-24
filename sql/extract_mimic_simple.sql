-- MIMIC-IV: simplified extraction for ICU transfer validation
WITH ortho_proc AS (
  SELECT hadm_id, MIN(chartdate) AS proc_date
  FROM mimiciv_hosp.procedures_icd
  WHERE (icd_version=10 AND (icd_code LIKE '0SR%' OR icd_code LIKE '0SQ%' OR icd_code LIKE '0SB%'
                          OR icd_code LIKE '0SH%' OR icd_code LIKE '0SP%'))
     OR (icd_version=9 AND (icd_code LIKE '815%' OR icd_code LIKE '816%' OR icd_code LIKE '79%'))
  GROUP BY hadm_id
),
base AS (
  SELECT a.subject_id, a.hadm_id, a.admittime, a.dischtime, a.hospital_expire_flag,
         pt.gender, pt.dod,
         pt.anchor_age + (EXTRACT(YEAR FROM a.admittime)::int - pt.anchor_year) AS age,
         op.proc_date,
         ROW_NUMBER() OVER (PARTITION BY a.subject_id ORDER BY a.admittime) AS rn
  FROM mimiciv_hosp.admissions a
  JOIN ortho_proc op ON op.hadm_id = a.hadm_id
  JOIN mimiciv_hosp.patients pt ON pt.subject_id = a.subject_id
),
cohort AS (
  SELECT * FROM base WHERE rn=1 AND age>=65 AND age<=110
),
-- Weight from chartevents
wt AS (
  SELECT DISTINCT ON (c.subject_id) c.subject_id, ce.valuenum AS weight_kg
  FROM cohort c
  JOIN mimiciv_icu.chartevents ce ON ce.subject_id = c.subject_id
  WHERE ce.itemid IN (762,763,7660,7693,8189,41904,42021,42198,42225,226512,226531,224639,224642,227442)
    AND ce.valuenum BETWEEN 20 AND 400
  ORDER BY c.subject_id, ce.charttime ASC
),
-- Height from chartevents
ht AS (
  SELECT DISTINCT ON (c.subject_id) c.subject_id, ce.valuenum AS height_cm
  FROM cohort c
  JOIN mimiciv_icu.chartevents ce ON ce.subject_id = c.subject_id
  WHERE ce.itemid IN (920,1395,226707,226730,459,2636,3116)
    AND ce.valuenum BETWEEN 80 AND 230
  ORDER BY c.subject_id, ce.charttime ASC
),
-- Albumin
alb AS (
  SELECT DISTINCT ON (le.hadm_id) le.hadm_id, le.valuenum AS albumin
  FROM mimiciv_hosp.labevents le
  JOIN cohort c ON c.hadm_id = le.hadm_id
  WHERE le.itemid IN (50862,53085,52022,53138)
    AND le.valuenum BETWEEN 0.5 AND 6.5
  ORDER BY le.hadm_id, le.charttime DESC
),
-- Neutrophil absolute count
neut AS (
  SELECT DISTINCT ON (le.hadm_id) le.hadm_id, le.valuenum AS neutrophil_abs
  FROM mimiciv_hosp.labevents le
  JOIN cohort c ON c.hadm_id = le.hadm_id
  WHERE le.itemid = 52075 AND le.valuenum BETWEEN 0 AND 50
  ORDER BY le.hadm_id, le.charttime DESC
),
-- Lymphocyte absolute count
lymph AS (
  SELECT DISTINCT ON (le.hadm_id) le.hadm_id, le.valuenum AS lymph_abs
  FROM mimiciv_hosp.labevents le
  JOIN cohort c ON c.hadm_id = le.hadm_id
  WHERE le.itemid IN (51133,53157) AND le.valuenum BETWEEN 0.05 AND 15
  ORDER BY le.hadm_id, le.charttime DESC
),
-- Platelet
plt AS (
  SELECT DISTINCT ON (le.hadm_id) le.hadm_id, le.valuenum AS platelet
  FROM mimiciv_hosp.labevents le
  JOIN cohort c ON c.hadm_id = le.hadm_id
  WHERE le.itemid = 51265 AND le.valuenum BETWEEN 5 AND 2000
  ORDER BY le.hadm_id, le.charttime DESC
),
-- Hemoglobin
hb AS (
  SELECT DISTINCT ON (le.hadm_id) le.hadm_id, le.valuenum AS hb
  FROM mimiciv_hosp.labevents le
  JOIN cohort c ON c.hadm_id = le.hadm_id
  WHERE le.itemid IN (50811,51222,51640) AND le.valuenum BETWEEN 2 AND 25
  ORDER BY le.hadm_id, le.charttime DESC
),
-- CRP
crp AS (
  SELECT DISTINCT ON (le.hadm_id) le.hadm_id, le.valuenum AS crp
  FROM mimiciv_hosp.labevents le
  JOIN cohort c ON c.hadm_id = le.hadm_id
  WHERE le.itemid = 50889 AND le.valuenum BETWEEN 0.1 AND 500
  ORDER BY le.hadm_id, le.charttime DESC
),
-- Hip fracture
hip AS (
  SELECT DISTINCT hadm_id FROM mimiciv_hosp.diagnoses_icd
  WHERE (icd_version=10 AND icd_code LIKE 'S72%')
     OR (icd_version=9 AND (icd_code LIKE '820%' OR icd_code LIKE '821%'))
),
-- ICU transfer
icu AS (SELECT DISTINCT hadm_id, 1 AS icu_admit FROM mimiciv_icu.icustays),
final AS (
  SELECT c.subject_id, c.hadm_id, c.age, c.gender AS sex,
         CASE WHEN wt.weight_kg > 0 AND ht.height_cm > 0
              THEN wt.weight_kg / POWER(ht.height_cm/100.0, 2) END AS bmi,
         wt.weight_kg, ht.height_cm,
         alb.albumin, lymph.lymph_abs, neut.neutrophil_abs, plt.platelet, hb.hb, crp.crp,
         -- GNRI calculation
         CASE WHEN alb.albumin BETWEEN 1.0 AND 6.0 AND ht.height_cm BETWEEN 100 AND 220
                   AND wt.weight_kg BETWEEN 20 AND 250
              THEN alb.albumin*10 + 41.7 * (wt.weight_kg / POWER(ht.height_cm/100.0, 2)) / (22 * POWER(ht.height_cm/100.0, 2))
              END AS gnri,
         -- NLR calculation
         CASE WHEN neut.neutrophil_abs IS NOT NULL AND lymph.lymph_abs IS NOT NULL
                   AND neut.neutrophil_abs BETWEEN 0 AND 50 AND lymph.lymph_abs BETWEEN 0.05 AND 15
                   AND lymph.lymph_abs > 0
              THEN neut.neutrophil_abs / lymph.lymph_abs END AS nlr,
         CASE WHEN h.hadm_id IS NOT NULL THEN 1 ELSE 0 END AS hip_fracture,
         COALESCE(i.icu_admit,0) AS icu_admit,
         c.hospital_expire_flag AS inhosp_death,
         CASE WHEN c.dod IS NOT NULL AND (c.dod - c.admittime::date) BETWEEN 0 AND 30 THEN 1 ELSE 0 END AS death_30d,
         ROUND(EXTRACT(EPOCH FROM (c.dischtime - c.admittime))/86400.0, 2) AS los_days
  FROM cohort c
  LEFT JOIN wt ON wt.subject_id=c.subject_id
  LEFT JOIN ht ON ht.subject_id=c.subject_id
  LEFT JOIN alb ON alb.hadm_id=c.hadm_id
  LEFT JOIN neut ON neut.hadm_id=c.hadm_id
  LEFT JOIN lymph ON lymph.hadm_id=c.hadm_id
  LEFT JOIN plt ON plt.hadm_id=c.hadm_id
  LEFT JOIN hb ON hb.hadm_id=c.hadm_id
  LEFT JOIN crp ON crp.hadm_id=c.hadm_id
  LEFT JOIN hip h ON h.hadm_id=c.hadm_id
  LEFT JOIN icu i ON i.hadm_id=c.hadm_id
)
SELECT 'mimiciv' AS db, subject_id, hadm_id, age, sex, bmi, weight_kg, height_cm,
       albumin, lymph_abs, neutrophil_abs, platelet, hb, crp, gnri, nlr,
       0 AS charlson, hip_fracture, icu_admit, inhosp_death, death_30d, los_days
FROM final
WHERE gnri IS NOT NULL OR nlr IS NOT NULL;