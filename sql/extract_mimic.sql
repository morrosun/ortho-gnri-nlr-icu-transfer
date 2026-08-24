-- MIMIC-IV: elderly (>=65) orthopaedic surgery cohort for ICU transfer validation
-- Primary outcome: ICU transfer (icu_admit)
-- Secondary outcomes: in-hospital mortality, 30-day mortality
-- Variables: GNRI, NLR, clinical covariates
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
hip AS (
  SELECT DISTINCT hadm_id FROM mimiciv_hosp.diagnoses_icd
  WHERE (icd_version=10 AND icd_code LIKE 'S72%')
     OR (icd_version=9 AND (icd_code LIKE '820%' OR icd_code LIKE '821%'))
),
labsel AS (
  SELECT le.hadm_id, le.charttime, le.valuenum, c.proc_date,
         CASE WHEN le.itemid IN (50862,53085,52022,53138) AND le.valuenum BETWEEN 0.5 AND 6.5 THEN 'alb'
              WHEN le.itemid IN (51133,53157) AND le.valuenum BETWEEN 0.05 AND 15 THEN 'lymph'
              WHEN le.itemid = 52075 AND le.valuenum BETWEEN 0 AND 50 THEN 'neut'
              WHEN le.itemid = 51265 AND le.valuenum BETWEEN 5 AND 2000 THEN 'platelet'
              WHEN le.itemid IN (50811,51222,51640) AND le.valuenum BETWEEN 2 AND 25 THEN 'hb'
              WHEN le.itemid = 50889 AND le.valuenum BETWEEN 0.1 AND 500 THEN 'crp' END AS analyte
  FROM mimiciv_hosp.labevents le
  JOIN cohort c ON c.hadm_id = le.hadm_id
  WHERE le.itemid IN (50862,53085,52022,53138,51133,53157,52075,51265,50811,51222,51640,50889)
    AND le.valuenum IS NOT NULL
),
preop AS (
  SELECT hadm_id, analyte, valuenum,
         ROW_NUMBER() OVER (PARTITION BY hadm_id, analyte ORDER BY charttime DESC) AS rn
  FROM labsel WHERE charttime::date <= proc_date
),
anylab AS (
  SELECT hadm_id, analyte, valuenum,
         ROW_NUMBER() OVER (PARTITION BY hadm_id, analyte ORDER BY charttime ASC) AS rn
  FROM labsel
),
lab_piv AS (
  SELECT c.hadm_id,
         COALESCE(pa.valuenum, aa.valuenum) AS albumin,
         COALESCE(pl.valuenum, al.valuenum) AS lymph_abs,
         COALESCE(pn.valuenum, an.valuenum) AS neutrophil_abs,
         COALESCE(pp.valuenum, ap.valuenum) AS platelet,
         COALESCE(ph.valuenum, ah.valuenum) AS hb,
         COALESCE(pc.valuenum, ac.valuenum) AS crp
  FROM cohort c
  LEFT JOIN preop pa ON pa.hadm_id=c.hadm_id AND pa.analyte='alb'  AND pa.rn=1
  LEFT JOIN anylab aa ON aa.hadm_id=c.hadm_id AND aa.analyte='alb' AND aa.rn=1
  LEFT JOIN preop pl ON pl.hadm_id=c.hadm_id AND pl.analyte='lymph' AND pl.rn=1
  LEFT JOIN anylab al ON al.hadm_id=c.hadm_id AND al.analyte='lymph' AND al.rn=1
  LEFT JOIN preop pn ON pn.hadm_id=c.hadm_id AND pn.analyte='neut' AND pn.rn=1
  LEFT JOIN anylab an ON an.hadm_id=c.hadm_id AND an.analyte='neut' AND an.rn=1
  LEFT JOIN preop pp ON pp.hadm_id=c.hadm_id AND pp.analyte='platelet' AND pp.rn=1
  LEFT JOIN anylab ap ON ap.hadm_id=c.hadm_id AND ap.analyte='platelet' AND ap.rn=1
  LEFT JOIN preop ph ON ph.hadm_id=c.hadm_id AND ph.analyte='hb' AND ph.rn=1
  LEFT JOIN anylab ah ON ah.hadm_id=c.hadm_id AND ah.analyte='hb' AND ah.rn=1
  LEFT JOIN preop pc ON pc.hadm_id=c.hadm_id AND pc.analyte='crp' AND pc.rn=1
  LEFT JOIN anylab ac ON ac.hadm_id=c.hadm_id AND ac.analyte='crp' AND ac.rn=1
),
bmi AS (
  SELECT subject_id, result_value::numeric AS bmi,
         ROW_NUMBER() OVER (PARTITION BY subject_id ORDER BY chartdate DESC) AS rn
  FROM mimiciv_hosp.omr
  WHERE result_name='BMI (kg/m2)' AND result_value ~ '^[0-9.]+$'
    AND result_value::numeric BETWEEN 10 AND 80
),
-- weight/height from chartevents for GNRI
wt_scan AS (
  SELECT DISTINCT ON (c.subject_id) c.subject_id, ce.valuenum AS weight_kg
  FROM cohort c
  JOIN mimiciv_icu.chartevents ce ON ce.subject_id = c.subject_id
  WHERE ce.itemid IN (762,763,7660,7693,8189,41904,42021,42198,42225,226512,226531,224639,224642,227442)
    AND ce.valuenum BETWEEN 20 AND 400
  ORDER BY c.subject_id, ce.charttime ASC
),
ht_scan AS (
  SELECT DISTINCT ON (c.subject_id) c.subject_id, ce.valuenum AS height_cm
  FROM cohort c
  JOIN mimiciv_icu.chartevents ce ON ce.subject_id = c.subject_id
  WHERE ce.itemid IN (920,1395,226707,226730,459,2636,3116)
    AND ce.valuenum BETWEEN 80 AND 230
  ORDER BY c.subject_id, ce.charttime ASC
),
wh_val AS (
  SELECT COALESCE(wt.subject_id, ht.subject_id) AS subject_id, wt.weight_kg, ht.height_cm
  FROM wt_scan wt FULL OUTER JOIN ht_scan ht ON ht.subject_id = wt.subject_id
),
icu AS (SELECT DISTINCT hadm_id, 1 AS icu_admit FROM mimiciv_icu.icustays),
final AS (
  SELECT c.subject_id, c.hadm_id, c.age, c.gender AS sex,
         b.bmi,
         wv.weight_kg, wv.height_cm,
         l.albumin, l.lymph_abs, l.neutrophil_abs, l.platelet, l.hb, l.crp,
         -- GNRI calculation
         CASE WHEN l.albumin IS NOT NULL AND wv.height_cm IS NOT NULL AND wv.weight_kg IS NOT NULL
                   AND l.albumin BETWEEN 1.0 AND 6.0 
                   AND wv.height_cm BETWEEN 100 AND 220
                   AND wv.weight_kg BETWEEN 20 AND 250
              THEN l.albumin*10 + 41.7 * (wv.weight_kg / POWER(wv.height_cm/100.0, 2)) / (22 * POWER(wv.height_cm/100.0, 2))
              END AS gnri,
         -- NLR calculation
         CASE WHEN l.neutrophil_abs IS NOT NULL AND l.lymph_abs IS NOT NULL
                   AND l.neutrophil_abs BETWEEN 0 AND 50 AND l.lymph_abs BETWEEN 0.05 AND 15
              THEN l.neutrophil_abs / l.lymph_abs END AS nlr,
         ch.charlson_comorbidity_index AS charlson,
         CASE WHEN h.hadm_id IS NOT NULL THEN 1 ELSE 0 END AS hip_fracture,
         COALESCE(i.icu_admit,0) AS icu_admit,
         c.hospital_expire_flag AS inhosp_death,
         CASE WHEN c.dod IS NOT NULL AND (c.dod - c.admittime::date) BETWEEN 0 AND 30 THEN 1 ELSE 0 END AS death_30d,
         ROUND(EXTRACT(EPOCH FROM (c.dischtime - c.admittime))/86400.0, 2) AS los_days
  FROM cohort c
  LEFT JOIN lab_piv l ON l.hadm_id=c.hadm_id
  LEFT JOIN bmi b ON b.subject_id=c.subject_id AND b.rn=1
  LEFT JOIN wh_val wv ON wv.subject_id=c.subject_id
  LEFT JOIN hip h ON h.hadm_id=c.hadm_id
  LEFT JOIN icu i ON i.hadm_id=c.hadm_id
  LEFT JOIN mimiciv_derived.charlson ch ON ch.hadm_id=c.hadm_id
)
SELECT 'mimiciv' AS db, subject_id, hadm_id, age, sex, bmi, weight_kg, height_cm,
       albumin, lymph_abs, neutrophil_abs, platelet, hb, crp, gnri, nlr,
       charlson, hip_fracture, icu_admit, inhosp_death, death_30d, los_days
FROM final
WHERE gnri IS NOT NULL OR nlr IS NOT NULL;