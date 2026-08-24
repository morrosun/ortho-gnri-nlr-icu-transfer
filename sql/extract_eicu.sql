-- eICU-CRD: elderly (>=65) orthopaedic ICU cohort for ICU transfer validation
-- Primary outcome: ICU transfer (all patients are in ICU by definition)
-- Secondary outcomes: in-hospital mortality
-- Variables: GNRI, NLR, clinical covariates
-- Note: eICU patients are already in ICU, so we use ICU length of stay as outcome
WITH ortho_stay AS (
  SELECT DISTINCT patientunitstayid FROM eicu_crd.diagnosis
  WHERE diagnosisstring ILIKE 'surgery|orthopedic surgery%'
     OR diagnosisstring ILIKE '%bone fracture(s)|%lower extremity%'
     OR diagnosisstring ILIKE '%bone fracture(s)|pelvis%'
     OR diagnosisstring ILIKE '%bone fracture(s)|acetabulum%'
     OR diagnosisstring ILIKE '%dislocation|hip%'
),
base AS (
  SELECT p.patientunitstayid, p.uniquepid,
         CASE WHEN p.age = '> 89' THEN 90 ELSE p.age::numeric END AS age,
         p.gender AS sex,
         CASE WHEN p.admissionheight > 0 THEN p.admissionheight END AS admissionheight,
         CASE WHEN p.admissionweight > 0 THEN p.admissionweight END AS admissionweight,
         p.hospitaldischargestatus, p.hospitaldischargeoffset, p.unitdischargeoffset,
         p.hospitaladmitoffset,
         ROW_NUMBER() OVER (PARTITION BY p.uniquepid ORDER BY p.hospitaladmitoffset) AS rn
  FROM eicu_crd.patient p
  JOIN ortho_stay o ON o.patientunitstayid = p.patientunitstayid
),
cohort AS (SELECT * FROM base WHERE rn=1 AND age>=65),
hip AS (
  SELECT DISTINCT patientunitstayid FROM eicu_crd.diagnosis
  WHERE diagnosisstring ILIKE '%hip%' OR diagnosisstring ILIKE '%femur%'
),
labsel AS (
  SELECT l.patientunitstayid, l.labresultoffset, l.labresult,
         CASE WHEN l.labname='albumin' THEN 'alb'
              WHEN l.labname='WBC x 1000' THEN 'wbc'
              WHEN l.labname='-lymphs' THEN 'lymphpct'
              WHEN l.labname='-polys' THEN 'polypct'
              WHEN l.labname='Hgb' THEN 'hb'
              WHEN l.labname='platelets x 1000' THEN 'platelet'
              WHEN l.labname IN ('CRP','CRP-hs') THEN 'crp' END AS analyte,
         ROW_NUMBER() OVER (PARTITION BY l.patientunitstayid, l.labname ORDER BY l.labresultoffset ASC) AS rn
  FROM eicu_crd.lab l
  JOIN cohort c ON c.patientunitstayid = l.patientunitstayid
  WHERE l.labname IN ('albumin','WBC x 1000','-lymphs','-polys','Hgb','platelets x 1000','CRP','CRP-hs')
    AND l.labresultoffset BETWEEN 0 AND 1440
),
lab_piv AS (
  SELECT patientunitstayid,
         MAX(CASE WHEN analyte='alb' THEN labresult END) AS albumin,
         MAX(CASE WHEN analyte='wbc' THEN labresult END) AS wbc,
         MAX(CASE WHEN analyte='lymphpct' THEN labresult END) AS lymph_pct,
         MAX(CASE WHEN analyte='polypct' THEN labresult END) AS polys_pct,
         MAX(CASE WHEN analyte='hb' THEN labresult END) AS hb,
         MAX(CASE WHEN analyte='platelet' THEN labresult END) AS platelet,
         MAX(CASE WHEN analyte='crp' THEN labresult END) AS crp
  FROM labsel WHERE rn=1 GROUP BY patientunitstayid
),
final AS (
  SELECT c.patientunitstayid, c.uniquepid, c.age, c.sex,
         c.admissionweight AS weight_kg,
         c.admissionheight AS height_raw,
         CASE WHEN c.admissionweight>0 AND c.admissionheight>0 THEN
           c.admissionweight / POWER((CASE WHEN c.admissionheight BETWEEN 36 AND 100
                                           THEN c.admissionheight*2.54 ELSE c.admissionheight END)/100.0, 2)
         END AS bmi,
         CASE WHEN l.albumin BETWEEN 1.0 AND 6.0 THEN l.albumin END AS albumin,
         CASE WHEN l.wbc>0 AND l.lymph_pct BETWEEN 0 AND 100
              THEN l.wbc*l.lymph_pct/100.0 END AS lymph_abs,
         CASE WHEN l.wbc>0 AND l.polys_pct BETWEEN 0 AND 100
              THEN l.wbc*l.polys_pct/100.0 END AS neutrophil_abs,
         CASE WHEN l.platelet>0 AND l.platelet<2000 THEN l.platelet END AS platelet,
         l.hb,
         CASE WHEN l.crp BETWEEN 0.1 AND 500 THEN l.crp END AS crp,
         -- NLR calculation
         CASE WHEN l.wbc>0 AND l.polys_pct BETWEEN 0 AND 100 AND l.lymph_pct BETWEEN 0 AND 100
                   AND l.lymph_pct > 0
              THEN (l.wbc*l.polys_pct/100.0) / (l.wbc*l.lymph_pct/100.0)
              WHEN l.polys_pct>0 AND l.lymph_pct>0 THEN l.polys_pct/l.lymph_pct END AS nlr,
         -- GNRI calculation
         CASE WHEN l.albumin BETWEEN 1.0 AND 6.0 AND c.admissionheight>0 AND c.admissionweight>0
              THEN l.albumin*10 + 41.7 * (c.admissionweight / POWER((CASE WHEN c.admissionheight BETWEEN 36 AND 100
                                           THEN c.admissionheight*2.54 ELSE c.admissionheight END)/100.0, 2)) / (22 * POWER((CASE WHEN c.admissionheight BETWEEN 36 AND 100
                                           THEN c.admissionheight*2.54 ELSE c.admissionheight END)/100.0, 2))
              END AS gnri,
         NULL::numeric AS charlson,
         CASE WHEN h.patientunitstayid IS NOT NULL THEN 1 ELSE 0 END AS hip_fracture,
         1 AS icu_admit,
         CASE WHEN c.hospitaldischargestatus='Expired' THEN 1 ELSE 0 END AS inhosp_death,
         NULL::int AS death_30d,
         ROUND((c.hospitaldischargeoffset - c.hospitaladmitoffset)/1440.0, 2) AS los_days,
         -- ICU-specific outcomes
         ROUND(c.unitdischargeoffset/1440.0, 2) AS icu_los_days
  FROM cohort c
  LEFT JOIN lab_piv l ON l.patientunitstayid=c.patientunitstayid
  LEFT JOIN hip h ON h.patientunitstayid=c.patientunitstayid
)
SELECT 'eicu' AS db, patientunitstayid AS subject_id, patientunitstayid AS hadm_id, age, sex, bmi,
       weight_kg, height_raw, albumin, lymph_abs, neutrophil_abs, platelet, hb, crp, nlr, gnri,
       charlson, hip_fracture, icu_admit, inhosp_death, death_30d, los_days, icu_los_days
FROM final 
WHERE (gnri IS NOT NULL OR nlr IS NOT NULL)
  AND (bmi IS NULL OR bmi BETWEEN 10 AND 80);