-- INSPIRE: elderly (>=65) orthopaedic surgery cohort for ICU transfer validation
-- Primary outcome: ICU transfer (icu_admit)
-- Secondary outcomes: in-hospital mortality, 30-day mortality
-- Variables: GNRI, NLR, clinical covariates
WITH os65 AS (
  SELECT o.*, orin_time::numeric AS orin, admission_time::numeric AS adm_t,
         discharge_time::numeric AS dis_t, opstart_time::numeric AS opstart,
         allcause_death_time::numeric AS acdt,
         ROW_NUMBER() OVER (PARTITION BY subject_id ORDER BY orin_time::numeric) AS rn
  FROM inspire.operations o
  WHERE department='OS' AND age>=65
),
cohort AS (SELECT * FROM os65 WHERE rn=1),
preop_lab AS (
  SELECT c.op_id, l.item_name, l.value::numeric AS val, l.chart_time::numeric AS ct, c.orin,
         ROW_NUMBER() OVER (PARTITION BY c.op_id, l.item_name ORDER BY l.chart_time::numeric DESC) AS rn
  FROM cohort c
  JOIN inspire.labs l ON l.subject_id = c.subject_id
  WHERE l.item_name IN ('albumin','lymphocyte','wbc','hb','crp','total_protein','seg','platelet')
    AND l.value ~ '^[0-9.]+$'
    AND l.chart_time::numeric <= c.orin
    AND l.chart_time::numeric >= c.orin - 90*1440
),
ly AS (SELECT op_id, ct AS ly_t, val AS lymph_pct FROM preop_lab WHERE item_name='lymphocyte' AND rn=1),
wb_near AS (
  SELECT p.op_id, p.val AS wbc_near, ABS(p.ct - ly.ly_t) AS dt,
         ROW_NUMBER() OVER (PARTITION BY p.op_id ORDER BY ABS(p.ct - ly.ly_t)) AS rn
  FROM preop_lab p JOIN ly ON ly.op_id=p.op_id
  WHERE p.item_name='wbc' AND ABS(p.ct - ly.ly_t) <= 1440
),
wb_any AS (SELECT op_id, val AS wbc_any FROM preop_lab WHERE item_name='wbc' AND rn=1),
lab_piv AS (
  SELECT c.op_id,
         MAX(CASE WHEN p.item_name='albumin' THEN p.val END) AS albumin,
         MAX(CASE WHEN p.item_name='hb' THEN p.val END) AS hb,
         MAX(CASE WHEN p.item_name='crp' THEN p.val END) AS crp,
         MAX(CASE WHEN p.item_name='total_protein' THEN p.val END) AS total_protein,
         MAX(CASE WHEN p.item_name='seg' THEN p.val END) AS seg_pct,
         MAX(CASE WHEN p.item_name='platelet' THEN p.val END) AS platelet
  FROM cohort c LEFT JOIN preop_lab p ON p.op_id=c.op_id AND p.rn=1
  GROUP BY c.op_id
),
hip AS (
  SELECT DISTINCT c.op_id
  FROM cohort c JOIN inspire.diagnosis d ON d.subject_id=c.subject_id
  WHERE d.icd10_cm LIKE 'S72%'
    AND ABS(d.chart_time::numeric - c.orin) <= 45*1440
),
final AS (
  SELECT c.subject_id, c.op_id, c.age, c.sex,
         CASE WHEN c.weight::numeric BETWEEN 20 AND 250 AND c.height::numeric BETWEEN 100 AND 220
              THEN c.weight::numeric END AS weight_kg,
         CASE WHEN c.height::numeric BETWEEN 100 AND 220 THEN c.height::numeric END AS height_cm,
         CASE WHEN c.weight::numeric>0 AND c.height::numeric BETWEEN 100 AND 220
              THEN c.weight::numeric / POWER(c.height::numeric/100.0, 2) END AS bmi,
         c.asa, c.emop, c.antype,
         lp.albumin,
         ly.lymph_pct,
         COALESCE(wn.wbc_near, wa.wbc_any) AS wbc,
         CASE WHEN COALESCE(wn.wbc_near, wa.wbc_any)>0 AND ly.lymph_pct BETWEEN 0 AND 100
              THEN COALESCE(wn.wbc_near, wa.wbc_any)*ly.lymph_pct/100.0 END AS lymph_abs,
         CASE WHEN COALESCE(wn.wbc_near, wa.wbc_any)>0 AND lp.seg_pct BETWEEN 0 AND 100
              THEN COALESCE(wn.wbc_near, wa.wbc_any)*lp.seg_pct/100.0 END AS neutrophil_abs,
         lp.platelet,
         lp.hb, lp.crp, lp.total_protein,
         -- NLR calculation
         CASE WHEN lp.seg_pct>0 AND ly.lymph_pct>0 THEN lp.seg_pct/ly.lymph_pct END AS nlr,
         -- GNRI calculation
         CASE WHEN lp.albumin BETWEEN 1.0 AND 6.0 AND c.height::numeric BETWEEN 100 AND 220
                   AND c.weight::numeric BETWEEN 20 AND 250
              THEN lp.albumin*10 + 41.7 * (c.weight::numeric / POWER(c.height::numeric/100.0, 2)) / (22 * POWER(c.height::numeric/100.0, 2))
              END AS gnri,
         NULL::numeric AS charlson,
         CASE WHEN h.op_id IS NOT NULL THEN 1 ELSE 0 END AS hip_fracture,
         CASE WHEN c.icuin_time IS NOT NULL AND c.icuin_time<>'' THEN 1 ELSE 0 END AS icu_admit,
         CASE WHEN c.inhosp_death_time IS NOT NULL AND c.inhosp_death_time<>'' THEN 1 ELSE 0 END AS inhosp_death,
         CASE WHEN c.acdt IS NOT NULL AND (c.acdt - c.orin) BETWEEN 0 AND 43200 THEN 1 ELSE 0 END AS death_30d,
         ROUND((c.dis_t - c.adm_t)/1440.0, 2) AS los_days
  FROM cohort c
  LEFT JOIN lab_piv lp ON lp.op_id=c.op_id
  LEFT JOIN ly ON ly.op_id=c.op_id
  LEFT JOIN wb_near wn ON wn.op_id=c.op_id AND wn.rn=1
  LEFT JOIN wb_any wa ON wa.op_id=c.op_id
  LEFT JOIN hip h ON h.op_id=c.op_id
)
SELECT 'inspire' AS db, subject_id, op_id AS hadm_id, age, sex, bmi, weight_kg, height_cm,
       asa, emop, antype, albumin, lymph_abs, neutrophil_abs, platelet, hb, crp, total_protein, nlr, gnri,
       charlson, hip_fracture, icu_admit, inhosp_death, death_30d, los_days
FROM final
WHERE gnri IS NOT NULL OR nlr IS NOT NULL;