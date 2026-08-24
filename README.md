# GNRI + NLR Postoperative ICU Transfer Risk — Analysis Code & Data

Preoperative Geriatric Nutritional Risk Index (GNRI) and neutrophil-to-lymphocyte
ratio (NLR) for predicting postoperative ICU transfer (interpreted as the need for
ICU-level monitoring) in older (≥65 years) orthopedic surgical patients.

- **Development cohort:** MIMIC-IV + INSPIRE, n = 5,466 (ICU transfer n = 685, 12.5%).
- **Independent validation:** single-center matched case-control study, n = 395 (77 ICU transfers).
- **Model:** common multivariable logistic regression (GNRI, NLR, age, sex, BMI).
  Development AUC 0.712 (95% CI 0.689–0.734); fixed-coefficient validation AUC 0.726
  (95% CI 0.661–0.790).

This repository contains the analysis code and a **de-identified derived cohort**
needed to reproduce the results. A free, offline **web calculator** (EN/ZH) implements
the final model.

## Repository structure

```
ortho-gnri-nlr-icu-transfer/
├── sql/                     # Data extraction from PhysioNet (MIMIC-IV/INSPIRE/eICU)
│   ├── extract_mimic.sql
│   ├── extract_inspire.sql
│   ├── extract_eicu.sql
│   ├── merge_all_databases.sql
│   └── extract_mimic_simple.sql
├── R/                       # Analysis scripts (paths use a BASE variable; see header)
│   ├── 01_multi_factor_analysis.R     # Public-cohort Table 1, logistic, ROC, calibration, DCA
│   ├── 02_integrate_cross_db.R        # Cross-database integration & summary
│   ├── 03_local_validation.R          # Local 395-case case-control validation
│   ├── 04_redraw_figs34.R             # Manuscript figures (forest / stratified rates)
│   ├── 05_redraw_figs34_en.R
│   └── 06_supplement_figures.R
├── data/
│   └── deriv/
│       └── cohort_analysis.csv        # De-identified derived cohort (NO subject_id/hadm_id)
├── index.html               # Web calculator (English)
├── index_zh.html            # Web calculator (Chinese)
├── .nojekyll
├── manuscript/              # Submitted manuscript (V13) — our own work, not patient data
│   ├── paper_en_ICUtransfer_v13.docx
│   ├── paper_zh_ICUtransfer_v13.docx
│   └── supplementary_en_v12.docx
├── README.md
├── LICENSE
└── .gitignore
```

## Data availability & PhysioNet DUA

- The **raw data** (MIMIC-IV, INSPIRE, eICU) are accessed under PhysioNet credentialed
  access and are **not redistributed** here. Obtain them from
  [PhysioNet](https://physionet.org) after completing the required training/DUA.
- `data/deriv/cohort_analysis.csv` is a **de-identified derived cohort** prepared for
  reproducibility: direct identifiers (`subject_id`, `hadm_id`) have been removed; only
  analysis-ready variables (demographics, laboratories, computed GNRI/NLR, outcomes) are
  retained. It remains governed by the respective PhysioNet data-use agreements; users
  must also hold PhysioNet access.
- The **local 395-case validation data** is hospital patient data and is **not included**.
  The validation script (`R/03_local_validation.R`) reads a user-supplied
  `data/local_validation.csv`; request it from the authors under a data-use agreement.

## Reproducibility

- R 4.6.0 (or compatible). Required packages: `pROC`, `rms`, `showtext`, `sysfonts`
  (figure scripts also use `ggplot2`/`grid` where noted).
- Each `R/*.R` script begins with a `BASE` variable; set it to the repository root (or
  use `here::here()`). Inputs are read from `data/deriv/cohort_analysis.csv`; outputs
  (tables, figures) are written to `results/` and `figures/`.

## Prediction model (web calculator)

```
logit(p) = 1.767 − 0.093·GNRI + 0.027·NLR + 0.012·age + 0.364·male + 0.196·BMI
p = 1 / (1 + exp(−logit))
```
- `male = 1`, `female = 0`; age in years; BMI in kg/m²; GNRI and NLR per unit.
- GNRI = 1.489·albumin(g/L) + 41.7·(weight / ideal-weight), ideal weight by Lorentz.
- Web calculator (GitHub Pages): https://&lt;OWNER&gt;.github.io/ortho-gnri-nlr-icu-transfer/
  (English `index.html`; Chinese `index_zh.html`). Pure frontend, no dependencies, works offline.

## License

Code, calculator, and derived data are released under the **MIT License**. The manuscript
files are provided for reference and are governed by the journal's copyright upon publication.

## Citation

> Wang K, et al. Preoperative GNRI and NLR for predicting postoperative ICU transfer in
> older orthopedic surgical patients: a multicenter database derivation and independent
> validation. *[Journal, Year]*.
> Analysis code & data: https://doi.org/&lt;ZENODO-DOI&gt;

Archived with Zenodo: https://doi.org/&lt;ZENODO-DOI&gt; (concept) ·
https://doi.org/&lt;ZENODO-VERSION-DOI&gt; (this version).
