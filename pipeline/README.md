# Polygenic liability to bipolar disorder and plasma proteins in UK Biobank — analysis code

Analysis code for a manuscript in preparation. The code relates a precomputed bipolar
disorder polygenic risk score (BD-PRS) to Olink plasma protein levels from the UK Biobank Pharma Proteomics Project,
following the analysis plan in `config/analysis_plan.yaml`. The same code also reproduces the staged procedure of the
earlier MSc dissertation analysis, for comparison.

**This repository contains no UK Biobank data.** The data are available to approved researchers from UK Biobank
(https://www.ukbiobank.ac.uk/).
- The code runs inside an approved environment.
- It writes aggregate results only, one row per assay and model. Counts below 10, and their complements, are
  suppressed.
- `R/export.R` refuses to write a table that carries participant ids.
- These checks are a safety net. They do not replace UK Biobank's own rules for releasing results.

The `synthetic/` folder simulates data in the same layouts, so you can run and test the code without access to UK
Biobank.

## What the code does

| Step | File | Output |
| --- | --- | --- |
| Read inputs; check ids for duplicates or missing values; check that assay columns are numeric; recode variables as set in the data map | `R/io.R` | — |
| Apply exclusions in the planned order: withdrawals, ancestry, missing BD-PRS, missing core covariates, then relatedness (seeded random removal of one member of each related pair) | `R/cohort.R` | `flow.json` (Figure 1) |
| Calculate per-assay missingness; exclude assays with more than 10% missing values (an assay with exactly 10% is kept) | `R/cohort.R` | `figure_s1_protein_missingness.csv` |
| **Primary analysis.** Each assay is regressed by OLS on the BD-PRS and the core covariates. Bonferroni correction across all assays analysed; Benjamini–Hochberg as a secondary criterion; partial R² from nested models; coefficients also per SD of each assay; Spearman ρ reported descriptively | `R/models.R` | `primary.csv`, `table_s2_primary_all_assays.csv` |
| Extended Models 1–3, fitted on available cases and on a common sample, plus the restricted sample. These run on the selected assays only and are descriptive | `R/secondary.R` | `models.csv`, `table_s4_secondary_models.csv` |
| Sensitivity analyses: HC3 standard errors, centred quadratic age, and, where the data allow, eGFR, schizophrenia PRS and PRS × sex. The eGFR analysis is required if a selected assay is listed in `egfr_required_for` | `R/secondary.R` | `table_s3_sensitivity.csv` |
| A withdrawal-compliant reimplementation of the dissertation's staged procedure, following its original code (`UKB.R`). It uses the thresholds as typed there: Spearman P < 1×10⁻⁴, then 0.0021, 0.0025, 0.002631579 and 0.004545455, also dropping P = 0. Smoking and alcohol are recoded as in the original. The restricted sample is the inner join with the mental-illness file, uses the Extended Model 3 formula, and keeps participants with unknown history. Covariates aliased by `lm()` are recorded. The only intentional difference from the original is that withdrawn participants are removed | `R/legacy.R` | `table_s5_legacy_*.csv` |
| Participant characteristics and covariate missingness, with small cells suppressed | `R/export.R` | `table1_characteristics.csv`, `table_s1_covariate_missingness.csv` |
| Named values quoted in the manuscript | `R/export.R` | `result_keys.json` |
| Figures 1–3, drawn from the exported files | `figures/make_figures.py` | `Figure1_flow`, `Figure2_primary`, `Figure3_forest` (PDF, PNG) |

**Provenance.** Each run writes `PROVENANCE.txt` and `provenance.json`. They record:
- the full git commit and any tag;
- the checksums of every file the run depends on: the plan, the data map, the code lists, the manifests, the R code
  and `environment.yml`;
- the checksums of the inputs and outputs;
- the relatedness seed and the R random-number settings.
- A run is labelled `FROZEN` only when all of these hold:
  - the plan is marked FROZEN, with an ISO date and the names of the authors who agreed it;
  - every item under `decisions:` in the plan is resolved;
  - all of those files are tracked by git;
  - the working tree is clean.
- Every other run is `DRAFT` or `SYNTHETIC`, and its figures carry a watermark.
- An existing, non-empty output directory is never overwritten.

## Data maps

| File | Use |
| --- | --- |
| `config/data_map_glasgow.yaml` | The dataset the research group supplied for the dissertation (`UKB.csv`, with the column names used in the original code; proteins in columns 34–761) and the group's `anyMentalIllness` file. The dataset was already restricted to unrelated White British or White European participants. Withdrawn participants have negative ids. |
| `config/data_map.yaml` | A template for a UK Biobank RAP export (`p<field>_i<instance>_a<array>` columns, withdrawal list, kinship file, assay manifest, code lists). Lines marked `VERIFY` must be checked against the actual extract. |

## Requirements

`environment.yml` lists the versions used to develop and test the code: R 4.4.3 with data.table, sandwich, jsonlite,
yaml, digest and testthat, plus Python 3 with matplotlib.

```
conda env create -p ./env -f environment.yml
```

## Run on simulated data (no UK Biobank access needed)

```
Rscript synthetic/make_synthetic.R                    # writes synthetic/glasgow/ and synthetic/rap/
Rscript run_all.R --synthetic --layout glasgow        # writes outputs/synthetic_glasgow/
python3 figures/make_figures.py outputs/synthetic_glasgow outputs/synthetic_glasgow/figures
Rscript tests/run_tests.R                             # unit tests, plus end-to-end tests on both layouts
```

## Run on the UK Biobank data

1. Put the data files named in the data map into `data/`, which is git-ignored.
2. Resolve every item under `decisions:` in `config/analysis_plan.yaml`, and every `VERIFY` in the data map. These
   include:
   - **D-A**: the Extended Model 2 blood-pressure variables;
   - **D-B**: whether participants with an unknown psychiatric history go in the restricted sample;
   - in `data_map_glasgow.yaml`, `proteins: columns: first/last`, the names of the first and last assay columns.
     A real run stops until these are set.
3. When all authors have agreed the plan:
   - set `status: FROZEN`, `frozen_date` and `frozen_by`;
   - commit, and tag the commit (for example `plan-v1.0`).
4. Run `Rscript run_all.R --map config/data_map_glasgow.yaml --out outputs/<run_id>`, then
   `python3 figures/make_figures.py outputs/<run_id> outputs/<run_id>/figures`.
5. Compare `table_s5_legacy_counts.csv` with the counts reported in the dissertation.
   - Known reasons the counts can differ: withdrawn participants are removed here, and the dissertation's
     restricted-sample file could not be verified.
   - Find out the reason for any other difference before using results.
6. The package is an aggregate package prepared for UK Biobank's output checking. It is not proof of compliance:
   follow UK Biobank's rules before copying anything out of the secure environment.

## Citation and licence

This code is shared for review while the manuscript is in preparation. No licence has been granted yet, so all
rights are reserved; please contact the author before reusing it. A licence and citation will be added when the
paper is submitted.
